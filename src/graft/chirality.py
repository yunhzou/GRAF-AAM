"""Coordinate/index-orientation selection inside current compressed AAM families.

This is a separate, opt-in post-processing stage. It does not assign CIP labels
or minimize RMSD. The saved action program and the selected concrete signed bond
changes remain authoritative. No atom bijections or group closures are expanded.
"""
from dataclasses import dataclass, field
from itertools import combinations
from functools import lru_cache
import math
import time

import numpy as np

from .alignment.index_chirality import (
    _simplex_measure,
)


@dataclass(frozen=True)
class ChiralityConfig:
    graph_floor: float = 0.2
    orientation_tolerance: float = 0.1
    group_orientation_tolerance: float = 0.0
    mode: str = 'all'  # 'mutable' is an explicit historical relaxation
    high_coordinate: str = 'strict'  # 'maximal' is an explicit heuristic relaxation
    high_coordinate_scope: str = 'union'  # only affects explicit maximal relaxation
    seconds: float | None = None

    def __post_init__(self):
        for name in ('graph_floor', 'orientation_tolerance', 'group_orientation_tolerance'):
            value = getattr(self, name)
            if not math.isfinite(value) or value < 0 or (name == 'graph_floor' and value == 0):
                raise ValueError(f'invalid {name}')
        if self.mode not in ('mutable', 'all'):
            raise ValueError("mode must be 'mutable' or 'all'")
        if self.high_coordinate not in ('maximal', 'strict'):
            raise ValueError("high_coordinate must be 'maximal' or 'strict'")
        if self.high_coordinate_scope not in ('selected_family', 'union'):
            raise ValueError("high_coordinate_scope must be 'selected_family' or 'union'")
        if self.seconds is not None and (not math.isfinite(self.seconds) or self.seconds < 0):
            raise ValueError('seconds must be nonnegative or None')


@dataclass(frozen=True)
class ChiralWitness:
    status: str  # allowed / relaxed / forbidden / unknown; saved family union only
    mapping: dict | None = None
    family_id: int | None = None
    actions: tuple = ()
    diagnostics: dict = field(default_factory=dict)


class _OutsideSavedRelation(ValueError):
    pass


class _BudgetExpired(Exception):
    pass


@lru_cache(maxsize=2048)
def _action_orbits(action):
    """Single-atom reachable sets; these are rejection bounds, never witnesses."""
    kind, data = action
    if kind == 'pool':
        group = frozenset(data)
        return {a: group for a in data}
    # Connected components of generator edges are exactly single-atom orbits.
    # This does not enumerate group elements or permit independent orbit swaps.
    neighbors = {}
    for generator in data:
        for a, b in enumerate(generator):
            if a != b:
                neighbors.setdefault(a, set()).add(b)
                neighbors.setdefault(b, set()).add(a)
    result = {}
    for atom in neighbors:
        if atom in result:
            continue
        todo, component = [atom], {atom}
        while todo:
            for other in neighbors[todo.pop()]:
                if other not in component:
                    component.add(other)
                    todo.append(other)
        orbit = frozenset(component)
        result.update((a, orbit) for a in component)
    return result


@lru_cache(maxsize=512)
def _program_domains(actions, degree):
    """Reachable single-atom images for an ordered product, ignoring coupling."""
    images = [frozenset((p,)) for p in range(degree)]
    for action in actions:
        table = _action_orbits(action)
        cache = {}
        for a, domain in enumerate(images):
            if domain not in cache:
                cache[domain] = frozenset().union(*(table.get(p, (p,)) for p in domain))
            images[a] = cache[domain]
    return tuple(images)


class _Workspace:
    def __init__(self, decoded, candidate, config, assignments=None):
        import z3
        self.z3 = z3
        self.assignments = dict(assignments or {})
        self.decoded, self.candidate, self.config = decoded, candidate, config
        self.problem = decoded.catalogue.problem
        self.start = time.perf_counter()
        self.deadline = math.inf if config.seconds is None else self.start + config.seconds
        self.r, self.p = self.problem.reactant, self.problem.product
        self.nr = [tuple(map(int, np.flatnonzero(row >= config.graph_floor))) for row in self.r.wbo]
        self.np = [frozenset(map(int, np.flatnonzero(row >= config.graph_floor))) for row in self.p.wbo]
        # Concrete event equality keeps the chosen reaction core fixed. This is
        # intentionally stronger than the decoder's source-symmetry event class.
        self.pattern = decoded.index.describe([candidate.mapping[i] for i in range(decoded.index.n)])
        self.family_order = list(dict.fromkeys([*candidate.family_ids, *range(len(decoded.catalogue.families))]))
        self.compiled = {}
        self.compilations = 0
        self.forced_rejections = 0
        self.first_forced_conflict = None
        self.bound_rejections = 0
        self.measures, self.mutable = {}, {}
        self.hard_rules = set()
        self.checks = 0
        self.refinements = 0
        self.mutability_queries = 0
        self.families_checked = set()

    def budget(self):
        if time.perf_counter() >= self.deadline:
            raise _BudgetExpired('time budget')

    def measure(self, side, center, neighbors, tol):
        self.budget()
        canonical = tuple(sorted(neighbors))
        key = side, center if len(neighbors) == 3 else None, canonical, tol
        if key not in self.measures:
            endpoint = self.r if side == 'r' else self.p
            self.measures[key] = _simplex_measure(endpoint.coordinates, center, canonical, tol)
        measure = self.measures[key]
        odd = sum(a > b for i, a in enumerate(neighbors) for b in neighbors[i + 1:]) % 2
        return measure.sign * (-1 if odd else 1), abs(measure.normalized)

    def compile(self, fid):
        self.budget()
        if fid not in self.compiled:
            from .family_query import compile_path
            from .event_patterns import _event_model, _same_event_pattern
            family = self.decoded.catalogue.families[fid]
            compiled = compile_path(family.as_path(self.problem), self.problem, {},
                                    source_atoms=(), complete_reference=False)
            compiled.solver.add(*(compiled.values[a][0] == p for a, p in self.assignments.items()))
            terms, total, _ = _event_model(compiled, self.decoded.index, self.deadline)
            compiled.solver.add(_same_event_pattern(terms, total, self.pattern))
            self.compiled[fid] = compiled
            self.compilations += 1
        return self.compiled[fid]

    def check(self, compiled, constraints):
        self.budget()
        solver = compiled.solver
        solver.push()
        try:
            solver.add(*constraints)
            remaining = self.deadline - time.perf_counter()
            if remaining <= 0:
                raise _BudgetExpired('encoding budget')
            solver.set(timeout=0 if math.isinf(remaining) else max(1, int(remaining * 1000)))
            self.checks += 1
            status = solver.check()
            if status == self.z3.unknown:
                raise _BudgetExpired(solver.reason_unknown())
            return compiled.realize(solver.model()) if status == self.z3.sat else None
        finally:
            solver.pop()

    def shell(self, mapping, center):
        return tuple(a for a in self.nr[center] if mapping[a] in self.np[mapping[center]])

    def scope(self, values, center, neighbors, pc, targets, exact_shell):
        z3 = self.z3
        conditions = [values[center][0] == pc]
        conditions += [z3.Or(*(values[a][0] == p for p in targets)) for a in neighbors]
        if exact_shell:
            # Necessary when 4/5-coordinate frames cross an active-edge floor.
            conditions += [z3.Not(z3.Or(*(values[a][0] == p for p in self.np[pc])))
                           for a in self.nr[center] if a not in neighbors]
        return z3.And(*conditions)

    def is_mutable(self, mapping, center, neighbors):
        if self.config.mode == 'all':
            return True
        pc = mapping[center]
        targets = tuple(sorted(mapping[a] for a in neighbors))
        images = tuple(mapping[a] for a in neighbors)
        key = center, neighbors, pc, targets, images
        if key in self.mutable:
            return self.mutable[key]
        self.mutability_queries += 1
        for fid in self.family_order:
            self.budget()
            family = self.decoded.catalogue.families[fid]
            reachable = _program_domains(family.actions, len(family.mapping))
            source = dict(family.mapping)
            domains = {a: reachable[source[a]] for a in (center, *self.nr[center])}
            target_set = frozenset(targets)
            if (pc not in domains[center]
                    or any(not (domains[a] & target_set) for a in neighbors)
                    or any(not (domains[a] - self.np[pc]) for a in self.nr[center] if a not in neighbors)
                    or all((domains[a] & target_set) <= {mapping[a]} for a in neighbors)):
                # A necessary condition failed. Bounds can only reject;
                # every positive answer still comes from the full solver.
                self.bound_rejections += 1
                continue
            c = self.compile(fid)
            scope = self.scope(c.values, center, neighbors, pc, targets, True)
            changed = self.z3.Or(*(c.values[a][0] != mapping[a] for a in neighbors))
            if self.check(c, [scope, changed]) is not None:
                self.mutable[key] = True
                return True
        self.mutable[key] = False
        return False

    def rule(self, center, neighbors, mapping, tol, exact_shell):
        # One learned local rule covers all 6/24 ligand orderings at once.
        # Target-set geometry is cached once; permutation parity gives its sign.
        pc = mapping[center]
        targets = tuple(sorted(mapping[a] for a in neighbors))
        sr, _ = self.measure('r', center, neighbors, tol)
        sp, _ = self.measure('p', pc, targets, tol)
        assert sr and sp
        return center, neighbors, pc, targets, sr != sp, exact_shell

    def expression(self, compiled, rule):
        center, neighbors, pc, targets, odd, exact = rule
        values = compiled.values
        parity = self.z3.BoolVal(False)
        for i, a in enumerate(neighbors):
            for b in neighbors[i + 1:]:
                parity = self.z3.Xor(parity, values[a][0] > values[b][0])
        return self.z3.Implies(self.scope(values, center, neighbors, pc, targets, exact), parity == odd)

    def ordinary(self, mapping):
        rules, rows = [], []
        for center in range(len(self.nr)):
            neighbors = self.shell(mapping, center)
            if len(neighbors) not in (3, 4):
                continue
            tol = self.config.orientation_tolerance
            sr, _ = self.measure('r', center, neighbors, tol)
            sp, _ = self.measure('p', mapping[center], tuple(mapping[a] for a in neighbors), tol)
            state = 'preserved' if sr == sp and sr else 'undefined'
            if sr and sp and sr != sp:
                if self.is_mutable(mapping, center, neighbors):
                    state = 'violation'
                    rules.append(self.rule(center, neighbors, mapping, tol, True))
                else:
                    state = 'fixed_orientation_change'
            rows.append(dict(center=center, neighbors=neighbors, source_sign=sr, target_sign=sp, status=state))
        return rules, rows

    def soft_frames(self, reference):
        frames = []
        for center, neighbors in enumerate(self.nr):
            # Include centers high-coordinate at either endpoint. Frames whose
            # ligands do not persist under a particular mapping are inactive.
            high_target = any(self.r.elements[center] == self.p.elements[p] and len(ns) > 4
                              for p, ns in enumerate(self.np))
            if len(neighbors) <= 4 and not high_target:
                continue
            sizes = (3, 4) if len(neighbors) > 4 else (3,)
            for size in sizes:
                tol = self.config.group_orientation_tolerance if size == 3 else self.config.orientation_tolerance
                for shell in combinations(neighbors, size):
                    sr, ar = self.measure('r', center, shell, tol)
                    if not sr:
                        continue
                    _, ap = self.measure('p', reference[center], tuple(reference[a] for a in shell), tol)
                    frames.append((min(ar, ap), center, shell, tol))
        return sorted(frames, key=lambda x: (-x[0], x[1:]))

    def frame_diagnostics(self, mapping, frames):
        rows = []
        for _robustness, center, shell, tol in frames:
            sr, _ = self.measure('r', center, shell, tol)
            sp, _ = self.measure('p', mapping[center], tuple(mapping[a] for a in shell), tol)
            persistent = (self.high_active(mapping, center) and
                          all(mapping[a] in self.np[mapping[center]] for a in shell))
            status = ('inactive' if not persistent else
                      'undefined' if not sr or not sp else
                      'preserved' if sr == sp else 'violation')
            rows.append(dict(center=center, neighbors=shell, source_sign=sr, target_sign=sp, status=status))
        return rows

    def high_active(self, mapping, center):
        # Potential targets only control which frames are generated. The actual
        # mapped center must be high-coordinate at at least one endpoint.
        return len(self.nr[center]) > 4 or len(self.np[mapping[center]]) > 4

    def soft_violations(self, mapping, frames):
        rules = []
        for _robustness, center, shell, tol in frames:
            if not self.high_active(mapping, center) or not all(mapping[a] in self.np[mapping[center]] for a in shell):
                continue
            sr, _ = self.measure('r', center, shell, tol)
            sp, _ = self.measure('p', mapping[center], tuple(mapping[a] for a in shell), tol)
            if sr and sp and sr != sp:
                rules.append(self.rule(center, shell, mapping, tol, False))
        return rules

    def locate(self, result):
        """Certify which saved family contains a geometrically valid witness."""
        if result.get('family_id') is not None:
            return result
        for fid in self.family_order:
            if not self.can_contain(fid, dict(result['mapping'])):
                continue
            compiled = self.compile(fid)
            realized = self.check(compiled, [compiled.values[a][0] == p
                                           for a, p in dict(result['mapping']).items()])
            if realized is not None:
                return dict(realized, family_id=fid)
            self.compiled.pop(fid, None)
        raise _OutsideSavedRelation('candidate witness is not contained in its saved AAM families')

    def can_contain(self, fid, assignments):
        family = self.decoded.catalogue.families[fid]
        reachable = _program_domains(family.actions, len(family.mapping))
        source = dict(family.mapping)
        return all(a in source and p in reachable[source[a]] for a, p in assignments.items())

    def forced_conflict(self, fid, frames):
        """Reject only when singleton reachable images force an orientation flip.

        The domains overapproximate the ordered correlated program. A singleton
        is therefore fixed in every realization. This can prove infeasibility,
        never feasibility, and does not enumerate mappings or group elements.
        """
        family = self.decoded.catalogue.families[fid]
        reachable = _program_domains(family.actions, len(family.mapping))
        domains = {a: reachable[p] for a, p in family.mapping}
        fixed = {a: next(iter(images)) for a, images in domains.items() if len(images) == 1}
        fixed.update(self.assignments)
        def conflict(kind, center, shell, tol):
            sr, _ = self.measure('r', center, shell, tol)
            sp, _ = self.measure('p', fixed[center], tuple(fixed[a] for a in shell), tol)
            return (dict(kind=kind, family_id=fid, center=center, neighbors=shell,
                         target_center=fixed[center], target_neighbors=tuple(fixed[a] for a in shell),
                         source_sign=sr, target_sign=sp) if sr and sp and sr != sp else None)
        if self.config.mode == 'all':
            for center, neighbors in enumerate(self.nr):
                if center not in fixed or any(a not in fixed for a in neighbors):
                    continue
                shell = tuple(a for a in neighbors if fixed[a] in self.np[fixed[center]])
                if len(shell) in (3, 4):
                    bad = conflict('ordinary', center, shell, self.config.orientation_tolerance)
                    if bad:
                        return bad
        for _robustness, center, shell, tol in frames:
            if center not in fixed or any(a not in fixed for a in shell):
                continue
            if not self.high_active(fixed, center) or not all(fixed[a] in self.np[fixed[center]] for a in shell):
                continue
            bad = conflict('high_coordinate', center, shell, tol)
            if bad:
                return bad
        return None

    def solve(self, frames, preferred=None, family_ids=None):
        soft_rules = set()
        if preferred is not None and all(dict(preferred['mapping']).get(a) == p for a, p in self.assignments.items()):
            mapping = dict(preferred['mapping'])
            bad, _ = self.ordinary(mapping)
            self.refinements += len(set(bad) - self.hard_rules)
            self.hard_rules.update(bad)
            soft_rules.update(self.soft_violations(mapping, frames))
            self.refinements += len(soft_rules)
            if not bad and not soft_rules:
                return preferred
        for fid in self.family_order if family_ids is None else family_ids:
            self.budget()
            self.families_checked.add(fid)
            if not self.can_contain(fid, self.assignments):
                self.bound_rejections += 1
                continue
            conflict = self.forced_conflict(fid, frames)
            if conflict is not None:
                self.forced_rejections += 1
                if self.first_forced_conflict is None:
                    self.first_forced_conflict = conflict
                self.compiled.pop(fid, None)
                continue
            compiled = self.compile(fid)
            while True:
                rules = self.hard_rules | soft_rules
                result = self.check(compiled, [self.expression(compiled, r) for r in rules])
                if result is None:
                    break
                mapping = dict(result['mapping'])
                bad, _ = self.ordinary(mapping)
                soft = self.soft_violations(mapping, frames)
                if not bad and not soft:
                    return dict(result, family_id=fid)
                new = (set(bad) | set(soft)) - rules
                assert new, 'local orientation refinement made no progress'
                self.refinements += len(new)
                self.hard_rules.update(bad)
                soft_rules.update(soft)
            self.compiled.pop(fid, None)  # release an exhausted family model
        return None


def _select(decoded, candidate, config=None, *, assignments=None, preferred_mapping=None):
    """Shared selector and fixed-assignment query implementation."""
    decoded._check_candidate(candidate)
    config = config or ChiralityConfig()
    w = _Workspace(decoded, candidate, config, assignments)
    accepted, excluded = [], []
    diagnostics = dict(policy='saved_family_index_orientation', event_scope='concrete_signed_edges',
                       mode=config.mode, high_coordinate=config.high_coordinate,
                       high_coordinate_scope=config.high_coordinate_scope,
                       graph_floor=config.graph_floor,
                       orientation_tolerance=config.orientation_tolerance,
                       group_orientation_tolerance=config.group_orientation_tolerance,
                       candidate_id=candidate.id,
                       guarantee='coordinate_orientation_not_CIP_or_clash_free',
                       search_scope='all_saved_families',
                       search_exhaustive=False)
    try:
        w.budget()
        proposed = dict(candidate.mapping if preferred_mapping is None else preferred_mapping)
        # Even an already-valid preferred witness must have a certificate in the
        # saved relation. Checking geometry alone could accept an external map.
        preferred = None
        if all(proposed.get(a) == p for a, p in w.assignments.items()):
            try:
                preferred = w.locate(dict(mapping=sorted(proposed.items()), family_id=None, actions=()))
            except _OutsideSavedRelation:
                pass
        frames = w.soft_frames(candidate.mapping)
        if config.high_coordinate == 'strict':
            accepted = frames
            result = w.solve(frames, preferred)
        else:
            result = w.solve([], preferred)
        if result is not None and config.high_coordinate != 'strict':
            family_ids = None
            if frames and config.high_coordinate_scope == 'selected_family':
                result = w.locate(result)
                family_ids = (result['family_id'],)
                diagnostics['high_coordinate_family_id'] = result['family_id']
            for frame in frames:
                w.budget()
                trial = w.solve([*accepted, frame], result, family_ids=family_ids)
                if trial is None:
                    excluded.append(frame)
                else:
                    accepted.append(frame)
                    result = trial
        if result is None:
            status, mapping = 'forbidden', None
            diagnostics['reason'] = 'no saved mapping satisfies the fixed assignments, concrete events, and requested orientation constraints'
            diagnostics['search_exhaustive'] = True
        else:
            status, mapping = 'allowed', dict(result['mapping'])
            _, ordinary = w.ordinary(mapping)
            diagnostics['ordinary_frames'] = ordinary
            diagnostics['high_coordinate_frames'] = w.frame_diagnostics(mapping, accepted)
            diagnostics['all_high_coordinate_frames'] = w.frame_diagnostics(mapping, frames)
            full_violations = [f for f in ordinary if f['status'] in ('violation', 'fixed_orientation_change')]
            full_violations += [f for f in diagnostics['all_high_coordinate_frames'] if f['status'] == 'violation']
            diagnostics['orientation_violations'] = full_violations
            diagnostics['strict_orientation_satisfied'] = not full_violations
            if full_violations:
                status = 'relaxed'
            assert w.decoded.index.describe([mapping[i] for i in range(w.decoded.index.n)])['events'] == w.pattern['events']
            assert not w.ordinary(mapping)[0] and not w.soft_violations(mapping, accepted)
    except (_BudgetExpired, TimeoutError) as exc:
        status, mapping, result = 'unknown', None, None
        diagnostics['reason'] = str(exc)
    if status == 'forbidden':
        diagnostics['strict_orientation_satisfied'] = False
    if config.high_coordinate != 'strict':
        diagnostics['search_scope'] = 'all_saved_families_then_' + config.high_coordinate_scope + '_relaxation'
    diagnostics['families_checked'] = len(w.families_checked)
    diagnostics['saved_family_count'] = len(decoded.catalogue.families)
    diagnostics.update(seconds=time.perf_counter() - w.start, solver_checks=w.checks,
                       compiled_families=w.compilations,
                       resident_family_models=len(w.compiled),
                       forced_orientation_rejections=w.forced_rejections,
                       first_forced_conflict=w.first_forced_conflict, local_refinements=w.refinements,
                       mutability_bound_rejections=w.bound_rejections,
                       local_geometry_evaluations=len(w.measures), mutability_queries=w.mutability_queries,
                       retained_high_coordinate_frames=accepted, reconfigured_high_coordinate_frames=excluded,
                       atom_bijections_enumerated=0)
    return ChiralWitness(status=status, mapping=mapping,
                         family_id=None if result is None else result.get('family_id'),
                         actions=() if result is None else tuple(result.get('actions', ())),
                         diagnostics=diagnostics)


def _validate_assignments(decoded, assignments, *, complete=False):
    fixed = dict(assignments)
    n = decoded.index.n
    if complete and set(fixed) != set(range(n)):
        raise ValueError('preferred_mapping must assign every source atom')
    if any(not isinstance(a, int) or not isinstance(p, int)
           or a not in range(n) or p not in range(n) for a, p in fixed.items()):
        raise ValueError('assignments require valid integer endpoint indices')
    problem = decoded.catalogue.problem
    if len(set(fixed.values())) != len(fixed) or any(
            problem.reactant.elements[a] != problem.product.elements[p] for a, p in fixed.items()):
        raise ValueError('assignments must be injective and element compatible')
    return fixed


def select_chiral_witness(decoded, candidate, config=None, *, preferred_mapping=None):
    """Return a certified coordinate-orientation-preserving AAM witness.

    Defaults require every defined ordinary persistent-shell orientation and
    every defined persistent high-coordinate simplex to be preserved. Search
    the whole saved family union, fixing concrete broken/formed edges. No frame
    is silently dropped. ``forbidden`` proves infeasibility within that saved
    relation; ``unknown`` is inconclusive. This is not a CIP/E-Z assignment,
    a guarantee of a physical pathway, or an interpolation-clash optimizer.

    ``preferred_mapping`` is an optional complete witness to retain if it is
    both in the saved relation and feasible; otherwise search normally. With
    no preference the decoded witness is tried first. Choosing one witness
    does not remove other feasible families: use ``query_chiral_witness``.

    Explicit historical mutable/maximal policies can return ``relaxed`` when
    the selected mapping violates full orientation preservation. Such results
    are never labeled ``allowed`` merely because the violated frames were
    omitted. Maximal-basis selection is heuristic and not a complete filter.
    """
    decoded._check_candidate(candidate)
    if preferred_mapping is not None:
        preferred_mapping = _validate_assignments(decoded, preferred_mapping, complete=True)
    return _select(decoded, candidate, config, preferred_mapping=preferred_mapping)


def query_chiral_witness(decoded, candidate, assignments, config=None):
    """Does a strict post-chirality witness extend these R->P assignments?

    Supply the entire mapping for exact membership; unspecified atoms may move
    jointly. Uses the same orientation predicate as default selection, with no
    selected-family restriction or witness-dependent frame basis. Partial
    assignments and saved anchors must coexist, and concrete events stay fixed.
    A negative answer covers only saved complete families, not all chemistry.
    """
    decoded._check_candidate(candidate)
    config = config or ChiralityConfig()
    if config.mode != 'all' or config.high_coordinate != 'strict':
        raise ValueError('subset queries require mode="all", high_coordinate="strict"')
    fixed = _validate_assignments(decoded, assignments)
    preferred = fixed if len(fixed) == decoded.index.n else None
    return _select(decoded, candidate, config, assignments=fixed, preferred_mapping=preferred)
