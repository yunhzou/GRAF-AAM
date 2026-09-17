"""Coordinate/index-orientation selection inside current compressed AAM families.

This is a separate, opt-in post-processing stage. It does not assign CIP labels
or minimize RMSD. The saved action program and the selected concrete signed bond
changes remain authoritative. No atom bijections or group closures are expanded.
"""
from dataclasses import dataclass, field
from itertools import combinations
import math
import time

import numpy as np

from .alignment.index_chirality import (
    _simplex_measure, fixed_mapping_aligned_rmsd,
)


@dataclass(frozen=True)
class ChiralityConfig:
    graph_floor: float = 0.2
    orientation_tolerance: float = 0.1
    group_orientation_tolerance: float = 0.0
    mode: str = 'mutable'  # historical index orientation; 'all' is stricter
    high_coordinate: str = 'maximal'  # maximal feasible basis, or strict
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
        if self.seconds is not None and (not math.isfinite(self.seconds) or self.seconds < 0):
            raise ValueError('seconds must be nonnegative or None')


@dataclass(frozen=True)
class ChiralWitness:
    status: str  # allowed / forbidden / unknown, within the saved family union
    mapping: dict | None = None
    family_id: int | None = None
    actions: tuple = ()
    fixed_mapping_rmsd: float | None = None
    diagnostics: dict = field(default_factory=dict)


class _BudgetExpired(Exception):
    pass


class _Workspace:
    def __init__(self, decoded, candidate, config):
        import z3
        self.z3 = z3
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
        self.measures, self.mutable = {}, {}
        self.hard_rules = set()
        self.checks = 0
        self.refinements = 0
        self.mutability_queries = 0

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
            terms, total, _ = _event_model(compiled, self.decoded.index, self.deadline)
            compiled.solver.add(_same_event_pattern(terms, total, self.pattern))
            self.compiled[fid] = compiled
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
            persistent = all(mapping[a] in self.np[mapping[center]] for a in shell)
            status = ('inactive_lost_connection' if not persistent else
                      'undefined' if not sr or not sp else
                      'preserved' if sr == sp else 'violation')
            rows.append(dict(center=center, neighbors=shell, source_sign=sr, target_sign=sp, status=status))
        return rows

    def soft_violations(self, mapping, frames):
        rules = []
        for _robustness, center, shell, tol in frames:
            if not all(mapping[a] in self.np[mapping[center]] for a in shell):
                continue
            sr, _ = self.measure('r', center, shell, tol)
            sp, _ = self.measure('p', mapping[center], tuple(mapping[a] for a in shell), tol)
            if sr and sp and sr != sp:
                rules.append(self.rule(center, shell, mapping, tol, False))
        return rules

    def solve(self, frames, preferred=None):
        soft_rules = set()
        if preferred is not None:
            mapping = dict(preferred['mapping'])
            bad, _ = self.ordinary(mapping)
            self.refinements += len(set(bad) - self.hard_rules)
            self.hard_rules.update(bad)
            soft_rules.update(self.soft_violations(mapping, frames))
            self.refinements += len(soft_rules)
            if not bad and not soft_rules:
                return preferred
        for fid in self.family_order:
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
        return None


def select_chiral_witness(decoded, candidate, config=None):
    """Select a coordinate-consistent witness for one decoded candidate.

    Preserve its *concrete* broken/formed edges, not only their symmetry class.
    Query all saved families, including families whose event decoding timed out.
    Ordinary persistent three/four-ligand orientations are hard when the saved
    relation permits a setwise ligand shuffle (``mode='mutable'``). Fixed
    orientation changes are reported, not rejected: genuine reaction inversions
    are possible. ``mode='all'`` enforces every defined ordinary orientation.

    Higher-coordinate frames use a robustness-ordered maximal feasible basis,
    with every excluded frame reported; ``high_coordinate='strict'`` makes all
    frames hard. This is a maximal basis, not a maximum-cardinality claim.
    Undefined/planar orientations and lost ligand connections are not constrained.

    No solution/count cap is added. ``seconds`` is an optional soft watchdog;
    use process isolation for a hard wall limit. Unknown never means forbidden.
    The returned RMSD is a proper fit of the selected fixed mapping, not a
    proof of the minimum RMSD over all chirality-valid mappings.
    """
    decoded._check_candidate(candidate)
    config = config or ChiralityConfig()
    w = _Workspace(decoded, candidate, config)
    accepted, excluded = [], []
    diagnostics = dict(policy='saved_family_index_orientation', event_scope='concrete_signed_edges',
                       mode=config.mode, high_coordinate=config.high_coordinate,
                       rmsd_optimized=False, graph_floor=config.graph_floor,
                       orientation_tolerance=config.orientation_tolerance,
                       group_orientation_tolerance=config.group_orientation_tolerance,
                       candidate_id=candidate.id)
    try:
        w.budget()
        preferred = dict(mapping=sorted(candidate.mapping.items()), family_id=None, actions=())
        result = w.solve([], preferred)
        if result is not None:
            frames = w.soft_frames(candidate.mapping)
            if config.high_coordinate == 'strict':
                result = w.solve(frames, result)
                accepted = frames if result is not None else []
            else:
                for frame in frames:
                    w.budget()
                    trial = w.solve([*accepted, frame], result)
                    if trial is None:
                        excluded.append(frame)
                    else:
                        accepted.append(frame)
                        result = trial
        if result is None:
            status, mapping, rmsd = 'forbidden', None, None
            diagnostics['reason'] = 'no saved mapping satisfies the requested orientation constraints'
        else:
            status, mapping = 'allowed', dict(result['mapping'])
            _, ordinary = w.ordinary(mapping)
            diagnostics['ordinary_frames'] = ordinary
            diagnostics['high_coordinate_frames'] = w.frame_diagnostics(mapping, accepted)
            assert w.decoded.index.describe([mapping[i] for i in range(w.decoded.index.n)])['events'] == w.pattern['events']
            assert not w.ordinary(mapping)[0] and not w.soft_violations(mapping, accepted)
            rmsd = fixed_mapping_aligned_rmsd(mapping, w.r.coordinates, w.p.coordinates)
    except (_BudgetExpired, TimeoutError) as exc:
        status, mapping, rmsd, result = 'unknown', None, None, None
        diagnostics['reason'] = str(exc)
    diagnostics.update(seconds=time.perf_counter() - w.start, solver_checks=w.checks,
                       compiled_families=len(w.compiled), local_refinements=w.refinements,
                       local_geometry_evaluations=len(w.measures), mutability_queries=w.mutability_queries,
                       retained_high_coordinate_frames=accepted, reconfigured_high_coordinate_frames=excluded,
                       atom_bijections_enumerated=0)
    return ChiralWitness(status, mapping, None if result is None else result.get('family_id'),
                         () if result is None else tuple(result.get('actions', ())), rmsd, diagnostics)
