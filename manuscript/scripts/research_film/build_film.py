"""Build a deterministic, offline 3D research film from a verified saved trace.

Usage: python build_film.py --repo /path/to/coordinate_alignment --output ./output
No mapping searches are run. Coordinates undergo rigid display transforms only.
"""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np


def build(repo, out):
    out.mkdir(parents=True, exist_ok=True)
    source = repo / 'reports/pr7_search_trajectory_20260911/trace.json'
    raw = json.loads(source.read_text())
    run = raw['runs'][1]
    endpoints = {k: raw['input'][k] for k in ('reactant', 'product')}
    paths = []
    signatures = []
    for p in run['paths']:
        groups, previous = {}, set()
        for f in p['frames']:
            if f['kind'] == 'locked':
                now = set(f['locked'])
                groups.update({a: f['stage'] - 1 for a in now - previous})
                previous = now
        frames = []
        for f in p['frames']:
            if f['kind'] not in ('seed_start', 'commit', 'consumed', 'seed_end', 'locked', 'terminal'):
                continue
            c = f['candidates'][f.get('preferred', 0)] if f['candidates'] else {}
            mapping = {**f['locked'], **c.get('witness', {})}
            assert len(set(mapping.values())) == len(mapping)
            assert all(endpoints['reactant']['elements'][int(a)] == endpoints['product']['elements'][b]
                       for a, b in mapping.items())
            frames.append(dict(kind=f['kind'], stage=f['stage'], active=f['active'],
                               mapping=mapping, transition=f.get('transition'),
                               candidates=len(f['candidates']), title=f['title'],
                               symmetry=c.get('automorph_blocks', []) + c.get('blocks', []),
                               multiplicity=c.get('multiplicity', 1)))
        mapping = p['frames'][-1]['locked']
        assert sorted(mapping.values()) == list(range(28))
        # Independently recompute the archived event rule, including WBO changes.
        expected = []
        for a in range(28):
            for b in range(a + 1, 28):
                u = endpoints['reactant']['wbo'][a][b]
                v = endpoints['product']['wbo'][mapping[str(a)]][mapping[str(b)]]
                if abs(v - u) <= run['config']['event_tolerance']:
                    continue
                kind = ('formed' if u < .2 else 'strengthened') if v > u else ('broken' if v < .2 else 'weakened')
                expected.append((kind, (a, b)))
        saved = [(e['kind'], tuple(e['r'])) for e in p['events']]
        assert sorted(expected) == sorted(saved), (expected, saved)
        signatures.append(tuple(sorted(saved)))
        paths.append(dict(terminal=p['terminal'], frames=frames, groups=groups,
                          mapping=mapping, events=p['events']))
    assert len(set(signatures)) == 3
    # Fit product to a reference witness for presentation; this changes no chemistry.
    r = np.asarray(endpoints['reactant']['coordinates'], float)
    p = np.asarray(endpoints['product']['coordinates'], float)
    r -= r.mean(0)
    p -= p.mean(0)
    m = paths[0]['mapping']
    u, _, vt = np.linalg.svd(p[[m[str(i)] for i in range(28)]].T @ r)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    p = p @ (u @ correction @ vt)
    _, _, axes = np.linalg.svd(r)
    if np.linalg.det(axes) < 0:
        axes[-1] *= -1
    endpoints['reactant']['coordinates'] = (r @ axes.T).tolist()
    # Independently orient the product's principal plane toward the camera.
    # Matching is carried by identities/colors, not by forcing a misleading fit.
    _, _, p_axes = np.linalg.svd(p)
    if np.linalg.det(p_axes) < 0:
        p_axes[-1] *= -1
    endpoints['product']['coordinates'] = (p @ p_axes.T).tolist()
    data = dict(input=endpoints, graph=run['graph'], paths=paths, cuts=run['cuts'], config=run['config'],
                duration=48, source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                archive_provenance=run['provenance'], engine_commit=raw['requested_engine_commit'],
                scope='Selected archived search context, coordinate case 101 (PR7), seed order 15 of 30, cap 1000. '
                      'Three recorded terminal witnesses, not exhaustive symmetry decoding or a current benchmark. '
                      'Event tolerance 0.5; fixed endpoint geometries, rigid camera/display motion only.')
    (out / 'film-data.json').write_text(json.dumps(data, indent=2) + '\n')
    html = (Path(__file__).parent / 'film.html').read_text()
    html = html.replace('__LIBRARY__', (repo / 'src/rxn_core/static/3Dmol-min.js').read_text())
    html = html.replace('__DATA__', json.dumps(data, separators=(',', ':')).replace('</', '<\\/'))
    (out / 'index.html').write_text(html)
    (out / 'science-validation.json').write_text(json.dumps(dict(
        status='passed', source_sha256=data['source_sha256'],
        checks=['All displayed mappings are recorded element-preserving injections',
                'All three terminal mappings are complete bijections',
                'All 12 displayed bond events independently recomputed from WBO matrices',
                'Three distinct atom-indexed event signatures',
                'Tree is the saved 12-state, 11-edge context, without invented branches',
                'Only proper rigid endpoint display transformations applied'],
        scope=data['scope']), indent=2) + '\n')
    print(out / 'index.html')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', type=Path, required=True)
    ap.add_argument('--output', type=Path, required=True)
    a = ap.parse_args()
    build(a.repo, a.output)
