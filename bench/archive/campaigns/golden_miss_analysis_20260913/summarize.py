"""Reduce final-output audits and explicitly separate diagnostic searches."""
from pathlib import Path
import json,hashlib,sys
from collections import Counter
HERE=Path(__file__).resolve().parent
REPO=Path(sys.argv[1]) if len(sys.argv)>1 else HERE.parents[1]
def read(p):return json.loads(p.read_text())
inputs={r['case']:r for r in read(HERE/'input-audit.json')}
original={r['case']:r for r in read(HERE/'original-reference-audit.json')['rows']}
comparators=read(REPO/'manuscript/evidence/competitors.json')['methods']
slap={r['case']:r['outcome'] for r in read(REPO/'manuscript/evidence/slap_sweep.json')['per_case']}
rows=[]
for case,inp in inputs.items():
 directions=[read(HERE/f'results/case{case}-{d}.json') for d in ['R_to_P','P_to_R']]
 assert all(d['checked_cuts']==d['expected_cuts'] for d in directions)
 pairs=sorted(set(int(k) for d in directions for k in d['cardinality_histogram']))
 if all(d['same_cardinality']==0 for d in directions):reason='extra_matched_pair'
 elif all(d['relaxed_target_orbit_pair_matches']==0 for d in directions):reason='atom_correspondence_excluded_by_relaxed_orbits'
 else:reason='joint_correspondence_excluded_by_full_verifier'
 assert original[case]['original_rdf_reference_verified']
 rows.append(dict(case=case,reference_pairs=inp['reference_pairs'],retained_pair_counts=pairs,classification=reason,reference_heavy_pair_changes=inp['reference_change_count'],all_directions_capped=all(d['cap_stops'] for d in directions),any_direction_capped=any(d['cap_stops'] for d in directions),terminals=sum(d['terminals'] for d in directions),source_choice_compatible_terminals=sum(d['source_orbit_selection_matches'] for d in directions),relaxed_orbit_compatible_terminals=sum(d['relaxed_target_orbit_pair_matches'] for d in directions),original_reference_verified=True,comparator_hits=[m['method'] for m in comparators if case in m['any_correct_cases']]+(['slap_sweep'] if slap[case]=='recovered' else [])))
assert Counter(r['classification'] for r in rows)==Counter(extra_matched_pair=4,atom_correspondence_excluded_by_relaxed_orbits=6,joint_correspondence_excluded_by_full_verifier=1)
probes=[read(p) for p in sorted((HERE/'probes').glob('*/result.json'))]
assert len(probes)==8 and all(p['evaluation']['reference_recovery']=='not_recovered' for p in probes)
traces=[read(HERE/f'traces/cause{case}.json') for case in [986,1285]]
assert all(t['native_python_graph_equal'] and all(not r['capped'] for r in t['rows']) for t in traces)
assert all(next(r for r in t['rows'] if r['policy']=='reference_constrained')['exact_reference_witnesses']>0 for t in traces)
assert traces[0]['diagnostic_label_swap']['ordinary_equivalent']
assert traces[1]['reference_selected_multicut']['reference_recovery']=='recovered'
summary=dict(denominator=1851,seed_count=10,branch_cap=100,iso_tolerance=1.0,nonrecovered=11,rows=rows,classifications=dict(Counter(r['classification'] for r in rows)),all_reference_conversions_verified=True,all_misses_have_some_cap_stops=all(r['any_direction_capped'] for r in rows),cases_recovered_by_any_comparator=[r['case'] for r in rows if r['comparator_hits']],diagnostic_probes=[dict(case=p['case'],direction=p['direction'],variant=p['variant'],reference_recovery=p['evaluation']['reference_recovery'],capped=bool(p['metrics']['subtree_branch_cap_count']),search_cpu_seconds=p['search_cpu_seconds']) for p in probes],causal_traces=traces,scope='Audit of final saved ten-seed families. Sensitivity and reference-directed feasibility probes are separate, excluded from benchmark scores. Strict absence is not proof of global unreachability or annotation error.',audit_cpu_seconds=sum(read(p)['cpu_seconds'] for p in (HERE/'results').glob('*.json')),audit_peak_mib=max(read(p)['peak_mib'] for p in (HERE/'results').glob('*.json')))
(HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:summary[k] for k in ['classifications','cases_recovered_by_any_comparator','audit_cpu_seconds','audit_peak_mib']}))
