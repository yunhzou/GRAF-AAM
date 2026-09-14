"""Validate the final manuscript's numerical sources and compiled references."""
from pathlib import Path
import hashlib,json,re
from pypdf import PdfReader
MAN=Path(__file__).resolve().parents[1];E=MAN/'evidence'
def read(n):return json.loads((E/n).read_text())
from build_molecule_figure import check_example
assert read('molecule_example.json') == json.loads(json.dumps(check_example()))
uncut=read('unswept.json')
assert uncut['denominator']==1851 and uncut['seed_count']==1
for name in ['graft','slap']:
 d=uncut['methods'][name]
 assert sum(d['counts'].values())==1851
 assert d['counts']['recovered']==len(set(d['recovered_cases']))
 assert d['sweep_recovered']==len(d['recovered_cases'])+len(d['added_by_sweep'])
competitors=read('competitors.json')
assert competitors['denominator']==1851 and competitors['rescored_with_current_evaluator'] and competitors['unchanged_case_outcomes']
assert len(competitors['methods'])==7
for d in competitors['methods']:
 assert d['total']==1851 and sum(d['statuses'].values())==1851
 assert d['first_correct']==len(set(d['first_correct_cases']))
 assert d['any_correct']==len(set(d['any_correct_cases']))
 assert set(d['first_correct_cases'])<=set(d['any_correct_cases'])<=set(range(1851))
 assert len(d['failed_cases'])==sum(n for k,n in d['statuses'].items() if k!='mapped')
 assert d['invalid_candidates']==sum(d['invalid_reasons'].values())
seed=read('seed_comparison.json');slap=read('slap_sweep.json');comp=read('competition_final.json');flat=read('final_dedup.json')
misses=read('golden_miss_analysis.json')
assert misses['nonrecovered']==11 and len(misses['rows'])==11
assert {r['case'] for r in misses['rows']}=={r['case'] for r in seed['methods']['seeds10']['per_case'] if r['outcome']=='not_recovered'}
assert sum(r['classification']=='extra_matched_pair' for r in misses['rows'])==4
assert all(r['retained_pair_counts']==[r['reference_pairs']+1] for r in misses['rows'] if r['classification']=='extra_matched_pair')
assert misses['all_reference_conversions_verified'] and misses['cases_recovered_by_any_comparator']==[871]
assert all(t['native_python_graph_equal'] for t in misses['causal_traces'])

assert seed['fresh'] and slap['fresh'] and comp['fresh']
keys=['seeds1','seeds2','seeds3','seeds10']
assert all(seed['methods'][k]['golden_cases']==1851 for k in keys)
assert all(sum(seed['methods'][k]['golden_outcomes'].values())==1851 for k in keys)
for k in keys:
 d=seed['methods'][k]
 assert len(d['per_case'])==1851 and {r['case'] for r in d['per_case']}==set(range(1851))
 assert all(r['outcome'] in {'recovered','not_recovered','unknown'} for r in d['per_case'])
 for status in ('recovered','not_recovered','unknown'):
  assert d['golden_outcomes'][status]==sum(r['outcome']==status for r in d['per_case'])
 if d['common_mean_cpu_seconds'] is not None:
  assert abs(d['common_mean_cpu_seconds']-sum(d['per_case'][c]['search_cpu_including_io'] for c in seed['common_case_indices'])/len(seed['common_case_indices']))<1e-8
  assert all(d['per_case'][c]['search_complete'] and d['per_case'][c]['search_hosts']==['Mac'] for c in seed['common_case_indices'])
assert sum(slap['outcomes'].values())==1851
assert slap['sweep_union_recovered']==slap['outcomes']['recovered']
assert comp['comparisons']['slap_sweep']['union_classes']==166 and comp['comparisons']['slap_sweep']['total_classes']==168
assert comp['comparisons']['slap_sweep']['union_complete_cases']==139
assert comp['comparisons']['native_slap']['union_classes']==155 and comp['comparisons']['native_slap']['total_classes']==160
assert comp['new_window_class_count']==36 and len(comp['new_window_cases'])==9
assert flat['old_branches']==237645 and flat['new_branches']==124641 and flat['flat_families']==236653
assert flat['old_median']==596.5 and flat['new_median']==397 and flat['workers']==3
assert all(r['complete'] for r in flat['per_case'])
assert len(flat['per_case'])==140 and flat['fresh_decode_all140'] is True
assert sum(flat['window_distribution'].values())==140
assert flat['flat_families']==sum(r['flat_saved_families'] for r in flat['per_case'])
assert flat['input_paths']==sum(r['input_paths'] for r in flat['per_case'])
sources=json.loads((E/'paper_sources.json').read_text())
assert sources['fresh_campaign_complete']
for row in sources['snapshots']:
 assert hashlib.sha256((E/row['snapshot']).read_bytes()).hexdigest()==row['sha256'],row['snapshot']
tex='\n'.join(p.read_text() for p in (MAN/'includes').glob('*.tex') if p.name in ['paper.tex','supplement.tex','include-abstract.tex'])
for obsolete in ['Adaptive no-sweep','Separate experimental versions','Holdout follow-ups','earlier publication engine','N_{\\mathrm{order\\ changed}}']:
 assert obsolete not in tex,obsolete
bib=(MAN/'references.bib').read_text()
for group in re.findall(r'\\cite\w*\{([^}]+)\}',tex):
 for key in group.split(','):assert re.search(r'@\w+\{'+re.escape(key)+',',bib),key
figs=re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}',tex)
assert set(figs)=={'figs/fig1_algorithm.pdf','figs/fig2_golden.pdf','figs/fig3_coordinate.pdf'}
for fig in figs:
 p=MAN/fig
 assert len(PdfReader(p).pages)==1
 for ext in ['svg','png']:assert p.with_suffix('.'+ext).is_file()
log=(MAN/'build/preprint.log').read_text()
assert not re.search(r'(?:Citation|Reference).*undefined|There were undefined|Overfull \\[hv]box|Missing character:',log), 'Inspect TeX log'
r=PdfReader(MAN/'manuscript.pdf');pages=[p.extract_text() for p in r.pages];alltext='\n'.join(pages)
assert '??' not in alltext
for value in [f"{seed['methods']['seeds1']['golden_recovery_percent']:.2f}",f"{seed['methods']['seeds3']['golden_recovery_percent']:.2f}",'124,641','237,645',f"{flat['wall_seconds']/60:.2f}",'166','168']:
 assert value in alltext,value
assert 'GRAFT and benchmarked mapping implementations on Golden' in alltext
assert 'First correct' not in alltext and 'two ten-order cases remain unresolved' not in alltext
comparison=(MAN/'includes/generated-competitor-table.tex').read_text()
assert comparison.count('GRAFT &')==5
assert comparison.count('SLAP,')==4
for fig in figs:
 assert 'GRAFT' in PdfReader(MAN/fig).pages[0].extract_text(), fig
for value in ['1,489','1,661','80.44','89.74']:assert value in alltext
for d in competitors['methods']:
 assert f"{d['any_correct']:,}" in alltext

direction=read('direction_recovery.json')
assert direction['no_search_reruns'] and direction['denominator']==1851
for key,d in direction['methods'].items():
 assert len(d['per_case'])==1851
 for group,g in d['groups'].items():
  rs=[r for r in d['per_case'] if group=='all' or r['cohort']==group]
  assert g['n']==len(rs)
  for policy in ['smaller_first','larger_first','bidirectional']:
   assert g[policy]['cases']==[r['case'] for r in rs if r[policy]=='recovered']
   assert g[policy]['recovered']==len(g[policy]['cases']) and g[policy]['unknown']==0
 assert d['groups']['all']['bidirectional']['recovered']==seed['methods'][key]['golden_outcomes']['recovered']
for value in ['1,803','1,779','97.41','96.11','Smaller-first','Larger-first','(default)']:assert value in alltext,value

timing=read('timing_comparison.json')
assert timing['searches_rerun']==0 and timing['common_case_indices']==seed['common_case_indices']
import statistics
for k,d in timing['methods'].items():
 rows=d['per_case'];assert [r['case'] for r in rows]==seed['common_case_indices']
 xs=sorted(r['cpu_seconds'] for r in rows);st=d['stats'];assert st['n']==1821
 assert abs(st['mean']-statistics.mean(xs))<1e-8 and abs(st['median']-statistics.median(xs))<1e-8
 # 1821 samples give an integer index under the stated linear percentile rule.
 assert abs(st['p95']-xs[1729])<1e-8
for prefix,key in [('graft1','seeds1'),('graft2','seeds2')]:
 assert abs(timing['methods'][prefix+'_bidirectional']['stats']['mean']-seed['methods'][key]['common_mean_cpu_seconds'])<1e-8
for v in ['1.141','1.377','2.518','3.329','1.490','4.819','95th pct.','Mean wall s']:assert v in alltext,v
for d in timing['default_comparators']:
 source=next(r for r in competitors['methods'] if r['method']==d['key'])
 assert d['calls']==source['statuses']['mapped']
 assert abs(d['mean_wall_seconds']*d['calls']-source['original_mapping_timing']['successful_mapping_wall_sum_seconds'])<1e-8
coord=timing['coordinate_decoding']
assert len(coord['per_case'])==140
assert abs(sum(r['cpu_seconds'] for r in coord['per_case'])-flat['cpu_seconds'])<1e-7
assert abs(coord['stats']['mean']-flat['cpu_seconds']/140)<1e-8
for value in ['11.01','1.22','52.67']:assert value in alltext,value
result=dict(status='passed',pages=len(pages),figures=figs,source_checks=True,references_resolved=True,no_overfull_boxes=True,
 manual_visual_review_required=True,scope='Numerical and build validation; visual review is recorded separately. GRAFT and the SLAP sweep use final-source campaigns; archived default-comparator outputs were rescored with the current evaluator. Completeness and limits are recorded in the evidence.',
 manuscript_sha256=hashlib.sha256((MAN/'manuscript.pdf').read_bytes()).hexdigest())
(MAN/'build/artifact-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
