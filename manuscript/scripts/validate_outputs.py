"""Validate the final manuscript's numerical sources and compiled references."""
from pathlib import Path
import hashlib,json,re
from pypdf import PdfReader
MAN=Path(__file__).resolve().parents[1];E=MAN/'evidence'
def read(n):return json.loads((E/n).read_text())
seed=read('seed_comparison.json');slap=read('slap_sweep.json');comp=read('competition_final.json');flat=read('final_dedup.json')
assert seed['fresh'] and slap['fresh'] and comp['fresh']
keys=['seeds1','seeds2','seeds3','seeds10']
assert all(seed['methods'][k]['golden_cases']==1851 for k in keys)
assert all(sum(seed['methods'][k]['golden_outcomes'].values())==1851 for k in keys)
for k in keys:
 d=seed['methods'][k]
 assert len(d['per_case'])==1851 and {r['case'] for r in d['per_case']}==set(range(1851))
 assert all(r['outcome'] in {'recovered','not_recovered','unknown'} for r in d['per_case'])
 assert d['golden_outcomes']['recovered']==sum(r['outcome']=='recovered' for r in d['per_case'])
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
result=dict(status='passed',pages=len(pages),figures=figs,source_checks=True,references_resolved=True,no_overfull_boxes=True,
 manual_visual_review_required=True,scope='Numerical and build validation; visual review is recorded separately. All paper benchmark configurations were rerun; completeness and limits are recorded in the evidence.',
 manuscript_sha256=hashlib.sha256((MAN/'manuscript.pdf').read_bytes()).hexdigest())
(MAN/'build/artifact-validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
