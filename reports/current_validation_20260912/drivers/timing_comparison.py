import json,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
v=ROOT/'hpc-return/first-pass'
def read(p):return json.loads(p.read_text()) if p.exists() else {}
rows=[];excluded=[]
for c in range(1851):
 s=read(v/'slap-golden/evaluations'/f'{c}.json')
 values={};ok=s.get('complete',False)
 for seed in [1,2]:
  total=0
  for direction in ['R_to_P','P_to_R']:
   f=v/'runs/golden'/f'seed{seed}'/f'case{c}'/direction
   e=read(f/'search_execution.json');r=read(f/'search.json')
   if e.get('status')!='passed' or e.get('platform')=='linux' or e.get('host'):ok=False
   total+=r.get('cpu_seconds',0)
  values[f'aam_seed{seed}']=total
 if ok:rows.append(dict(case=c,slap=s['workflow_cpu_excluding_io'],**values))
 else:excluded.append(c)
metrics={k:dict(mean_seconds=statistics.mean(r[k] for r in rows),median_seconds=statistics.median(r[k] for r in rows),total_cpu_minutes=sum(r[k] for r in rows)/60) for k in ['aam_seed1','aam_seed2','slap']}
result=dict(cases=len(rows),excluded_cases=excluded,metrics=metrics,scope='Same-reaction Mac cohort: both AAM directions completed on Mac at seed settings 1 and 2; SLAP mapping complete. AAM search CPU includes finalizing/compressing its saved output. SLAP workflow CPU includes graph construction, mapping and export but excludes IO. Reference-membership verification excluded for both. SLAP uses both directions, binary/weighted modes, uncut and all single-edge cuts. Not a comparison of identical output products; higher-seed mixed-host timings excluded.')
(ROOT/'hpc/same-mac-timing.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
