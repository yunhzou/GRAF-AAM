"""Summarize a mixed-host continuation without blending machine timings."""
import sys,json,collections,statistics,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));from run import read,save

def optional(p):return read(p) if p.exists() else {}
tasks=read(ROOT/'hpc/tasks.json');records=[optional(ROOT/'hpc/tasks'/f'{i}.json') for i in range(len(tasks))]
missing=[i for i,r in enumerate(records) if not r]
failures=[dict(index=i,**r) for i,r in enumerate(records) if r and r['status']!='passed']
golden={}
for seed in [1,2,3,10]:
 rows=[]
 for case in range(1851):
  directions=[]
  for direction in ['R_to_P','P_to_R']:
   f=ROOT/'runs/golden'/f'seed{seed}'/f'case{case}'/direction
   se=optional(f/'search_execution.json');ee=optional(f/'score_execution.json');ev=optional(f/'evaluation.json')
   outcome=ev.get('reference_recovery','unknown') if ee.get('status')=='passed' else 'unknown'
   directions.append(dict(direction=direction,outcome=outcome,search_status=se.get('status','missing'),score_status=ee.get('status','missing'),host=ee.get('host',se.get('host','Mac'))))
  outcomes=[r['outcome'] for r in directions]
  result='recovered' if 'recovered' in outcomes else 'not_recovered' if outcomes==['not_recovered','not_recovered'] else 'unknown'
  rows.append(dict(case=case,outcome=result,directions=directions))
 golden[str(seed)]=dict(counts=dict(collections.Counter(r['outcome'] for r in rows)),rows=rows)
slap=[]
for i,t in enumerate(tasks):
 if t['kind']=='slap_score':
  c=t['case'];r=optional(ROOT/'slap-golden/evaluations'/f'{c}.json')
  slap.append(dict(case=c,outcome=r.get('reference_recovery','unknown') if records[i].get('status')=='passed' else 'unknown'))
summary=dict(created=time.time(),tasks=len(tasks),completed=len(tasks)-len(missing),missing_task_indices=missing,resource_or_execution_failures=failures,golden=golden,slap=dict(counts=dict(collections.Counter(r['outcome'] for r in slap)),rows=slap),timing_scope='Mac and Linux measurements are separate; no pooled speedup estimate. Per-attempt execution records retain host and runtime hashes.')
audit=[optional(ROOT/'coordinate-audit'/f'case{c}.json') for c in range(140)]
audit_tasks=[records[i] for i,t in enumerate(tasks) if t['kind']=='coordinate_audit']
if all(r.get('status')=='passed' for r in audit_tasks) and all(audit):
 result=dict(cases=140,per_case=audit,old_branches=sum(r['old_literal_branches'] for r in audit),new_branches=sum(r['final_branches'] for r in audit),old_median=statistics.median(r['old_literal_branches'] for r in audit),new_median=statistics.median(r['final_branches'] for r in audit),new_window_class_count=sum(len(r['new_classes']) for r in audit),new_window_cases=[r['case'] for r in audit if r['new_classes']],scope='Fresh archive audit on Linux; archive identities translated by path only with original hashes verified.')
 save(ROOT/'coordinate-audit/summary.json',result);summary['coordinate_audit']={k:v for k,v in result.items() if k!='per_case'}
save(ROOT/'hpc/result-summary.json',summary)
print(json.dumps({k:v for k,v in summary.items() if k not in ['golden','slap','resource_or_execution_failures']},indent=2))
print('Golden:',{s:v['counts'] for s,v in golden.items()},'SLAP:',summary['slap']['counts'])
if missing:raise SystemExit(2)

# Compact return bundle; large raw search archives remain in the HPC project.
import tarfile,hashlib
paths=set()
for p in (ROOT/'runs/golden').glob('seed*/case*/*/*'):
 if p.is_file() and p.suffix in ['.json','.log']:paths.add(p)
for base in [ROOT/'coordinate-audit',ROOT/'slap-golden/evaluations',ROOT/'hpc/tasks',ROOT/'hpc/prior-attempts',ROOT/'hpc/late-cuts']:
 paths.update(p for p in base.rglob('*') if p.is_file() and p.suffix in ['.json','.log'])
paths.update(p for p in (ROOT/'hpc').iterdir() if p.is_file() and p.suffix in ['.json','.py','.sbatch','.md'])
paths.add(ROOT/'checkpoint-verification-tests.json')
archive=ROOT/'hpc/results.tar.gz'
with tarfile.open(archive,'w:gz',compresslevel=1) as tar:
 for p in sorted(paths):tar.add(p,arcname=str(p.relative_to(ROOT)),recursive=False)
digest=hashlib.sha256()
with archive.open('rb') as f:
 for block in iter(lambda:f.read(8*1024*1024),b''):digest.update(block)
(ROOT/'hpc/results.sha256').write_text(digest.hexdigest()+'  results.tar.gz\n')
print('Return bundle:',archive,'bytes:',archive.stat().st_size,flush=True)
