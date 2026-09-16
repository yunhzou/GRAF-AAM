import sys,json,os
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'hpc'))
import worker
assert os.environ['SLURM_JOB_PARTITION']=='cpu_short'
tasks=worker.read(ROOT/'hpc/tasks.json')
atoms={r['case']:r['atoms'] for r in worker.read(ROOT/'inputs.json')['records']}
indices=[]
ordinary=[i for i,t in enumerate(tasks) if t['kind']=='aam' and t['seed']==3 and atoms[t['case']]>=20 and not (ROOT/'runs/golden'/f'seed{t["seed"]}'/f'case{t["case"]}'/t['direction']/'search.json').exists()]
indices.append(min(ordinary,key=lambda i:atoms[tasks[i]['case']]))
indices.extend(next(i for i,t in enumerate(tasks) if t['kind']==kind and t['case']==1) for kind in ['slap_score','coordinate_audit'])
rows=[]
for i in indices:
 row=worker.task_run(i,tasks[i]);assert row['status']=='passed',row;rows.append(dict(index=i,task=tasks[i],status=row['status']))
worker.save(ROOT/'hpc/smoke-proof.json',dict(status='passed',rows=rows,scope='Actual HPC wrapper executed one new AAM search/evaluation, one saved SLAP verification, and one relocated coordinate archive audit. Results retained and skipped by the array.'))
print(json.dumps(rows),flush=True)
