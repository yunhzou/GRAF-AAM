"""Slurm dispatch around unchanged benchmark children; one CPU per worker."""
import sys,os,json,time,shutil,socket,hashlib,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import run
read,save,sha=run.read,run.save,run.sha
ENV=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',PYTHONHASHSEED='0',PYTHONDONTWRITEBYTECODE='1')
def execute(command,log,limit=300,memory=3072):
 guard=run.MemoryGuard(memory,memory,128)
 log.parent.mkdir(parents=True,exist_ok=True)
 start=time.perf_counter()
 with log.open('w') as f:status,peak=run.bounded_process(command,ENV,f,limit,guard)
 return dict(status=status,peak_mib=peak,elapsed_seconds=time.perf_counter()-start,host=socket.gethostname(),job_id=os.environ.get('SLURM_JOB_ID'),platform='linux',memory_limit_mib=memory,watchdog_seconds=limit,runtime_sha256=sha(ROOT/'hpc/runtime.json'),engine_sha256=sha(ROOT/'hpc/linux-engine.json'))
def backup(folder,phase,attempt):
 dest=ROOT/'hpc/prior-attempts'/folder.relative_to(ROOT)/str(attempt)/phase
 dest.mkdir(parents=True,exist_ok=True)
 for name in [f'{phase}_execution.json',f'{phase}.log','search.json' if phase=='search' else 'evaluation.json']:
  if (folder/name).exists():shutil.copy2(folder/name,dest/name)
 return str(dest)
def aam(task):
 seed,case,direction=task['seed'],task['case'],task['direction']
 folder=ROOT/'runs/golden'/f'seed{seed}'/f'case{case}'/direction
 folder.mkdir(parents=True,exist_ok=True);attempts=[]
 def state(phase):
  p=folder/f'{phase}_execution.json'
  return read(p) if p.exists() else {}
 # Verify completed raw cuts first if a prior search was interrupted.
 previous=state('search')
 if previous.get('status') not in [None,'passed'] and (folder/'cuts/manifest.json').exists():
  row=execute([sys.executable,str(ROOT/'score_checkpoints.py'),str(seed),str(case),direction],folder/'checkpoint-score-hpc.log',memory=6144)
  row.update(phase='checkpoint_score',previous_attempt=backup(folder,'score','before-checkpoint'))
  attempts.append(row)
  p=folder/'checkpoint-verification.json'
  if row['status']=='passed' and p.exists() and read(p).get('reference_recovery')=='recovered':
   save(folder/'evaluation.json',read(p));save(folder/'score_execution.json',row)
   return dict(task=task,status='passed',scope='Positive cut witness; interrupted search remains incomplete',attempts=attempts)
 for phase,target in [('search','search.json'),('score','evaluation.json')]:
  prior=state(phase)
  if prior.get('status')=='passed' and (folder/target).exists():continue
  # One initial attempt plus at most two checkpoint continuations; prior Mac attempt counts.
  start_attempt=1 if prior else 0
  for attempt in range(start_attempt,3):
   old=backup(folder,phase,attempt)
   row=execute([sys.executable,str(ROOT/'run.py'),'child','golden',str(seed),str(case),direction,phase],folder/f'{phase}.log',memory=3072 if attempt==0 else 6144)
   row.update(phase=phase,attempt=attempt,previous_attempt=old)
   save(folder/f'{phase}_execution.json',row);attempts.append(row)
   if row['status']=='passed':break
   if row['status']=='failed':break # Code errors require inspection, not blind repeats.
  if row['status']!='passed':return dict(task=task,status=row['status'],attempts=attempts)
 return dict(task=task,status='passed',attempts=attempts)
def coordinate_child(case):
 import audit_coordinate as audit
 original=audit.read
 local=read(ROOT/'hpc/transfer-plan.json')['local_root']
 def relocated_read(path):
  value=original(path)
  if path.name=='identity.json' and 'decoded-optimized' in path.parts:
   value=dict(value,archives={str(ROOT/Path(k).relative_to(local)):v for k,v in value['archives'].items()})
  return value
 audit.read=relocated_read
 # The setup job already verifies the exact cut-finalization self-test.
 if case==0:
  real_run=subprocess.run
  def skip_completed_test(command,*args,**kwargs):
   if len(command)>=3 and command[1]==str(ROOT/'score_checkpoints.py') and command[2]=='test':
    proof=read(ROOT/'checkpoint-verification-tests.json');assert proof['status']=='passed' and proof['driver_sha256']==sha(ROOT/'score_checkpoints.py')
    return subprocess.CompletedProcess(command,0)
   return real_run(command,*args,**kwargs)
  subprocess.run=skip_completed_test
 audit.child(case)
def task_run(index,task):
 dest=ROOT/'hpc/tasks'/f'{index}.json'
 if dest.exists():return read(dest)
 if task['kind']=='aam':result=aam(task)
 elif task['kind']=='slap_score':
  c=task['case'];row=execute([sys.executable,str(ROOT/'slap_golden.py'),'evaluate',str(c)],ROOT/'slap-golden/evaluations'/f'{c}.hpc.log')
  result=dict(task=task,**row)
 elif task['kind']=='coordinate_audit':
  c=task['case'];row=execute([sys.executable,__file__,'coordinate_child',str(c)],ROOT/'coordinate-audit'/f'case{c}.hpc.log',memory=6144)
  result=dict(task=task,**row)
 else:raise ValueError(task)
 save(dest,result);return result
if __name__=='__main__':
 if sys.argv[1]=='coordinate_child':coordinate_child(int(sys.argv[2]))
 else:
  assert os.environ.get('SLURM_JOB_PARTITION')=='cpu_short'
  assert int(os.environ.get('SLURM_CPUS_PER_TASK','1'))==1
  assert read(ROOT/'hpc/build-proof.json')['status']=='passed'
  rank=int(sys.argv[1]);n=int(sys.argv[2]);tasks=read(ROOT/'hpc/tasks.json')
  for index in range(rank,len(tasks),n):
   result=task_run(index,tasks[index]);print(json.dumps(dict(index=index,status=result['status'],task=tasks[index])),flush=True)
