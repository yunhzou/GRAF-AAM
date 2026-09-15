"""Aggregate every attempted case, retaining resource limits and paired timing."""
import sys,json,collections,statistics,hashlib,gzip
from pathlib import Path
import run
S=Path(__file__).resolve().parent

def stats(values):
 values=sorted(values)
 if not values:return None
 def pct(q):
  i=(len(values)-1)*q;lo=int(i);hi=min(lo+1,len(values)-1);return values[lo]+(values[hi]-values[lo])*(i-lo)
 return dict(n=len(values),sum_cpu_seconds=sum(values),mean_cpu_seconds=statistics.mean(values),median_cpu_seconds=statistics.median(values),p95_cpu_seconds=pct(.95))

def aggregate():
 present={int(p.name) for p in (S/'results').iterdir() if p.is_dir() and p.name.isdigit()} if (S/'results').exists() else set()
 configs=run.configurations();out={c['key']:dict(config=c,per_case=[]) for c in configs};machines={};done=[];phases=collections.Counter()
 for case in range(1851):
  root=S/'results'/str(case)
  if case in present and (root/'done.json').exists():done.append(case)
  if case in present and (root/'machine.json').exists():
   m=run.read(root/'machine.json');cpu=json.loads(m['cpu_info']).get('lscpu',[]) if m.get('cpu_info','').startswith('{') else []
   model=next((r['data'] for r in cpu if r['field']=='Model name:'),'unknown');machines[case]=dict(host=m['host'],model=model,threads=m['threads'])
  for cfg in configs:
   dirs={}
   for d in ['R_to_P','P_to_R']:
    p=run.folder(case,cfg['key'],d);files={n:run.read(p/n) for n in ['search-execution.json','score-execution.json','search.json','evaluation.json'] if case in present and (p/n).exists()}
    execution=files.get('search-execution.json',{});search=files.get('search.json',{});score=files.get('evaluation.json',{})
    phases[cfg['key'],execution.get('status','missing')]+=1
    dirs[d]=dict(search_status=execution.get('status','missing'),score_status=files.get('score-execution.json',{}).get('status','missing'),search_complete=bool(execution.get('status')=='passed' and search.get('search_complete')),reference_recovery=score.get('reference_recovery','unknown'),api_cpu_seconds=execution.get('api_cpu_seconds'),worker_cpu_seconds=execution.get('worker_cpu_seconds'),capped=search.get('capped'),host=execution.get('host'),peak_mib=execution.get('peak_mib'),compute=search.get('compute'),stages=search.get('stages'))
   verdicts=[r['reference_recovery'] for r in dirs.values()]
   verdict='recovered' if 'recovered' in verdicts else 'not_recovered' if all(v=='not_recovered' for v in verdicts) else 'unknown'
   row=dict(case=case,reference_recovery=verdict,directions=dirs,search_complete=all(r['search_complete'] for r in dirs.values()),api_cpu_seconds=sum(r['api_cpu_seconds'] for r in dirs.values()) if all(r['api_cpu_seconds'] is not None for r in dirs.values()) else None,worker_cpu_seconds=sum(r['worker_cpu_seconds'] for r in dirs.values()) if all(r['worker_cpu_seconds'] is not None for r in dirs.values()) else None)
   out[cfg['key']]['per_case'].append(row)
 common=[i for i in range(1851) if all(v['per_case'][i]['search_complete'] for v in out.values())]
 models=collections.Counter(m['model'] for m in machines.values());hosts=collections.Counter(m['host'] for m in machines.values())
 for key,v in out.items():
  v['outcomes']=dict(collections.Counter(r['reference_recovery'] for r in v['per_case']))
  v['completed_search_cases']=sum(r['search_complete'] for r in v['per_case'])
  v['all_attempted_api_cpu']=stats([r['api_cpu_seconds'] for r in v['per_case'] if r['api_cpu_seconds'] is not None])
  v['all_attempted_worker_cpu']=stats([r['worker_cpu_seconds'] for r in v['per_case'] if r['worker_cpu_seconds'] is not None])
  v['paired_completed_api_cpu']=stats([v['per_case'][i]['api_cpu_seconds'] for i in common])
  v['paired_directions']={d:stats([v['per_case'][i]['directions'][d]['api_cpu_seconds'] for i in common]) for d in ['R_to_P','P_to_R']}
 result=dict(status='all_tasks_processed' if len(done)==1851 else 'in_progress',processed_cases=len(done),case_denominator=1851,common_completed_cases=common,methods=out,cpu_models=dict(models),hosts=dict(hosts),same_cpu_model=len(models)==1 and 'unknown' not in models,all_threads_one=all(all(v=='1' for v in m['threads'].values()) for m in machines.values()),source_manifest_sha256=run.sha(S/'manifest.json'),timing_scope='Fresh API search including graph preparation and native output persistence, excluding imports and verification. All-attempted CPU includes interrupted calls with recorded baselines; paired completed timing uses identical reactions across all twelve configurations. Worker CPU additionally includes startup and imports. No reference-guided early stopping of search.')
 run.save(S/'progress.json',dict(status=result['status'],processed_cases=len(done),common_completed_cases=len(common),cpu_models=dict(models),methods={k:dict(outcomes=v['outcomes'],search_complete=v['completed_search_cases'],api=v['paired_completed_api_cpu']) for k,v in out.items()}))
 with gzip.open(S/'controlled-results.json.gz','wt') as f:json.dump(result,f)
 print(json.dumps(run.read(S/'progress.json')),flush=True)
 return result
if __name__=='__main__':aggregate()
