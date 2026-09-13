"""Rescore saved predictions; no mapping algorithms, model loads, or searches."""
from pathlib import Path
import os,sys,json,hashlib,time,signal
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
REPO=Path('/Users/yunhengz/Desktop/AAM Writing')
V=Path(__file__).resolve().parent
sys.path[:0]=[str(REPO/'bench'),str(REPO/'src')]
METHODS=['rxnmapper','localmapper','chython','slap_binary','slap_weighted','indigo','rdt']
VERSIONS=['0.4.3','0.1.5','2.18 / chython-rxnmap 2.0','1.0.0 binary','1.0.0 weighted','1.46.0','4.0.0']

def worker(method):
 # The five-minute cap also covers signature construction on unusual graphs.
 signal.alarm(300)
 from rdkit import RDLogger
 RDLogger.DisableLog('rdApp.*')
 from golden_competitors import signatures
 data=json.loads((V/'archive-metadata.json').read_text())
 audit=[json.loads(s) for s in (V.parent/'current-validation/golden-data/audit.jsonl').read_text().splitlines()]
 old={r['index']:r for r in data[method+'_evaluation.json']}
 start=time.perf_counter();cpu=time.process_time();records=[];changes=[]
 for row in audit:
  i=row['index'];raw=data[f'{method}/{i}.json'];record={'index':i,'status':raw['status']}
  assert raw['index']==i and raw['method']==method
  if raw['status']=='mapped':
   endpoints,expected=signatures(row['mapped_reaction']);hits=[];invalid=[]
   for j,candidate in enumerate(raw['candidates']):
    try:
     actual_endpoints,actual=signatures(candidate['mapped_rxn'])
     if actual_endpoints!=endpoints:invalid.append({'candidate':j,'error':'Endpoint chemistry changed'});continue
     if actual==expected:hits.append(j)
    except Exception as e:invalid.append({'candidate':j,'error':str(e)})
   record.update(first_correct=0 in hits,any_correct=bool(hits),matching_candidate_indices=hits,candidates=len(raw['candidates']),invalid=invalid)
  if record!=old[i]:changes.append({'index':i,'old':old[i],'current':record})
  records.append(record)
 result={'method':method,'records':records,'changes':changes,'elapsed_seconds':time.perf_counter()-start,'cpu_seconds':time.process_time()-cpu}
 (V/f'{method}-recheck.json').write_text(json.dumps(result,indent=2)+'\n')
 signal.alarm(0)
 return {k:v for k,v in result.items() if k!='records'}

if __name__=='__main__':
 start=time.perf_counter()
 with ProcessPoolExecutor(max_workers=4) as pool:
  for result in pool.map(worker,METHODS):print(json.dumps(result),flush=True)
 print(json.dumps({'total_elapsed_seconds':time.perf_counter()-start}),flush=True)
