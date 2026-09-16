from paired_control import *
import statistics

def measure(a,trial,mode):
 t=time.process_time();r=compete_fragments(a,CompetitionConfig(proposal_mode=mode,operation_budget=100000,depth_limit=1,seconds=270),trial=trial);search=time.process_time()-t
 assert r.pending==0
 cat=FinalBranchCatalogue(a.problem)
 for i,repair in enumerate(r.repairs):cat.add_aam(repair,str(i))
 idx=SignedEventIndex(a.problem);t=time.process_time()
 for f in cat.families:
  z=extract_path_events(f.as_path(a.problem),a.problem,idx,max_events=8 if trial['case']==64 else 5,seconds=20,max_patterns=None)
  assert z['complete']
 return dict(search_cpu=search,decode_cpu=time.process_time()-t,families=len(cat.families),completion_calls=r.counts['completion_calls'])
rows=[]
for case in [6,11,59,64,101]:
 x=json.loads((S/f'paired-control/{case}/result.json').read_text());trial=dict(x['trial'],case=case);a=read_aam_checkpoint(base/f'inputs/case{case}.pkl.gz')
 for mode in ['extended','prefix']:measure(a,trial,mode)
 samples={m:[] for m in ['extended','prefix']}
 for repeat in range(3):
  for mode in (['extended','prefix'] if repeat%2==0 else ['prefix','extended']):samples[mode].append(measure(a,trial,mode))
 result=dict(case=case,samples=samples,medians={m:{k:statistics.median(x[k] for x in samples[m]) for k in ['search_cpu','decode_cpu','families','completion_calls']} for m in samples})
 save(S/f'paired-control/{case}/timing.json',result);rows.append(result);print(case,result['medians'],flush=True)
save(S/'paired-control/timing.json',rows)
