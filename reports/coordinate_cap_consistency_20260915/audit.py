"""Validate paired configurations and decoded records before reporting."""
import json,hashlib
from pathlib import Path
S=Path(__file__).resolve().parent
d=json.loads((S/'comparison.json').read_text());seen=set();checked=0
for row in d['per_case']:
 assert row['case'] not in seen;seen.add(row['case'])
 configs=[];competition=[]
 for cap in ('100','2000'):
  r=row['caps'][cap];s=r.get('search')
  if s:
   cfg=dict(s['config']);assert cfg.pop('branch_limit')==int(cap);configs.append(cfg)
  if r.get('competition'):competition.append(r['competition']['config'])
  x=r['decoded']
  if not x:continue
  assert x['cap']==int(cap) and x['case']==row['case']
  assert x['completed_families']<=x['families']
  assert x['complete_saved_full_window']==(x['completed_families']==x['families'])
  assert x['no_full_families']==(x['families']==0)
  totals=[]
  for key,p in x['patterns'].items():
   assert key==p['id'];assert len(p['mapping'])==len(set(p['mapping']))
   assert p['total']==len(p['events']['broken'])+len(p['events']['formed'])
   assert p['total']<=x['window'];totals.append(p['total'])
  assert x['minimum']==min(totals,default=None)
  assert x['minimum_ids']==sorted(k for k,p in x['patterns'].items() if p['total']==x['minimum'])
  assert x['minimum_proven']==(x['complete_saved_full_window'] and bool(totals))
  if x['no_full_families']:assert not x['patterns']
  checked+=1
 if len(configs)==2:assert configs[0]==configs[1]
 if len(competition)==2:assert competition[0]==competition[1]
assert len(seen)==d['processed']
proof=dict(status='passed',processed=d['processed'],decoded_records_checked=checked,comparison_sha256=hashlib.sha256((S/'comparison.json').read_bytes()).hexdigest(),scope='Paired configurations differ only in branch limit; competition budgets agree; stored event counts, unique mapping images, minima and completeness flags are internally consistent. Witnesses were independently rescored against original endpoint matrices by run.py.')
(S/'audit.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))
