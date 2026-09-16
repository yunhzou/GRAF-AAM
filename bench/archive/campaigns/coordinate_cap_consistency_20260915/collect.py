"""Collect stage records and distinguish empty searches from unfinished decoding."""
from pathlib import Path
import json
S=Path(__file__).resolve().parent

def classify(d):
 d['missing_minimum_counts']=[r['case'] for r in d['per_case'] if not r['both_resolved']]
 d['unresolved']=[r['case'] for r in d['per_case'] if any(not r['caps'][c]['decoded'] or not r['caps'][c]['decoded']['complete_saved_full_window'] for c in ['100','2000'])]
 d['no_full_mapping']={c:[r['case'] for r in d['per_case'] if r['caps'][c]['decoded'] and r['caps'][c]['decoded']['no_full_families']] for c in ['100','2000']}
 for r in d['per_case']:
  a,b=[r['caps'][c]['decoded'] for c in ['100','2000']]
  if a and b and a['complete_saved_full_window'] and b['complete_saved_full_window']:
   r['cap2000_minimum_patterns_missing_at_cap100']=sorted(set(b['minimum_ids'])-set(a['patterns']))
   r['cap2000_window_patterns_missing_at_cap100']=sorted(set(b['patterns'])-set(a['patterns']))
 d['cap100_missing_minimum_patterns']=[r['case'] for r in d['per_case'] if r.get('cap2000_minimum_patterns_missing_at_cap100')]
 return d

def main():
 from aggregate import aggregate
 d=aggregate()
 for r in d['per_case']:
  root=S/'results'/str(r['case']);attempt=root/json.loads((root/'selected.json').read_text())['attempt']
  for cap,item in r['caps'].items():
   for stage in ['search','competition']:
    p=attempt/cap/(stage+'.json');item[stage]=json.loads(p.read_text()) if p.exists() else None
 (S/'comparison.json').write_text(json.dumps(classify(d))+'\n')
if __name__=='__main__':main()
