from concurrent.futures import ThreadPoolExecutor
from run import launch,save,S
jobs=[(c,'local') for c in [6,11,59,64,101]]+[(c,'prefix') for c in [6,11,59,64,101]]+[(11,'extended')]
with ThreadPoolExecutor(max_workers=4) as p:rows=list(p.map(lambda j:launch(*j,'search'),jobs))
assert all(x['status']=='passed' for x in rows),rows
with ThreadPoolExecutor(max_workers=4) as p:rows+=list(p.map(lambda j:launch(*j,'decode'),[(c,'prefix') for c in [6,11,59,64,101]]+[(59,'local'),(11,'extended')]))
save(S/'execution.json',rows)
