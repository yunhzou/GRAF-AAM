from concurrent.futures import ThreadPoolExecutor
from run import launch,save,S
if __name__ == '__main__':
 with ThreadPoolExecutor(max_workers=2) as p: rows=list(p.map(lambda c:launch(c,'prefix512','search'),[59,64]))
 with ThreadPoolExecutor(max_workers=2) as p: rows+=list(p.map(lambda c:launch(c,'prefix512','decode'),[r['case'] for r in rows if r['status']=='passed']))
 save(S/'budget-execution.json',rows)
