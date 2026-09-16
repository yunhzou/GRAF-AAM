from pathlib import Path
import json,hashlib,statistics,math
S=Path(__file__).resolve().parent;W=S.parent/'current-validation';rows={i:dict(case=i,cpu_seconds=0.,call_wall_seconds=0.,attempts=[]) for i in range(140)};hashes={}
for a in [1,2]:
 root=W/f'decode-optimized-attempts/{a}';x=json.loads((root/'execution.json').read_text())
 for r in x['rows']:
  c=r['case'];p=root/f'case{c}/result.json';d=json.loads(p.read_text());hashes[str(p.relative_to(W))]=hashlib.sha256(p.read_bytes()).hexdigest();rows[c]['cpu_seconds']+=d['cpu_seconds'];rows[c]['call_wall_seconds']+=d['wall_seconds'];rows[c]['attempts'].append(a)
x=sorted(r['cpu_seconds'] for r in rows.values());q=(len(x)-1)*.95;l=math.floor(q);u=math.ceil(q)
st=dict(n=140,mean=statistics.mean(x),median=statistics.median(x),p95=x[l]+(x[u]-x[l])*(q-l),total=sum(x))
flat=json.loads((S/'manuscript/evidence/final_dedup.json').read_text());assert abs(st['total']-flat['cpu_seconds'])<1e-7
out=dict(scope='Per-reaction complete catalogue and event-window decoding CPU, summed across bounded continuations; includesjournaling, excludesmatching andcompetition. Summedcall wall is not campaign latency.',stats=st,per_case=list(rows.values()),source_hashes=hashes)
p=S/'timing_comparison.json';t=json.loads(p.read_text());t['coordinate_decoding']=out;p.write_text(json.dumps(t,separators=(',',':'))+'\n');print(st)
