import os,sys,json,platform,hashlib,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import run
if len(sys.argv)>1 and sys.argv[1]=='parity':
 run.ROOT=ROOT/'hpc/parity'
 run.ROOT.mkdir(exist_ok=True)
 (run.ROOT/'golden-inputs').symlink_to(ROOT/'golden-inputs',target_is_directory=True) if not (run.ROOT/'golden-inputs').exists() else None
 case=int(sys.argv[2]);run.child('golden',2,case,'R_to_P','search');run.child('golden',2,case,'R_to_P','score')
 raise SystemExit(0)
from worker import execute
sources={}
for name in ['engine','optimized']:
 for p in (ROOT/name).rglob('*'):
  if p.is_file() and p.suffix in ['.py','.cpp','.h','.c','.so'] and not any(x in p.parts for x in ['build','__pycache__']):sources[str(p.relative_to(ROOT))]=run.sha(p)
for rel,digest in run.read(ROOT/'engine.json')['sources'].items():
 if not rel.endswith('.so'):assert run.sha(ROOT/'engine'/rel)==digest,rel
for rel,digest in run.read(ROOT/'optimized-engine.json')['changed_files'].items():
 assert run.sha(ROOT/'optimized'/rel)==digest,rel
run.save(ROOT/'hpc/linux-engine.json',dict(platform=platform.platform(),python=sys.version,sources=sources,source_commit='3a9ef9a8aec3e129ea5f4ee67df47ce2c25ade3c',mac_manifest_sha256=run.sha(ROOT/'engine.json'),optimized_manifest_sha256=run.sha(ROOT/'optimized-engine.json')))
rows=[]
fields=['reference_recovery','representative_recovery','candidate_terminals','capped','best_target_heavy_coverage','best_target_all_atom_coverage']
for case in [0,1]:
 result=execute([sys.executable,__file__,'parity',str(case)],ROOT/'hpc'/f'parity-{case}.log',memory=6144)
 assert result['status']=='passed',result
 a=run.read(ROOT/'runs/golden/seed2'/f'case{case}'/'R_to_P/evaluation.json')
 b=run.read(ROOT/'hpc/parity/runs/golden/seed2'/f'case{case}'/'R_to_P/evaluation.json')
 differences={k:[a.get(k),b.get(k)] for k in fields if a.get(k)!=b.get(k)}
 assert not differences,differences
 rows.append(dict(case=case,status='passed',fields={k:b.get(k) for k in fields},execution=result))
run.save(ROOT/'hpc/build-proof.json',dict(status='passed',source_hash_verification=True,parity=rows,checkpoint_test=run.read(ROOT/'checkpoint-verification-tests.json'),scope='Native Linux build, focused tests, unchanged frozen source hashes, two Mac/Linux search-and-reference parity checks. Timing remains platform-specific.'))
print('Linux build and parity verified.',flush=True)
