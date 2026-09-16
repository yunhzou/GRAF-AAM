"""Assemble final manuscript evidence from the completed mixed-host campaign."""
import json,gzip,hashlib,shutil,collections,statistics,tarfile
from pathlib import Path
V=Path(__file__).resolve().parent;R=V/'hpc-return/final';M=V/'paper-stage';E=M/'evidence';OUT=V/'report-publish'
cache={}
with tarfile.open(R/'results.tar.gz','r:gz') as archive:
 for member in archive:
  if member.isfile() and member.name.endswith('.json'):
   cache[member.name]=archive.extractfile(member).read()
print('Indexed archived metadata:',len(cache),flush=True)
def cached(p):
 try:return cache.get(str(p.relative_to(R)))
 except ValueError:return None
def exists(p):
 try:return str(p.relative_to(R)) in cache
 except ValueError:return p.exists()
def read(p):
 data=cached(p)
 return json.loads(data if data is not None else p.read_bytes())
def sha(p):
 data=cached(p)
 return hashlib.sha256(data if data is not None else p.read_bytes()).hexdigest()
def save(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2)+'\n')
def gz(p,d):
 with gzip.open(p,'wt') as f:json.dump(d,f,separators=(',',':'))
summary=read(R/'hpc/result-summary.json');assert summary['completed']==summary['tasks']==6706 and not summary['missing_task_indices']
assert summary['golden']['10']['counts']==dict(recovered=1840,not_recovered=9,unknown=2)
timing=read(V/'hpc/same-mac-timing.json');common=sorted(set(range(1851))-set(timing['excluded_cases']));assert len(common)==1821
OUT.mkdir(exist_ok=True)
methods={};direction_records=[]
for seed in [1,2,3,10]:
 rows=[]
 for item in summary['golden'][str(seed)]['rows']:
  c=item['case'];cpu=0;complete=True;hosts=set()
  for direction in ['R_to_P','P_to_R']:
   f=R/'runs/golden'/f'seed{seed}'/f'case{c}'/direction
   files={n:read(f/n) for n in ['search.json','evaluation.json','search_execution.json','score_execution.json','checkpoint-verification.json'] if exists(f/n)}
   se=files.get('search_execution.json',{});complete &= se.get('status')=='passed';cpu+=files.get('search.json',{}).get('cpu_seconds',0)
   hosts.add('Linux' if se.get('host') or se.get('platform')=='linux' else 'Mac')
   direction_records.append(dict(seed=seed,case=c,direction=direction,files=files,metadata_sha256={n:sha(f/n) for n in files}))
  rows.append(dict(item,search_complete=complete,search_cpu_including_io=cpu,search_hosts=sorted(hosts)))
 outcomes={k:summary['golden'][str(seed)]['counts'].get(k,0) for k in ['recovered','not_recovered','unknown']}
 cost=statistics.mean(rows[c]['search_cpu_including_io'] for c in common) if seed in [1,2] else None
 if cost is not None:
  assert all(rows[c]['search_complete'] and rows[c]['search_hosts']==['Mac'] for c in common)
  assert abs(cost-timing['metrics'][f'aam_seed{seed}']['mean_seconds'])<1e-8
 methods[f'seeds{seed}']=dict(golden_cases=1851,golden_outcomes=outcomes,golden_recovery_percent=100*outcomes['recovered']/1851,common_mean_cpu_seconds=cost,per_case=rows,incomplete_search_cases=[r['case'] for r in rows if not r['search_complete']])
 print('Assembled Golden seed',seed,flush=True)
seed=dict(algorithm_commit='3a9ef9a8aec3e129ea5f4ee67df47ce2c25ade3c',methods=methods,common_case_indices=common,fresh=True,denominator=1851,same_host_timing=timing,timing_scope=timing['scope'],configuration='Both directions; uncut and single-edge sweeps; iso 1; floor 0.2; cap 100; root seed 42; 1/2/3/10 orderings; no competition.')
save(E/'seed_comparison.json',seed)
slaprows=summary['slap']['rows'];counts={k:summary['slap']['counts'].get(k,0) for k in ['recovered','not_recovered','unknown']};hits={r['case'] for r in slaprows if r['outcome']=='recovered'}
slap=dict(sweep_union_recovered=counts['recovered'],outcomes=counts,max_recovery_if_unknowns_succeed=counts['recovered']+counts['unknown'],fresh=True,cases=1851,per_case=slaprows,upstream_commit='ea248fd9494f52f4865193e87a98cc92c62b5f9e',configuration='Both directions, binary/weighted modes, uncut and every single-edge deletion; original initialization and explicit H.',paired_recovery={})
for k,m in methods.items():
 h={r['case'] for r in m['per_case'] if r['outcome']=='recovered'}
 slap['paired_recovery'][k]=dict(both=len(h&hits),aam_only=sorted(h-hits),slap_only=sorted(hits-h),conservative_count_advantage=len(h)-counts['recovered']-counts['unknown'],neither_verified=sorted(set(range(1851))-(h|hits)))
save(E/'slap_sweep.json',slap)
import finalize_evidence as fe
oldread=fe.read
fe.read=lambda p:read(R/'coordinate-audit/summary.json') if p==V/'coordinate-audit/summary.json' else oldread(p)
comp,flat=fe.coordinate()
scopes={'seed_comparison.json':'Final Golden reference verification, including positive witnesses in saved cuts; separate same-Mac timing cohort','slap_sweep.json':'Fresh pinned SLAP Golden sweep and reference verification','competition_final.json':'One-seed coordinate search, competition and complete event-window decoding; fresh comparator witnesses','final_dedup.json':'Fresh final representation audit and complete event-window decoding on all 140 cases'}
save(E/'paper_sources.json',dict(algorithm_commit=seed['algorithm_commit'],frozen_search_commit=read(V/'engine.json')['commit'],version_proof='reports/current_validation_20260912/version-proof.json',fresh_campaign_complete=True,snapshots=[dict(snapshot=n,source='reports/current_validation_20260912/'+n,scope=scope,sha256=sha(E/n)) for n,scope in scopes.items()],limits='Coordinate final pipeline has only a one-seed run. Golden ten-seed has two unresolved reactions. Timing across Mac and Linux is not pooled.'))
for n in [*scopes,'paper_sources.json']:shutil.copy2(E/n,OUT/n)
gz(OUT/'golden-direction-records.json.gz',direction_records);del direction_records
gz(OUT/'slap-golden-evaluations.json.gz',[read(R/'slap-golden/evaluations'/f'{c}.json') for c in range(1851)])
gz(OUT/'coordinate-decoded.json.gz',[read(V/'decoded-optimized'/f'case{c}'/'result.json') for c in range(140)])
gz(OUT/'coordinate-archive-identities.json.gz',[dict(case=c,**read(V/'decoded-optimized'/f'case{c}'/'identity.json')) for c in range(140)])
gz(OUT/'coordinate-slap-comparison.json.gz',read(V/'fresh-coordinate-comparison.json'))
gz(OUT/'hpc-task-executions.json.gz',[read(p) for p in sorted((R/'hpc/tasks').glob('*.json'))])
gz(OUT/'late-cut-verification.json.gz',[dict(folder=p.parent.name,execution=read(p),proof=read(p.parent/'new-proof.json') if exists(p.parent/'new-proof.json') else None) for p in sorted((R/'hpc/late-cuts').glob('*/execution.json'))])
for n in ['engine.json','optimized-engine.json','version-proof.json','environment.json','inputs.json','completed-event-comparison.json','slap-case90-tie-audit.json','optimized-tests.log']:
 shutil.copy2(V/n,OUT/n)
for src,name in [('golden-data/manifest.json','dataset.json'),('slap-golden/manifest.json','slap-golden-protocol.json')]:shutil.copy2(V/src,OUT/name)
for src,name in [('hpc/runtime.json','linux-runtime.json'),('hpc/linux-engine.json','linux-engine.json'),('hpc/build-proof.json','linux-build-proof.json'),('hpc/smoke-proof.json','linux-smoke-proof.json'),('checkpoint-verification-tests.json','checkpoint-verification-tests.json'),('coordinate-audit/summary.json','coordinate-audit.json'),('hpc/submission.json','hpc-submission.json'),('hpc/cpu-allocation-correction.json','cpu-allocation-correction.json'),('hpc/result-summary.json','final-summary.json'),('hpc/first-pass-summary.json','first-pass-summary.json')]:
 shutil.copy2(R/src,OUT/name)
shutil.copy2(V/'hpc/same-mac-timing.json',OUT/'same-mac-timing.json')
drivers=OUT/'drivers';drivers.mkdir(exist_ok=True)
for n in ['run.py','competition.py','decode_optimized.py','slap_golden.py','slap_coordinate.py','audit_coordinate.py','score_checkpoints.py','finalize_hpc_paper.py']:shutil.copy2(V/n,drivers/n)
for n in ['worker.py','aggregate.py','build_proof.py','smoke.py','verify_remaining.py','array.sbatch','build.sbatch','verify_remaining.sbatch','aggregate.sbatch']:
 shutil.copy2(R/'hpc'/n,drivers/n)
(OUT/'README.md').write_text('''# Final-source validation, September 2026

The algorithm source is commit 3a9ef9a. The search/evaluation source is unchanged from eb1a7d1; the final commit adds exact adjacent nested-subgroup absorption during final-family decoding. The source hashes and version proof record this boundary.

Golden includes all 1,851 original records, dataset revision 793475e. Search uses both directions, floor 0.2, iso tolerance 1, cap 100, root seed 42, uncut plus single-edge sweeps, and 1/2/3/10 seed orderings. Competition is disabled. Reference labels enter only verification. Recovered / not recovered / unknown counts are 1834/17/0, 1834/17/0, 1837/14/0, and 1840/9/2. The pinned SLAP sweep returns 1796/52/3. These are reference-family recoveries, not top-1 accuracy. The remaining ten-seed unknowns are cases 590 and 1358.

The coordinate evaluation uses all 140 original WBO inputs, forward one-seed search at cap 2,000, and the 128-completion competition policy. The final two-seed coordinate pipeline was not measured. Complete decoding covers all 236,653 saved flat families within fixed event windows, producing 336 reaction-specific classes. Coverage of returned SLAP-sweep minimum-event classes is 166/168, with every compared class recovered in 139/140 reactions. Native SLAP coverage is 155/160. This measures comparator-output coverage, not chemical accuracy. XYZ-derived comparator adjacency is restored and verified; XYZ preparation itself is not rerun. SLAP hydrogen label families are not exhaustively decoded.

Final decoding takes 25.70 recorded CPU-minutes and 11.60 elapsed minutes on three Mac workers, including continued passes and certificate journals, excluding preceding search and competition. The 1,821-case same-Mac timing cohort gives mean search/workflow CPU seconds of 2.518 (one seed), 4.953 (two seeds), and 4.819 (SLAP sweep). Reference verification is excluded. AAM search includes checkpoint IO; SLAP graph/mapping/export instrumentation excludes IO. Three- and ten-seed mixed-host times are not used for speedup claims.

The remaining campaign ran on NRT cpu_short. All 6,706 work items were processed; the initial direction-level statuses include 72 memory limits and 9 timeouts. These are not counts of unresolved reactions because an opposite direction can verify recovery. Separate positive-only checks of saved cuts recovered ten of twelve remaining unknown reactions, without new search. Only verified positive membership is credited from partial archives; absence is never inferred from an interrupted search. Each subprocess retains a five-minute watchdog; continuations use at most 6 GiB sampled RSS.

NRT allocates two CPUs per requested single-CPU worker. The initial 300-worker release therefore temporarily allocated 600 CPUs. Concurrency was corrected to 150 workers / 300 allocated CPUs, and 142 excess active workers were requeued with completed task results and raw cut checkpoints retained. The correction record is included. The main release-to-summary duration was 42.96 minutes; Slurm records about 47.40 CPU-hours for the main array. These totals include multiple configurations, checks and retries, and must not be compared to one configuration's successful search CPU.

Compressed records retain configurations, archive hashes, explicit verification witnesses, decoder patterns, initial resource limits, and subsequent cut proofs. All 140 coordinate audits and full-window decoder checks passed. Linux checks comprise 52 decoder tests, eight native tests, two Mac/Linux search/reference parity cases, and end-to-end AAM/SLAP/audit checks. The original Mac archives and large HPC raw checkpoints remain in their campaign directories, identified in the evidence; they are not duplicated in Git. The exact drivers are included, with their original import layout and run paths. Relocation requires adapting paths and regenerating manifests in a new run directory; do not run the Slurm scripts without inspecting their account and project paths.

The manuscript contains the final algorithm and measured results. Operational details and development audits belong in this report. The coordinate data's original WBO provenance and redistribution terms still require author confirmation before publishing the underlying dataset.
''')
save(OUT/'SHA256SUMS.json',{str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='SHA256SUMS.json'})
print('Final paper snapshots and report prepared:',OUT,flush=True)
