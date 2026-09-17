import os,sys,time,json,gzip,dataclasses,subprocess,concurrent.futures,argparse
from pathlib import Path
from itertools import combinations
import numpy as np
REPO=Path(os.environ.get('GRAFT_AUDIT_REPO', Path(__file__).resolve().parents[3]))
sys.path.insert(0,str(REPO/'src'))
from graft.artifacts import read_aam_checkpoint
from graft.final_branches import FinalBranchCatalogue
from graft.event_patterns import SignedEventIndex
from graft.postprocessing import DecodedEvents,EventCandidate,EventDecodeConfig
from graft.chirality import ChiralityConfig,select_chiral_witness,query_chiral_witness
from graft.family_query import query_path
ROOT=Path(os.environ.get('GRAFT_AUDIT_WORKSPACE', Path.cwd()))
OUT=ROOT/'outputs/chirality-strict-audit-20260917'

def sign(xyz,center,shell,tol):
    # Independent determinant calculation: same mathematical normalization,
    # without the production measurement/parity helpers or solver predicates.
    pts=np.asarray(xyz[list(shell)],dtype=np.longdouble)
    if len(shell)==3:
        a,b,c=pts-np.asarray(xyz[center],dtype=np.longdouble)
        scale=np.sqrt(np.dot(a,a)*np.dot(b,b)*np.dot(c,c))
    else:
        a,b,c=pts[1:]-pts[0]
        scale=np.prod([np.dot(pts[i]-pts[j],pts[i]-pts[j]) for i,j in combinations(range(4),2)])**np.longdouble(.25)
    if scale==0:return 0
    pos=[a[0]*b[1]*c[2],a[1]*b[2]*c[0],a[2]*b[0]*c[1]]
    neg=[a[2]*b[1]*c[0],a[1]*b[0]*c[2],a[0]*b[2]*c[1]]
    det=sum(pos)-sum(neg);eps=np.finfo(np.longdouble).eps
    # Conservative roundoff guard, including both affine and center-relative
    # edge scales. Near-degenerate frames are also reported as undefined.
    base=np.sqrt(np.dot(a,a)*np.dot(b,b)*np.dot(c,c))
    bound=16*eps/(1-16*eps)*max(sum(abs(v) for v in pos+neg),base)
    if abs(det)<=bound or abs(det/scale)<=tol:return 0
    return 1 if det>0 else -1

def geometry(problem,mapping):
    r,p=problem.reactant,problem.product
    nr=[tuple(int(a) for a in np.flatnonzero(row>=.2)) for row in r.wbo]
    np_=[set(int(a) for a in np.flatnonzero(row>=.2)) for row in p.wbo]
    rows=[]
    for center,neighbors in enumerate(nr):
        persistent=tuple(a for a in neighbors if mapping[a] in np_[mapping[center]])
        frames=[]
        if len(persistent) in (3,4):frames.append(('ordinary',persistent,.1))
        if len(neighbors)>4 or len(np_[mapping[center]])>4:
            for size in ((3,4) if len(neighbors)>4 else (3,)):
                frames.extend(('high_coordinate',shell,0 if size==3 else .1) for shell in combinations(persistent,size))
        for kind,shell,tol in frames:
            sr=sign(r.coordinates,center,shell,tol);sp=sign(p.coordinates,mapping[center],tuple(mapping[a] for a in shell),tol)
            status='undefined' if not(sr and sp) else 'preserved' if sr==sp else 'violation'
            rows.append(dict(kind=kind,center=center,neighbors=shell,source_sign=sr,target_sign=sp,status=status))
    return dict(violations=[r for r in rows if r['status']=='violation'],defined=sum(r['status']!='undefined' for r in rows),undefined=sum(r['status']=='undefined' for r in rows))

def dump(path,record):
    path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(record,indent=2));tmp.replace(path)

def one(i):
    start=time.monotonic();folder=ROOT/f'outputs/final-end-to-end-optimized-20260916/case{i}';out=OUT/f'case{i}.json'
    saved=json.load(open(folder/'candidates.json'));settings=json.load(open(folder/'result.json'))
    aam=read_aam_checkpoint(folder/'cuts/aam.pkl.gz');cat=FinalBranchCatalogue(aam.problem).add_aam(aam);problem=aam.problem;del aam
    candidates=[EventCandidate(c['id'],{int(a):b for a,b in c['mapping'].items()},{k:tuple(tuple(b) for b in v) for k,v in c['events'].items()},c['total'],c['family_ids']) for c in saved['candidates']]
    dec=DecodedEvents(candidates,cat,SignedEventIndex(problem,threshold=settings['decode_config']['threshold'],metal_threshold=settings['decode_config']['metal_threshold']),saved['family_reports'],EventDecodeConfig(**settings['decode_config']))
    old=json.load(open(ROOT/f'work/chirality-comparison/old-{problem.name}.json'))
    record=dict(case=i,name=problem.name,families=len(cat.families),status='running',rows=[],old=[])
    for m in old['mechanisms']:
        mapping={int(a):b for a,b in m['mapping_RP'].items()};geo=geometry(problem,mapping)
        events=dec.index.describe([mapping[a] for a in range(dec.index.n)])['events'];candidate=next((c for c in candidates if c.events==events),None)
        o=dict(id=m['id'],geometry=geo,matching_candidate_id=None if candidate is None else candidate.id)
        if candidate is not None and not geo['violations']:
            q=query_chiral_witness(dec,candidate,mapping,ChiralityConfig(seconds=60))
            o['strict_membership']=dataclasses.asdict(q)
            assert q.status in ('allowed','forbidden','unknown')
        record['old'].append(o)
    minimum=min(c.total for c in candidates);minimum_candidates=[c for c in candidates if c.total==minimum]
    record['expected_candidates']=len(minimum_candidates);dump(out,record)
    for ci,c in enumerate(minimum_candidates):
        wall,cpu=time.perf_counter(),time.process_time();res=select_chiral_witness(dec,c,ChiralityConfig(seconds=240))
        row=dict(candidate=ci,id=c.id,total_events=c.total,result=dataclasses.asdict(res),cpu_seconds=time.process_time()-cpu,wall_seconds=time.perf_counter()-wall,original_geometry=geometry(problem,c.mapping))
        record['rows'].append(row);dump(out,record)
        if res.status=='allowed':
            row['selected_geometry']=geometry(problem,res.mapping)
            assert not row['selected_geometry']['violations'],row
            assert dec.index.describe([res.mapping[a] for a in range(dec.index.n)])['events']==c.events
            status,witness=query_path(cat.families[res.family_id].as_path(problem),problem,res.mapping,source_atoms=tuple(range(dec.index.n)),timeout_ms=60000)
            assert status=='recovered' and dict(witness['mapping'])==res.mapping
            row['independent_membership_and_events_verified']=True
        elif res.status=='forbidden':
            assert res.diagnostics['search_exhaustive']
        print(i,ci,res.status,round(row['cpu_seconds'],3),flush=True);dump(out,record)
    record.update(status='complete',wall_seconds=time.monotonic()-start);dump(out,record)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--case',type=int);ap.add_argument('--cases',nargs='*',type=int);args=ap.parse_args()
    if args.case is not None:one(args.case)
    else:
        OUT.mkdir(exist_ok=True)
        def launch(i):
            existing=OUT/f'case{i}.json'
            if existing.exists() and json.loads(existing.read_text()).get('status')=='complete':return i,'cached'
            with open(OUT/f'case{i}.log','w') as log:
                try:
                    p=subprocess.run([sys.executable,__file__,'--case',str(i)],stdout=log,stderr=subprocess.STDOUT,timeout=300)
                    status='complete' if p.returncode==0 else 'error'
                except subprocess.TimeoutExpired:status='watchdog'
            dump(OUT/f'execution{i}.json',dict(case=i,status=status));return i,status
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
            futures=[ex.submit(launch,i) for i in (args.cases if args.cases else range(140))]
            for f in concurrent.futures.as_completed(futures):print(*f.result(),flush=True)
