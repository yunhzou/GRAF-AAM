"""Small post-processing checks; no AAM searches or complete decoding reruns.

python bench/experiments/chirality_postprocessing/run.py --output /tmp/chirality
Uses one process, a 300-second watchdog per case, and no mapping/count caps.
"""
from pathlib import Path
import argparse,gzip,hashlib,json,platform,signal,statistics,time,subprocess,sys
import numpy as np
import z3  # load optional backend before stage timing
from graft import AAMProblem,MolecularEndpoint
from graft.final_branches import FinalBranchCatalogue,FinalFamily
from graft.event_patterns import SignedEventIndex
from graft.postprocessing import DecodedEvents,EventCandidate,EventDecodeConfig
from graft.chirality import select_chiral_witness
from graft.family_query import query_path
ROOT=Path(__file__).resolve().parents[3]

def ep(v):
    return MolecularEndpoint(tuple(v['elements']),np.asarray(v['coordinates']),np.asarray(v['wbo']))

def selected(catalogue,mappings):
    index=SignedEventIndex(catalogue.problem);candidates=[]
    for mapping in mappings:
        pattern=index.describe([mapping[a] for a in range(index.n)])
        candidates.append(EventCandidate(pattern['id'],mapping,pattern['events'],pattern['total'],[]))
    return DecodedEvents(candidates,catalogue,index,[],EventDecodeConfig())

def membership(decoded,mapping):
    for fid,family in enumerate(decoded.catalogue.families):
        status,_=query_path(family.as_path(decoded.catalogue.problem),decoded.catalogue.problem,mapping,
                            source_atoms=tuple(mapping),timeout_ms=300000)
        if status=='recovered':return fid
        if status=='unknown':raise TimeoutError('membership verification timed out')
    raise AssertionError('witness not in saved families')

def cases():
    for k in (1,8,24):
        base=np.array([[0,0,0],[-.9,-.9,-.9],[.9,.9,-.9],[.9,-.9,.9],[-.9,.9,.9]])
        r=np.vstack([base+[4*i,0,0] for i in range(k)]);p=r.copy();n=len(r);w=np.zeros((n,n));actions=[]
        for i in range(k):
            o=5*i;w[o,o+1:o+5]=1;w[o+1:o+5,o]=1;p[[o+3,o+4]]=p[[o+4,o+3]];actions.append(('pool',(o+2,o+3,o+4)))
        elements=('C','F','H','H','H')*k
        problem=AAMProblem(MolecularEndpoint(elements,r,w),MolecularEndpoint(elements,p,w))
        cat=FinalBranchCatalogue(problem)
        cat.families=[FinalFamily(tuple(enumerate(range(n))),tuple((a,b) for a in range(n) for b in range(a+1,n) if w[a,b]),tuple(actions),(.2,1))]
        d=selected(cat,[dict(enumerate(range(n)))])
        yield f'independent_tetrahedra_{k}',d,d.candidates[0],dict(scope='synthetic independent pool actions',represented_assignments=6**k)
    folder=ROOT/'examples/gold_rearrangement/data'
    problem=AAMProblem(ep(json.loads((folder/'reactant.json').read_text())),ep(json.loads((folder/'product.json').read_text())))
    rows=json.loads((folder/'selected-families.json').read_text());cat=FinalBranchCatalogue(problem)
    cat.families=[FinalFamily.from_record(row['family']) for row in rows.values()]
    d=selected(cat,[dict(row['certificate']['mapping']) for row in rows.values()])
    for key,candidate in zip(rows,d.candidates):
        yield 'gold_65_'+key,d,candidate,dict(scope='two published saved families, not full gold search catalogue')
    raw=json.loads(gzip.decompress((Path(__file__).parent/'verification-135.json.gz').read_bytes()))
    problem=AAMProblem(ep(raw['reactant']),ep(raw['product']))
    cat=FinalBranchCatalogue.from_record(problem,raw['catalogue']);d=selected(cat,[dict(cat.families[0].mapping)])
    yield 'structure_verification_135',d,d.candidates[0],dict(scope='complete saved unswept catalogue')
    data=json.loads((ROOT/'manuscript/animations/graft_research_preview/film-data.json').read_text())
    problem=AAMProblem(ep(data['source_input']['reactant']),ep(data['source_input']['product']))
    cat=FinalBranchCatalogue.from_record(problem,json.loads(gzip.decompress((ROOT/'reports/published_baseline_20260915/film/catalogue.json.gz').read_bytes())))
    d=selected(cat,[{int(a):b for a,b in row['mapping'].items()} for row in data['decoded_candidates']])
    for i,candidate in enumerate(d.candidates):
        yield f'golden_15_candidate_{i}',d,candidate,dict(scope='complete baseline catalogue; two published candidates')

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--repeats',type=int,default=3)
    parser.add_argument('--case',help=argparse.SUPPRESS);args=parser.parse_args()
    if args.repeats<1:parser.error('repeats must be positive')
    args.output.mkdir(parents=True,exist_ok=True)
    if args.case is None:
        names=['independent_tetrahedra_1','independent_tetrahedra_8','independent_tetrahedra_24',
               'gold_65_a','gold_65_b','structure_verification_135','golden_15_candidate_0','golden_15_candidate_1']
        rows=[];report={}
        for name in names:
            child=args.output/name
            try:
                subprocess.run([sys.executable,__file__,'--output',str(child),'--repeats',str(args.repeats),
                                '--case',name],check=True,timeout=300)
                report=json.loads((child/'results.json').read_text());rows.extend(report['results'])
            except subprocess.TimeoutExpired:
                rows.append(dict(case=name,status='watchdog',reason='parent terminated process after 300 seconds'))
        report.update(watchdog='parent subprocess hard timeout',watchdog_seconds=300,results=rows)
        (args.output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
        assert all(row['status']=='passed' for row in rows)
        return
    signal.signal(signal.SIGALRM,lambda *_: (_ for _ in ()).throw(TimeoutError('300-second process watchdog')))
    results=[]
    for name,decoded,candidate,scope in cases():
        if name != args.case:continue
        print('START',name,flush=True);signal.alarm(300)
        try:
            candidate.family_ids=[membership(decoded,candidate.mapping)]
            trials=[]
            for _ in range(args.repeats):
                wall,cpu=time.perf_counter(),time.process_time()
                result=select_chiral_witness(decoded,candidate)
                elapsed,processor=time.perf_counter()-wall,time.process_time()-cpu
                assert result.status=='allowed',result
                fid=membership(decoded,result.mapping)
                actual=decoded.index.describe([result.mapping[i] for i in range(decoded.index.n)])
                assert actual['events']==candidate.events
                trials.append(dict(wall_seconds=elapsed,cpu_seconds=processor,family_id=fid,
                                   changed_atoms=sum(candidate.mapping[a]!=p for a,p in result.mapping.items()),
                                   mapping=result.mapping,diagnostics=result.diagnostics))
            row=dict(case=name,status='passed',atoms=decoded.index.n,families=len(decoded.catalogue.families),
                     median_wall_seconds=statistics.median(t['wall_seconds'] for t in trials),
                     median_cpu_seconds=statistics.median(t['cpu_seconds'] for t in trials),trials=trials,**scope)
        except TimeoutError as exc:
            row=dict(case=name,status='watchdog',reason=str(exc),**scope)
        finally:signal.alarm(0)
        results.append(row);print('DONE',name,row['status'],row.get('median_wall_seconds'),flush=True)
        report=dict(workers=1,watchdog_seconds=300,searches_rerun=0,full_decoding_reruns=0,repeats=args.repeats,
                    python=platform.python_version(),z3=z3.get_version_string(),
                    timing_scope='Chirality call only; module imports, saved-input loading and independent membership rechecks excluded. Every trial starts a fresh chirality workspace.',
                    solver_sha256=hashlib.sha256((ROOT/'src/graft/chirality.py').read_bytes()).hexdigest(),results=results)
        (args.output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    assert all(row['status']=='passed' for row in results)

if __name__=='__main__':main()
