"""Replay saved holdout candidates through chirality and the existing interpolation viewer.
No search, decoding, symmetry enumeration, or RMSD ranking is performed.
"""
from pathlib import Path
import argparse, concurrent.futures, dataclasses, gzip, hashlib, json, os, subprocess, sys, time, traceback
import numpy as np
from graft.artifacts import read_aam_checkpoint
from graft.final_branches import FinalBranchCatalogue
from graft.event_patterns import SignedEventIndex
from graft.postprocessing import DecodedEvents, EventCandidate, EventDecodeConfig
from graft.chirality import select_chiral_witness, ChiralityConfig
from graft.family_query import query_path
from graft.alignment.interpolation import internal_coordinate_interpolation, proper_align_coordinates
from graft.viewers import comparison_document, reaction_html


def dump(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(data,separators=(',',':'))+'\n');temp.replace(path)


def verify(decoded,candidate,mapping,preferred=None):
    n=decoded.index.n
    assert sorted(mapping)==sorted(mapping.values())==list(range(n))
    assert decoded.index.describe([mapping[i] for i in range(n)])['events']==candidate.events
    order=list(dict.fromkeys(([preferred] if preferred is not None else [])+candidate.family_ids+list(range(len(decoded.catalogue.families)))))
    for fid in order:
        status,_=query_path(decoded.catalogue.families[fid].as_path(decoded.catalogue.problem),decoded.catalogue.problem,
                            mapping,source_atoms=tuple(range(n)),timeout_ms=300000)
        if status=='recovered':return fid
        if status=='unknown':raise TimeoutError('independent family-membership check was inconclusive')
    raise AssertionError('selected mapping is outside every saved family')


def interpolate(problem,mapping):
    r,p=problem.reactant,problem.product;n=len(r.elements)
    product=p.coordinates[[mapping[a] for a in range(n)]]
    rb={(a,b) for a in range(n) for b in range(a+1,n) if r.wbo[a,b]>=.2}
    pb={(a,b) for a in range(n) for b in range(a+1,n) if p.wbo[mapping[a],mapping[b]]>=.2}
    path=internal_coordinate_interpolation(r.coordinates,product,r.elements,bonded_pairs=rb|pb,
         persistent_bonded_pairs=rb&pb,reactant_bonded_pairs=rb,product_bonded_pairs=pb,n_frames=101)
    coords=np.asarray([f['coords'] for f in path['frames']]);assert np.isfinite(coords).all()
    assert np.allclose(coords[0],r.coordinates) and np.allclose(coords[-1],proper_align_coordinates(product,r.coordinates))
    counts=[f['clashes']['count'] for f in path['frames']]
    summary=dict(max_clashes=max(counts),clashing_frames=sum(c>0 for c in counts),sum_frame_clashes=sum(counts),
                 endpoint_clashes=[counts[0],counts[-1]],peak_frame=int(np.argmax(counts)),
                 minimum_radius_ratio=min(f['clashes']['minimum_radius_ratio'] for f in path['frames'] if f['clashes']['minimum_radius_ratio'] is not None))
    return path,summary


def record(problem,candidate,mapping,label):
    r,p=problem.reactant,problem.product
    events=[]
    for kind,bonds in candidate.events.items():
        for a,b in bonds:
            events.append(dict(kind=kind,r=[a,b],p=[mapping[a],mapping[b]],wbo=[float(r.wbo[a,b]),float(p.wbo[mapping[a],mapping[b]])]))
    return dict(name=label,mapping=sorted(mapping.items()),events=events,pattern=candidate.id)


def one(args):
    folder=args.source/f'case{args.case}';out=args.output/f'case{args.case}';out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter();saved=json.loads((folder/'candidates.json').read_text());settings=json.loads((folder/'result.json').read_text())
    aam=read_aam_checkpoint(folder/'cuts/aam.pkl.gz');problem=aam.problem
    cat=FinalBranchCatalogue(problem).add_aam(aam);del aam
    assert len(cat.families)==settings['completed_families']
    candidates=[EventCandidate(c['id'],{int(a):b for a,b in c['mapping'].items()},{k:tuple(tuple(b) for b in v) for k,v in c['events'].items()},c['total'],c['family_ids']) for c in saved['candidates']]
    decoded=DecodedEvents(candidates,cat,SignedEventIndex(problem,threshold=settings['decode_config']['threshold'],metal_threshold=settings['decode_config']['metal_threshold']),saved['family_reports'],EventDecodeConfig(**settings['decode_config']))
    minimum=min(c.total for c in candidates);available=len(candidates)
    if args.minimum_only:candidates=[c for c in candidates if c.total==minimum]
    rows=[];records=[];paths=[];diags=[]
    summary=dict(case=args.case,name=problem.name,atoms=len(problem.reactant.elements),families=len(cat.families),
        available_decoded_candidates=available,minimum_only=args.minimum_only,
        input_candidates=len(candidates),minimum_candidates=sum(c.total==minimum for c in candidates),
        saved_decode_complete=saved['complete'],source_run=str(args.source.resolve()),rows=rows,source_config=settings['config'],decode_config=settings['decode_config'])
    dump(out/'summary.json',summary)
    for ci,c in enumerate(candidates):
        print(f'case {args.case} candidate {ci+1}/{len(candidates)} chirality',flush=True)
        wall,cpu=time.perf_counter(),time.process_time();selected=select_chiral_witness(decoded,c,ChiralityConfig())
        timing=dict(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu)
        row=dict(candidate=ci,id=c.id,total_events=c.total,minimum=c.total==minimum,status=selected.status,
                 chirality=timing,original_mapping=c.mapping,selected=dataclasses.asdict(selected))
        rows.append(row);dump(out/'summary.json',summary)
        original_fid=verify(decoded,c,c.mapping)
        if selected.status=='allowed':
            row['verified_family_id']=verify(decoded,c,selected.mapping,selected.family_id)
            row['membership_verified']=True;row['concrete_events_preserved']=True
            row['changed_atoms']=sum(selected.mapping[a]!=b for a,b in c.mapping.items())
        row['original_verified_family_id']=original_fid
        wall,cpu=time.perf_counter(),time.process_time()
        before,before_stats=interpolate(problem,c.mapping);row['before']=before_stats
        label=f'C{ci+1} · {c.total} events'+(' · minimum' if c.total==minimum else '')
        if selected.status=='allowed' and row['changed_atoms']==0:
            row['after']=before_stats
            records.append(record(problem,c,c.mapping,label+' · unchanged, valid'))
            paths.append(before);diags.append(selected.diagnostics)
        else:
            records.append(record(problem,c,c.mapping,label+' · before chirality'))
            paths.append(before);diags.append(None)
            if selected.status=='allowed':
                after,after_stats=interpolate(problem,selected.mapping);row['after']=after_stats
                records.append(record(problem,c,selected.mapping,label+f" · corrected ({row['changed_atoms']} atoms)"))
                paths.append(after);diags.append(selected.diagnostics)
        row['interpolation']=dict(wall_seconds=time.perf_counter()-wall,cpu_seconds=time.process_time()-cpu)
        dump(out/'summary.json',summary)
    endpoints=[dict(elements=e.elements,coordinates=e.coordinates.tolist(),wbo=e.wbo.tolist()) for e in (problem.reactant,problem.product)]
    allowed=sum(r['status']=='allowed' for r in rows);changed=sum(r.get('changed_atoms',0)>0 for r in rows)
    note=f'{allowed}/{len(rows)} candidates satisfy the chirality policy; {changed} required shuffles. Select before/corrected to compare. All successful mappings passed independent saved-family and exact-event checks. '
    note+='The original 101-frame internal-coordinate interpolation is reused. Magenta marks nonbonded distances below 0.70 × summed covalent radii; endpoint bonds are excluded. This is a geometric diagnostic, not an optimized reaction path. '
    if allowed<len(rows):note+='Unresolved candidates show only their original decoded mapping; they are not chirality-corrected. '
    document=comparison_document(dict(index=args.case,name=problem.name,endpoints=endpoints,records=records,note=note))
    for mechanism,path,diag in zip(document['mechanisms'],paths,diags):
        mechanism['endpoint_interpolation']=path;mechanism['index_chirality']=diag
    # Prefer the first corrected/validated minimum candidate for initial display.
    document['default_mech_id']=next((m['id'] for m in document['mechanisms'] if m['index_chirality'] is not None),1)
    (out/'view.html').write_text(reaction_html(document))
    with gzip.open(out/'viewer-data.json.gz','wt') as f:json.dump(document,f,separators=(',',':'))
    summary.update(status='complete',elapsed_seconds=time.perf_counter()-start,viewer_records=len(records))
    dump(out/'summary.json',summary)
    print(f'DONE case {args.case} {allowed}/{len(rows)} allowed; {changed} changed',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--minimum-only',action='store_true');parser.add_argument('--case',type=int);parser.add_argument('--workers',type=int,default=4);parser.add_argument('--cases',type=int,nargs='+',default=list(range(140)))
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    if args.case is not None:
        one(args);return
    def job(i):
        out=args.output/f'case{i}';out.mkdir(exist_ok=True)
        existing=out/'summary.json'
        if existing.exists():
            previous=json.loads(existing.read_text())
            if (previous.get('status')=='complete'
                    and previous.get('minimum_only')==args.minimum_only
                    and previous.get('source_run')==str(args.source.resolve())):
                return dict(case=i,status='cached')
        start=time.perf_counter()
        with (out/'run.log').open('w') as log:
            try:
                child=subprocess.run([sys.executable,__file__,'--source',str(args.source),'--output',str(args.output),'--case',str(i)]+(['--minimum-only'] if args.minimum_only else []),stdout=log,stderr=subprocess.STDOUT,timeout=300)
                status='finished' if child.returncode==0 else 'error'
            except subprocess.TimeoutExpired:status='watchdog'
        result=dict(case=i,status=status,wall_seconds=time.perf_counter()-start)
        dump(out/'execution.json',result);print(json.dumps(result),flush=True);return result
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for r in pool.map(job,args.cases):
            results.append(r);dump(args.output/'campaign.json',dict(workers=args.workers,minimum_only=args.minimum_only,watchdog_seconds=300,searches_rerun=0,decoding_rerun=0,results=results))

if __name__=='__main__':main()
