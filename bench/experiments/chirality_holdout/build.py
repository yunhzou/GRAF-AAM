"""Package the existing 3D interpolation viewers into the original offline navigator."""
from pathlib import Path
import argparse,gzip,hashlib,importlib.util,json,shutil


def read(path):return json.loads(Path(path).read_text())
def dump(path,value):
    Path(path).parent.mkdir(parents=True,exist_ok=True);Path(path).write_text(json.dumps(value,indent=2)+'\n')


def build(repo,source,output):
    source=source.resolve();output=output.resolve();output.mkdir(parents=True,exist_ok=True)
    rows=sorted([read(f) for f in source.glob('case*/summary.json')],key=lambda x:x['case'])
    complete=[r for r in rows if r.get('status')=='complete'];batch=source/'viewer-batch'
    specs={k:[] for k in ['small','medium','large']}
    for row in complete:
        paired_rows=[r for r in row['rows'] if 'after' in r]
        trend=(' [worse]' if any(r['after']['max_clashes']>r['before']['max_clashes'] for r in paired_rows)
               else ' [improved]' if any(r['after']['max_clashes']<r['before']['max_clashes'] for r in paired_rows) else '')
        step=f"{row['case']:03d} · {row['name']}"+trend;i=row['case'];tier='small' if row['atoms']<50 else 'medium' if row['atoms']<100 else 'large'
        specs[tier].append(dict(source_index=i,step_id=step,atom_count=row['atoms']))
        with gzip.open(source/f'case{i}/viewer-data.json.gz','rt') as f:doc=json.load(f)
        dump(batch/'cases'/step/'summary.json',row)
        dump(batch/'cases'/step/'rp_stage.json',dict(mechanisms=doc['mechanisms']))
        v=batch/'views'/step/'view.html';v.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(source/f'case{i}/view.html',v)
    for tier,cases in specs.items():dump(batch/'manifests'/f'{tier}.json',dict(cases=cases))
    dump(batch/'inventory.json',dict(case_count=len(complete)))
    incomplete=[r['case'] for r in rows if r.get('status')!='complete']
    dump(batch/'batch_summary.json',dict(error_count=len(incomplete),total_elapsed_case_seconds=sum(r['elapsed_seconds'] for r in complete)))
    spec=importlib.util.spec_from_file_location('integrated_viewer',repo/'tools/build_integrated_batch_viewer.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    result=mod.build(batch,output/'viewer.html',compact=True)
    html=(output/'viewer.html').read_text().replace('AAM integrated viewer','Chirality and interpolation audit').replace(' mechanisms',' mapping views').replace(' mechanism${',' mapping view${')
    html=html.replace('violations=${c.violations}', 'paired before/after audit').replace('Mechanism:', 'Candidate / mapping:')
    (output/'viewer.html').write_text(html)
    all_candidates=[dict(case=r['case'],**c) for r in rows for c in r['rows']]
    paired=[r for r in all_candidates if 'after' in r]
    def summarize(items):
        return dict(candidates=len(items),changed=sum(c.get('changed_atoms',0)>0 for c in items),
           clear_before=sum(c['before']['max_clashes']==0 for c in items),clear_after=sum(c['after']['max_clashes']==0 for c in items),
           improved=sum(c['after']['max_clashes']<c['before']['max_clashes'] for c in items),
           worsened=sum(c['after']['max_clashes']>c['before']['max_clashes'] for c in items),
           unchanged_peak=sum(c['after']['max_clashes']==c['before']['max_clashes'] for c in items),
           mean_chirality_cpu_seconds=sum(c['chirality']['cpu_seconds'] for c in items)/max(1,len(items)),
           total_chirality_cpu_seconds=sum(c['chirality']['cpu_seconds'] for c in items),
           total_interpolation_cpu_seconds=sum(c['interpolation']['cpu_seconds'] for c in items))
    report=dict(expected_cases=140,complete_cases=len(complete),incomplete_cases=incomplete,
                processed_candidates=len(all_candidates),paired=summarize(paired),minimum=summarize([c for c in paired if c['minimum']]),
                unresolved=[dict(case=c['case'],candidate=c['candidate'],status=c['status'],reason=c['selected']['diagnostics'].get('reason')) for c in all_candidates if c['status']!='allowed'],
                regressions=[dict(case=c['case'],candidate=c['candidate'],minimum=c['minimum'],before=c['before'],after=c['after']) for c in paired if c['after']['max_clashes']>c['before']['max_clashes']],
                interpolation='existing graft.alignment.interpolation.internal_coordinate_interpolation, 101 frames',
                clash_threshold=.70,config=dict(chirality_mode='mutable',high_coordinate='maximal',high_coordinate_scope='selected_family',graph_floor=.2,orientation_tolerance=.1,minimum_only=all(r.get('minimum_only') for r in rows)),
                sources={f:hashlib.sha256((repo/f).read_bytes()).hexdigest() for f in ['src/graft/chirality.py','src/graft/alignment/interpolation.py','src/graft/static/reaction_viewer.html']})
    dump(output/'summary.json',report)
    with gzip.open(output/'case-results.json.gz','wt') as f:json.dump(rows,f,separators=(',',':'))
    print(json.dumps({**result,**report},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--repo',type=Path,required=True);p.add_argument('--source',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();build(a.repo,a.source,a.output)
