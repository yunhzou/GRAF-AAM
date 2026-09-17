"""Current-family chirality contract; tiny exhaustive sets are test oracles only."""
from itertools import permutations
import numpy as np
import pytest
from graft import AAMProblem, MolecularEndpoint, AAMSearchConfig, search_aam
from graft.final_branches import FinalBranchCatalogue, FinalFamily
from graft.postprocessing import decode_events, EventDecodeConfig
from graft.chirality import ChiralityConfig, select_chiral_witness
from graft.alignment.index_chirality import _simplex_measure


def tetra():
    xyz = np.array([[0,0,0],[-.9,-.9,-.9],[.9,.9,-.9],[.9,-.9,.9],[-.9,.9,.9]])
    w = np.zeros((5,5));w[0,1:]=1;w[1:,0]=1
    return ('C','F','H','H','H'),xyz,w


def decoded_family(elements, r, p, wr, wp=None, *, actions=(), mapping=None, extra=(), threshold=.5):
    n=len(elements);mapping=mapping or dict(enumerate(range(n)))
    problem=AAMProblem(MolecularEndpoint(elements,r,wr),MolecularEndpoint(elements,p,wr if wp is None else wp))
    catalogue=FinalBranchCatalogue(problem)
    required=tuple((a,b) for a in range(n) for b in range(a+1,n) if wr[a,b]>=.2 and problem.product.wbo[mapping[a],mapping[b]]>=.2)
    catalogue.families=[FinalFamily(tuple(sorted(mapping.items())),required,tuple(actions),(.2,1.0)),*extra]
    result=decode_events(catalogue,EventDecodeConfig(threshold=threshold,seconds_per_family=30))
    assert result.complete
    return result


def swapped(x, a=3,b=4):
    p=x.copy();p[[a,b]]=p[[b,a]];return p


def test_current_public_search_decode_then_chirality():
    e,x,w=tetra();problem=AAMProblem(MolecularEndpoint(e,x,w),MolecularEndpoint(e,swapped(x),w))
    aam=search_aam(problem,AAMSearchConfig(sweep_cuts=False,branch_limit=200),workers=1)
    raw=aam.graph.to_record();decoded=decode_events(aam);candidate=decoded.minimum_candidates[0]
    before=dict(candidate.mapping)
    result=select_chiral_witness(decoded,candidate)
    assert result.status=='allowed'
    assert decoded.query(candidate,result.mapping).status=='allowed'
    assert result.diagnostics['atom_bijections_enumerated']==0
    assert not any(f['status']=='violation' for f in result.diagnostics['ordinary_frames'])
    assert raw==aam.graph.to_record() and before==candidate.mapping


def test_repair_pool_with_one_parity_constraint(monkeypatch):
    e,x,w=tetra();decoded=decoded_family(e,x,swapped(x),w,actions=[('pool',(2,3,4))])
    # Force a known valid-family, wrong-handed representative; decoder is free
    # to choose another witness, but the public selector must repair this one.
    c=decoded.candidates[0];c.mapping=dict(enumerate(range(5)))
    import graft.alignment.index_chirality as legacy
    monkeypatch.setattr(legacy,'_generated_atom_permutations',lambda *a:pytest.fail('expanded actions'))
    selected=select_chiral_witness(decoded,c)
    assert selected.status=='allowed' and selected.mapping!=c.mapping
    assert _simplex_measure(swapped(x),0,tuple(selected.mapping[a] for a in range(1,5)),.1).sign==_simplex_measure(x,0,(1,2,3,4),.1).sign
    assert decoded.query(c,selected.mapping).status=='allowed'
    assert selected.diagnostics['compiled_families']==1
    assert selected.diagnostics['solver_checks']<=4


def test_no_chiral_work_for_already_valid_mapping():
    e,x,w=tetra();d=decoded_family(e,x,x,w,actions=[('pool',(2,3,4))]);c=d.candidates[0]
    c.mapping=dict(enumerate(range(5)))
    result=select_chiral_witness(d,c)
    assert result.status=='allowed' and result.mapping==c.mapping
    assert result.diagnostics['solver_checks']==0
    assert result.fixed_mapping_rmsd<1e-12


def test_fixed_inversion_is_reported_and_strict_mode_rejects():
    e,x,w=tetra();d=decoded_family(e,x,swapped(x),w);c=d.candidates[0]
    result=select_chiral_witness(d,c)
    assert result.status=='allowed'
    assert result.diagnostics['ordinary_frames'][0]['status']=='fixed_orientation_change'
    assert select_chiral_witness(d,c,ChiralityConfig(mode='all')).status=='forbidden'


def test_nontrivial_correlated_action_cannot_be_split():
    e,x,w=tetra();elements=e+e;xyz=np.vstack([x,x+[5,0,0]])
    weights=np.zeros((10,10));weights[:5,:5]=w;weights[5:,5:]=w
    perm=list(range(10));perm[3],perm[4]=4,3;perm[8],perm[9]=9,8
    d=decoded_family(elements,xyz,swapped(xyz),weights,actions=[('group',(tuple(perm),))])
    c=d.candidates[0];c.mapping=dict(enumerate(range(10)))
    assert select_chiral_witness(d,c).status=='forbidden'
    assert select_chiral_witness(d,c,ChiralityConfig(mode='all')).status=='forbidden'


@pytest.mark.parametrize('coordination',[4,5])
def test_persistent_shell_survives_coordination_change(coordination):
    e,x,w=tetra();e=e+('Cl',);x=np.vstack([x,[0,0,2]])
    wr=np.zeros((6,6));wp=wr.copy();wr[0,1:coordination+1]=1;wr[1:coordination+1,0]=1
    wp[0,1:10-coordination]=1;wp[1:10-coordination,0]=1
    d=decoded_family(e,x,swapped(x),wr,wp,actions=[('pool',(2,3,4))]);c=d.minimum_candidates[0];c.mapping=dict(enumerate(range(6)))
    selected=select_chiral_witness(d,c)
    assert selected.status=='allowed'
    assert selected.diagnostics['ordinary_frames'][0]['neighbors']==(1,2,3,4)
    assert selected.diagnostics['ordinary_frames'][0]['status']=='preserved'


def test_planar_endpoint_is_undefined_not_wrong_handed():
    e,x,w=tetra();p=x.copy();p[:,2]=0
    d=decoded_family(e,x,p,w,actions=[('pool',(2,3,4))])
    result=select_chiral_witness(d,d.candidates[0],ChiralityConfig(mode='all'))
    assert result.status=='allowed' and result.diagnostics['ordinary_frames'][0]['status']=='undefined'


def test_exact_event_edges_cannot_change_to_repair_orientation():
    e,x,w=tetra();p=swapped(x);wp=w.copy();wp[0,3]=wp[3,0]=.4
    d=decoded_family(e,x,p,w,wp,actions=[('pool',(3,4))]);c=d.candidates[0]
    c.mapping=dict(enumerate(range(5)))
    # Each permitted swap moves the broken edge to the other H. This has the
    # same event class but is outside the chosen concrete event/core.
    strict=select_chiral_witness(d,c,ChiralityConfig(mode='all'))
    assert strict.status=='forbidden'
    assert select_chiral_witness(d,c).diagnostics['ordinary_frames'][0]['status']=='fixed_orientation_change'


def test_union_queries_families_not_recorded_as_candidate_support():
    e,x,w=tetra();perm=tuple([0,1,2,4,3])
    extra=FinalFamily(tuple(enumerate(perm)),tuple((0,a) for a in range(1,5)),(),(.2,1))
    d=decoded_family(e,x,swapped(x),w,extra=[extra]);c=d.candidates[0];c.mapping=dict(enumerate(range(5)));c.family_ids=[0]
    result=select_chiral_witness(d,c)
    assert result.status=='allowed' and result.family_id==1
    assert result.mapping==dict(enumerate(perm))


def test_high_coordinate_reconfiguration_is_explicit():
    rng=np.random.default_rng(51);e=('Sc',)+('O',)*6;x=np.vstack([np.zeros(3),rng.normal(size=(6,3))]);w=np.zeros((7,7));w[0,1:]=.5;w[1:,0]=.5
    p=x.copy();p[:,0]*=-1
    d=decoded_family(e,x,p,w);c=d.candidates[0]
    result=select_chiral_witness(d,c)
    assert result.status=='allowed' and result.diagnostics['reconfigured_high_coordinate_frames']
    assert select_chiral_witness(d,c,ChiralityConfig(high_coordinate='strict')).status=='forbidden'


def test_watchdog_and_validation():
    e,x,w=tetra();d=decoded_family(e,x,x,w)
    assert select_chiral_witness(d,d.candidates[0],ChiralityConfig(seconds=0)).status=='unknown'
    with pytest.raises(ValueError):ChiralityConfig(seconds=-1)
    with pytest.raises(ValueError):ChiralityConfig(graph_floor=0)
    with pytest.raises(ValueError):select_chiral_witness(d,decode_events(d.catalogue).candidates[0])


def test_ordered_overlapping_actions_against_tiny_exhaustive_oracle():
    e,x,w=tetra();actions=[('pool',(2,3)),('pool',(3,4))]
    # Compose exactly two pools, not the closure S3 generated by their union.
    maps=[]
    for first in [(2,3),(3,2)]:
        for second in [(3,4),(4,3)]:
            a=dict(zip((2,3),first));b=dict(zip((3,4),second))
            maps.append({i:b.get(a.get(i,i),a.get(i,i)) for i in range(5)})
    for order in permutations((2,3,4)):
        p=x.copy();p[[2,3,4]]=x[list(order)]
        d=decoded_family(e,x,p,w,actions=actions);c=d.candidates[0]
        result=select_chiral_witness(d,c,ChiralityConfig(mode='all'))
        sign=_simplex_measure(x,0,(1,2,3,4),.1).sign
        valid=[m for m in maps if _simplex_measure(p,m[0],tuple(m[a] for a in (1,2,3,4)),.1).sign==sign]
        assert result.status==('allowed' if valid else 'forbidden')
        assert result.mapping in valid


@pytest.mark.parametrize("mode", ["all", "mutable"])
def test_hard_constraints_agree_with_event_filtered_exhaustive_oracle(mode):
    # Correlated ordered pool actions, multiple event classes and unequal WBOs.
    e,x,w=tetra();w[0,2]=w[2,0]=.4;w[0,3]=w[3,0]=.8
    maps=[]
    for a in [(2,3),(3,2)]:
        for b in [(3,4),(4,3)]:
            first=dict(zip((2,3),a));second=dict(zip((3,4),b))
            maps.append({i:second.get(first.get(i,i),first.get(i,i)) for i in range(5)})
    for order in permutations((2,3,4)):
        p=x.copy();p[[2,3,4]]=x[list(order)]
        d=decoded_family(e,x,p,w,actions=[('pool',(2,3)),('pool',(3,4))],threshold=.3)
        for candidate in d.candidates:
            result=select_chiral_witness(d,candidate,ChiralityConfig(mode=mode))
            sign=_simplex_measure(x,0,(1,2,3,4),.1).sign
            good=[m for m in maps if d.index.describe([m[i] for i in range(5)])['events']==candidate.events
                  and _simplex_measure(p,0,tuple(m[a] for a in (1,2,3,4)),.1).sign==sign]
            same_event=[m for m in maps if d.index.describe([m[i] for i in range(5)])['events']==candidate.events]
            if mode=='mutable' and len(same_event)==1:good=same_event
            assert result.status==('allowed' if good else 'forbidden')
            if good:assert result.mapping in good


def test_rotation_translation_and_atom_relabeling_invariance():
    e,x,w=tetra();p=swapped(x)
    q,_=np.linalg.qr(np.random.default_rng(13).normal(size=(3,3)))
    if np.linalg.det(q)<0:q[:,0]*=-1
    for order in [(0,1,2,3,4),(3,1,4,0,2)]:
        inverse={old:new for new,old in enumerate(order)}
        d=decoded_family(tuple(e[a] for a in order),x[list(order)],(p@q+3)[list(order)],w[np.ix_(order,order)],
                         actions=[('pool',tuple(inverse[i] for i in (2,3,4)))])
        result=select_chiral_witness(d,d.candidates[0],ChiralityConfig(mode='all'))
        assert result.status=='allowed'
        assert all(f['status']=='preserved' for f in result.diagnostics['ordinary_frames'])


def test_public_search_anchors_survive_chirality_repair():
    e,x,w=tetra();problem=AAMProblem(MolecularEndpoint(e,x,w),MolecularEndpoint(e,swapped(x),w))
    aam=search_aam(problem,AAMSearchConfig(sweep_cuts=False,branch_limit=200,anchors=((3,3),)),workers=1)
    d=decode_events(aam)
    result=select_chiral_witness(d,d.candidates[0],ChiralityConfig(mode='all'))
    assert result.status=='allowed' and result.mapping[3]==3
    assert d.query(d.candidates[0],result.mapping).status=='allowed'


def test_high_coordinate_retained_basis_is_satisfied_and_maximal():
    rng=np.random.default_rng(51);e=('Sc',)+('O',)*6;x=np.vstack([np.zeros(3),rng.normal(size=(6,3))]);w=np.zeros((7,7));w[0,1:]=.5;w[1:,0]=.5
    p=x.copy();p[:,0]*=-1
    g=list(range(7));g[1],g[2]=2,1
    d=decoded_family(e,x,p,w,actions=[('group',(tuple(g),))]);c=d.candidates[0]
    result=select_chiral_witness(d,c)
    assert result.status=='allowed'
    from graft.chirality import _Workspace
    workspace=_Workspace(d,c,ChiralityConfig())
    retained=result.diagnostics['retained_high_coordinate_frames']
    assert not workspace.soft_violations(result.mapping,retained)
    assert result.diagnostics['reconfigured_high_coordinate_frames']
    for frame in result.diagnostics['reconfigured_high_coordinate_frames']:
        # Only two possible atom actions in this independent test oracle.
        assert all(workspace.soft_violations(m,[*retained,frame]) for m in [dict(enumerate(range(7))),dict(enumerate(g))])


def test_solver_unknown_is_not_a_chirality_conflict(monkeypatch):
    e,x,w=tetra();d=decoded_family(e,x,swapped(x),w,actions=[('pool',(2,3,4))]);c=d.candidates[0];c.mapping=dict(enumerate(range(5)))
    import graft.chirality as module
    def unknown(*args,**kwargs):raise module._BudgetExpired('simulated solver unknown')
    monkeypatch.setattr(module._Workspace,'check',unknown)
    result=select_chiral_witness(d,c)
    assert result.status=='unknown' and result.mapping is None
