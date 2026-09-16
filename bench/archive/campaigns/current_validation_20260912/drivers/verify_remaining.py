"""Complete-coverage verification of existing cuts for unresolved reactions."""
import sys,os,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'hpc'))
import worker
rows=worker.read(ROOT/'hpc/remaining-cut-tasks.json');task=rows[int(sys.argv[1])]
seed,case,direction=task['seed'],task['case'],task['direction'];folder=ROOT/'runs/golden'/f'seed{seed}'/f'case{case}'/direction
out=ROOT/'hpc/late-cuts'/f'seed{seed}-case{case}-{direction}';out.mkdir(parents=True,exist_ok=True)
for name in ['evaluation.json','score_execution.json','checkpoint-verification.json']:
 if (folder/name).exists() and not (out/name).exists():shutil.copy2(folder/name,out/name)
row=worker.execute([sys.executable,str(ROOT/'score_checkpoints.py'),str(seed),str(case),direction],out/'run.log',memory=6144)
row.update(phase='checkpoint_score',task=task,previous_attempt=str(out),scope='Reference membership in saved cuts; absence requires the complete expected schedule.')
proof=folder/'checkpoint-verification.json'
if proof.exists():shutil.copy2(proof,out/'new-proof.json')
verdict=worker.read(proof).get('reference_recovery','unknown') if row['status']=='passed' and proof.exists() else 'unknown'
resolved=verdict in ('recovered','not_recovered')
if resolved:
 worker.save(folder/'evaluation.json',worker.read(proof));worker.save(folder/'score_execution.json',row)
row['reference_recovery']=verdict;worker.save(out/'execution.json',row)
print(row,flush=True)
