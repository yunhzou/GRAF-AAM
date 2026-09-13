"""Bundle only material used by the final preprint, with a file-hash manifest."""
from pathlib import Path
import hashlib,json,zipfile
MAN=Path(__file__).resolve().parents[1]
names=['manuscript.pdf','preprint.tex','preprint.bbl','references.bib','README.md','EDITORIAL_NOTES.md','requirements-build.txt','Makefile','natbib.sty']
names += ['assets/'+n for n in ['matterlab.cls','paper-preamble.tex','plainnat.bst','Optimistic.ttf','logo_matterlab.pdf','logo_ac.pdf','logo_uoft.pdf','logo_vector.png','nvidia-logo-vert.png']]
names += ['includes/'+n for n in ['include-abstract.tex','include-body.tex','include-appendix.tex','paper.tex','supplement.tex','generated-seed-table.tex','generated-decoder-table.tex','generated-results.tex']]
names += ['figs/'+stem+'.'+ext for stem in ['fig1_algorithm','fig2_golden','fig3_coordinate'] for ext in ['pdf','svg','png']]
names += ['evidence/'+n for n in ['paper_sources.json','seed_comparison.json','slap_sweep.json','competition_final.json','final_dedup.json','PAPER_EVIDENCE.md']]
names += ['scripts/'+n for n in ['build.sh','build_figures.py','validate_outputs.py','package_bundle.py']]
checks={n:hashlib.sha256((MAN/n).read_bytes()).hexdigest() for n in sorted(names)}
with zipfile.ZipFile(MAN/'manuscript_bundle.zip','w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for n in names:z.write(MAN/n,'manuscript/'+n)
 z.writestr('manuscript/SHA256SUMS.json',json.dumps(checks,indent=2)+'\n')
 for n in ['artifact-validation.json','visual-review.json']:
  if (MAN/'build'/n).is_file():
   record=json.loads((MAN/'build'/n).read_text())
   if record.get('manuscript_sha256')==checks['manuscript.pdf']:z.write(MAN/'build'/n,'manuscript/validation/'+n)
with zipfile.ZipFile(MAN/'manuscript_bundle.zip') as z:assert z.testzip() is None
print(f'Bundle: {len(names)} files, {(MAN/"manuscript_bundle.zip").stat().st_size/1e6:.2f} MB')
