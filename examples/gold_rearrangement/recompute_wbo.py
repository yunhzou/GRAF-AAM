"""Optional: regenerate the input bond weights with tblite==0.7.0 (GFN2-xTB).

Original supplied endpoint geometries are kept fixed. Requires numpy, RDKit,
and tblite. Set OMP_NUM_THREADS=1 for this small single-point calculation.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from rdkit import Chem
from tblite.interface import Calculator


def compute(xyz_path):
    rows = [s.split() for s in xyz_path.read_text().splitlines()[2:] if s.strip()]
    elements = [r[0] for r in rows]
    coordinates = np.array([[float(x) for x in r[1:4]] for r in rows])
    numbers = np.array([Chem.GetPeriodicTable().GetAtomicNumber(e) for e in elements])
    calculator = Calculator('GFN2-xTB', numbers, coordinates / .529177210903, charge=1, uhf=0)
    calculator.set('verbosity', 0)
    result = calculator.singlepoint()
    wbo = np.asarray(result.get('bond-orders')).squeeze()
    assert wbo.shape == (len(elements), len(elements))
    wbo = (wbo + wbo.T) / 2
    preprocessing = dict(symmetrized=True, raw_minimum=float(wbo.min()),
                         negative_values_clipped_to_zero=int((wbo < 0).sum()))
    wbo = np.maximum(wbo, 0)
    np.fill_diagonal(wbo, 0)
    return dict(elements=elements, coordinates=coordinates.tolist(), wbo=wbo.tolist(),
                source=xyz_path.name, source_sha256=hashlib.sha256(xyz_path.read_bytes()).hexdigest(),
                method='GFN2-xTB single point, tblite 0.7.0, supplied DFT geometry, charge +1, closed shell',
                energy_hartree=float(result.get('energy')), preprocessing=preprocessing)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('xyz', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    args.output.write_text(json.dumps(compute(args.xyz), indent=2) + '\n')
