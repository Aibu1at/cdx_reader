"""Keep the historical references immutable while checking a source-reviewed R7 edit.

Six old CDX/MOL-consensus references encode bars as isopropyl groups. Their original
TIFs and patent text contradict that interpretation (see R7_REFERENCE_REVIEW.md).
Assert the new ports first, then undo ONLY those source-addressed edits on a COPY
for the old graph comparison. This checks the untouched remainder, not chemical
correctness of that remainder. test_wavy_attachments checks the actual new output.
"""
import json
from pathlib import Path

from rdkit import Chem
from .crossed_reference import prior_comparison_copy

REVIEWED = json.loads((Path(__file__).parent / 'data/wavy_attachments/reviewed_legacy_cases.json').read_text())


def legacy_comparison_copy(mol, reference):
    mol = prior_comparison_copy(mol, reference)
    review = REVIEWED.get(reference['id'])
    if review is None:
        return mol
    assert reference['sha256'] == review['sha256']
    assert Chem.MolToSmiles(mol) == review['expected_repaired_smiles']
    rw = Chem.RWMol(mol)
    atoms = {a.GetProp('_CDX_ATOM_ID'): a for a in rw.GetAtoms() if a.HasProp('_CDX_ATOM_ID')}
    for centre, tails in review['bars'].items():
        atom = atoms[centre]
        assert atom.GetAtomicNum() == 0 and atom.GetDegree() == 1
        assert atom.GetProp('_CDX_ATTACHMENT_REPAIR') == 'wavy_two_tail_bar'
        assert json.loads(atom.GetProp('_CDX_ATTACHMENT_TAIL_IDS')) == tails
        assert all(tail not in atoms for tail in tails)
        atom.SetAtomicNum(6)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
        centre_index = atom.GetIdx()
        for tail in tails:
            source = review['source_nodes'][tail]
            assert source.get('Element', '6') == '6'
            new_atom = Chem.Atom(6)
            new_atom.SetProp('_CDX_ATOM_ID', tail)
            rw.AddBond(centre_index, rw.AddAtom(new_atom), Chem.BondType.SINGLE)
    rw.UpdatePropertyCache(strict=False)
    return rw.GetMol()
