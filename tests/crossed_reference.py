"""Check source-reviewed crossed ports, then reconstruct the preceding test view.

Frozen CDX/MOL-consensus answers are left intact. The comparison copy tests only
the untouched remainder; test_crossed_attachments checks the actual corrected graph.
"""
import json
from pathlib import Path

from rdkit import Chem

REVIEWED = json.loads((Path(__file__).parent / 'data/crossed_attachments/reviewed_legacy_cases.json').read_text())


def prior_comparison_copy(mol, reference):
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
        assert atom.GetProp('_CDX_ATTACHMENT_REPAIR') == 'wavy_crossed_stem_bar'
        assert json.loads(atom.GetProp('_CDX_ATTACHMENT_TAIL_IDS')) == tails
        assert len(tails) == 3 and all(t not in atoms for t in tails)
        atom.SetAtomicNum(6)
        atom.SetNumExplicitHs(0)
        atom.SetNoImplicit(False)
        atom.ClearProp('_CDX_ATTACHMENT_REPAIR')
        atom.ClearProp('_CDX_ATTACHMENT_TAIL_IDS')
        for tail in tails:
            source = review['source_nodes'][tail]
            assert source.get('Element', '6') == '6'
            endpoint = Chem.Atom(6)
            endpoint.SetProp('_CDX_ATOM_ID', tail)
            rw.AddBond(atom.GetIdx(), rw.AddAtom(endpoint), Chem.BondType.SINGLE)
    rw.UpdatePropertyCache(strict=False)
    return rw.GetMol()
