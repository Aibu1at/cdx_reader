"""Real binary CDX regressions with source provenance and stereo-aware checks."""
import hashlib
import json
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from cdx_reader import read_cdx_molecule
from cdx_reader._normalization import R_LABEL
from .r7_reference import legacy_comparison_copy

DATA = Path(__file__).parent / 'data' / 'patent_regressions'
REFERENCES = json.loads((DATA / 'references.json').read_text())

# These frozen chemical references were originally also controls for the old
# implicit selection from a reaction scheme. Selection now belongs to the agent:
# retain exactly the same expected molecule, addressed by its original object ID.
SELECTED_SOURCE_OBJECTS = {
    'US20210079003A1_C00690': '20109',
    'US20210017189A1_C00004': '302',
}


def canonical(mol):
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True)


def test_corpus_diversity_and_provenance():
    assert len(REFERENCES) >= 60
    assert len({r['patent'] for r in REFERENCES}) >= 20
    assert len({r['sha256'] for r in REFERENCES}) == len(REFERENCES)
    assert all(r['member'] and r['evidence']['details'] for r in REFERENCES)


@pytest.mark.parametrize('reference', REFERENCES, ids=lambda r: r['id'])
@pytest.mark.parametrize('input_mode', ['path', 'bytes'])
def test_real_patent_regression(reference, input_mode):
    path = DATA / (reference['id'] + '.cdx')
    data = path.read_bytes()
    assert hashlib.sha256(data).hexdigest() == reference['sha256']
    assert hashlib.sha256(path.with_suffix('.mol').read_bytes()).hexdigest() == reference['mol_sha256']
    expected = reference['expected']
    source = path if input_mode == 'path' else data
    if 'error' in expected:
        with pytest.raises(Exception, match=expected['error']):
            read_cdx_molecule(source)
        return
    if reference['id'] in SELECTED_SOURCE_OBJECTS:
        identifier = SELECTED_SOURCE_OBJECTS[reference['id']]
        mol = next(m for m in read_cdx_molecule(source, as_single=False)
                   if m is not None and m.GetProp('CDXML_FRAG_ID') == identifier)
    else:
        mol = read_cdx_molecule(source)
    if expected.get('none'):
        assert mol is None
        return
    assert mol is not None
    if 'smiles' in expected:
        # Six source-reviewed wavy bars were mistaken for carbons in these
        # frozen references. Check their actual new ports, then compare only
        # the untouched remainder in the old representation (no oracle rewrite).
        comparison = legacy_comparison_copy(mol, reference)
        assert canonical(comparison) == expected['smiles']
        assert Chem.MolToInchiKey(comparison) == expected['inchikey']
        assert rdMolDescriptors.CalcMolFormula(comparison) == expected['formula']
    if 'raw_smiles' in expected:
        clean = Chem.Mol(mol)
        for atom in clean.GetAtoms():
            atom.SetAtomMapNum(0)
            if (atom.HasProp('_CDX_ATOM_ID') and atom.GetProp('_CDX_ATOM_ID')
                    in expected.get('legacy_comparison_carbon_ids', [])):
                # Keep the old deferred-chemistry reference unchanged while
                # permitting only the newly verified R-label substitutions.
                # Their actual element, label and neighbors are checked below.
                atom.SetAtomicNum(6)
        clean.UpdatePropertyCache(strict=False)
        assert Chem.MolToSmiles(clean) == expected['raw_smiles']
        assert Chem.MolToInchiKey(clean) == expected['inchikey']
        assert rdMolDescriptors.CalcMolFormula(clean) == expected['formula']
    if 'variable_labels' in expected:
        atoms = {a.GetProp('_CDX_ATOM_ID'): a for a in mol.GetAtoms() if a.HasProp('_CDX_ATOM_ID')}
        for identifier, label in expected['variable_labels'].items():
            atom = atoms[identifier]
            assert atom.GetAtomicNum() == 0
            assert atom.GetProp('dummyLabel') == label
            assert atom.GetProp('_CDX_ATOM_LABEL').replace(' ', '') == label
            assert atom.GetAtomMapNum() == 0
            if 'variable_source_neighbors' in expected:
                assert sorted(neighbor.GetProp('_CDX_ATOM_ID') for neighbor in atom.GetNeighbors()) == expected['variable_source_neighbors'][identifier]
                assert not atom.HasProp('_MolFileRLabel')
    if 'source_aromatic_bonds' in expected:
        assert mol.GetBoolProp('_CDX_UNKEKULIZED_QUERY')
        atoms = {a.GetProp('_CDX_ATOM_ID'): a.GetIdx() for a in mol.GetAtoms()
                 if a.HasProp('_CDX_ATOM_ID')}
        for bond in expected['source_aromatic_bonds']:
            assert mol.GetBondBetweenAtoms(atoms[bond['begin']], atoms[bond['end']]).GetBondType() == Chem.BondType.AROMATIC
    if 'repaired_hydrogens' in expected:
        atoms = {a.GetProp('_CDX_ATOM_ID'): a for a in mol.GetAtoms() if a.HasProp('_CDX_ATOM_ID')}
        for identifier, count in expected['repaired_hydrogens'].items():
            assert atoms[identifier].GetNumExplicitHs() == count
    if 'tetrahedral_source_centers' in expected:
        atoms = {a.GetProp('_CDX_ATOM_ID'): a for a in mol.GetAtoms()
                 if a.HasProp('_CDX_ATOM_ID')}
        for identifier, (neighbors, expected_parity) in expected['tetrahedral_source_centers'].items():
            atom = atoms[identifier]
            ids = [neighbor.GetProp('_CDX_ATOM_ID') for neighbor in atom.GetNeighbors()]
            assert sorted(ids) == neighbors
            tag = atom.GetChiralTag()
            assert tag in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            inversions = sum(ids[i] > ids[j] for i in range(len(ids)) for j in range(i + 1, len(ids)))
            parity = (1 if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW else -1) * (-1) ** inversions
            assert parity == expected_parity
    if 'r_labels' in expected:
        found = {a.GetProp('_CDX_ATOM_ID'): a.GetProp('dummyLabel') for a in mol.GetAtoms()
                 if a.HasProp('dummyLabel') and R_LABEL.fullmatch(a.GetProp('dummyLabel'))}
        assert found == expected['r_labels']
        for atom in mol.GetAtoms():
            if not atom.HasProp('dummyLabel') or atom.GetProp('_CDX_ATOM_ID') not in expected['r_labels']:
                continue
            assert atom.GetAtomicNum() == 0
            assert atom.GetAtomMapNum() == 0
            number = atom.GetProp('dummyLabel')[1:]
            if number.isdigit():
                assert atom.GetIntProp('_MolFileRLabel') == int(number)
            else:
                assert not atom.HasProp('_MolFileRLabel')
        # Check MOL R-number exchange on copies of the verified R atoms. Some
        # complete query drawings contain other labels/bonds that MOL cannot
        # serialize; that is independent of numeric R-label preservation.
        r_atoms = Chem.RWMol()
        for atom in mol.GetAtoms():
            if atom.HasProp('_MolFileRLabel'):
                r_atoms.AddAtom(Chem.Atom(atom))
        reread = Chem.MolFromMolBlock(Chem.MolToMolBlock(r_atoms, forceV3000=True), sanitize=False)
        assert sorted(a.GetIntProp('_MolFileRLabel') for a in reread.GetAtoms()
                      if a.HasProp('_MolFileRLabel')) == sorted(int(v[1:]) for v in expected['r_labels'].values() if v[1:].isdigit())


@pytest.mark.parametrize('reference', [r for r in REFERENCES if 'r_labels' in r['expected']],
                         ids=lambda r: r['id'])
def test_r_labels_in_existing_list_mode(reference):
    mols = read_cdx_molecule(DATA / (reference['id'] + '.cdx'), as_single=False)
    found = {a.GetProp('_CDX_ATOM_ID'): a.GetProp('dummyLabel')
             for mol in mols if mol is not None for a in mol.GetAtoms()
             if a.HasProp('dummyLabel') and R_LABEL.fullmatch(a.GetProp('dummyLabel'))}
    assert found == reference['expected']['r_labels']
