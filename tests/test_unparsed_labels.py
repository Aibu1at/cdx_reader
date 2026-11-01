"""Real-deposit regressions for the boundary between chemistry and undecoded text."""
import hashlib
import json
from pathlib import Path
import pytest
from rdkit import Chem
from cdx_reader import read_cdx_molecule

DATA = Path(__file__).with_name('data') / 'unparsed_labels'


def fixture(name):
    path = DATA / name
    record = next(r for r in json.loads((DATA / 'manifest.json').read_text()) if r['file'] == name)
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record['sha256']
    return path


@pytest.mark.parametrize('as_bytes', [False, True])
def test_het_keeps_source_identity_but_is_not_a_declared_port(as_bytes):
    path = fixture('US20090221600A1_C00002.cdx')
    mol = read_cdx_molecule(path.read_bytes() if as_bytes else path)
    atom = next(a for a in mol.GetAtoms() if a.GetProp('_CDX_ATOM_ID') == '41')
    assert atom.GetAtomicNum() == 0
    assert atom.GetProp('_CDX_UNPARSED_LABEL') == 'het'
    assert atom.GetProp('_CDX_ATOM_LABEL') == 'het'


def test_s_annotation_is_a_typed_empty_object_not_methane():
    path = fixture('US20200352925A1_C00005.cdx')
    objects = read_cdx_molecule(path, as_single=False)
    assert [m.GetProp('CDXML_FRAG_ID') for m in objects] == ['73', '101']
    assert objects[1].GetNumAtoms() == 0
    assert objects[1].GetProp('_CDX_OBJECT_KIND') == 'annotation'
    assert json.loads(objects[1].GetProp('_CDX_ANNOTATIONS'))[0]['text'] == '(S)'
    assert Chem.MolToSmiles(read_cdx_molecule(path)) == Chem.MolToSmiles(objects[0])


def test_bullet_hydrochloride_stays_chemical_and_separate():
    objects = read_cdx_molecule(fixture('US20200352925A1_C00006.cdx'), as_single=False)
    chloride = next(m for m in objects if m.GetProp('CDXML_FRAG_ID') == '118')
    assert Chem.MolToSmiles(chloride) == 'Cl'
    assert not any(a.HasProp('_CDX_UNPARSED_LABEL') for a in chloride.GetAtoms())


def test_expansion_keeps_the_parent_object_identifier():
    objects = read_cdx_molecule(fixture('US20210122713A1_C00127.cdx'), as_single=False)
    assert [m.GetProp('CDXML_FRAG_ID') for m in objects] == ['14', '90']
    assert Chem.MolToSmiles(objects[1]) == 'O=C(O)C(F)(F)F'
    assert not any(a.GetAtomicNum() == 0 for a in objects[0].GetAtoms())


def test_opaque_nested_fragment_is_not_an_attachment_site():
    mol = read_cdx_molecule(fixture('US20110135650A1_C00661.cdx'))
    for identifier, label in [('36', 'SN3-1'), ('84', 'SN3-2')]:
        atom = next(a for a in mol.GetAtoms() if a.HasProp('_CDX_ATOM_ID') and a.GetProp('_CDX_ATOM_ID') == identifier)
        assert atom.GetAtomicNum() == 0
        assert atom.GetProp('_CDX_UNPARSED_LABEL') == label


def test_a_usable_nickname_graph_is_preserved():
    mol = read_cdx_molecule(fixture('US20040002022A1_C00026.cdx'))
    assert not any(a.HasProp('_CDX_UNPARSED_LABEL') for a in mol.GetAtoms())
    # Independently retained alkyl chains in this ChemDraw Nickname deposit.
    assert mol.HasSubstructMatch(Chem.MolFromSmiles('CCCCCCCC'))


@pytest.mark.parametrize('filename,label,identifier', [
    ('US20180022705A1_C00028.cdx', 'A', '70'),
    ('US20150152119A1_C00076.cdx', 'X', '47'),
])
def test_existing_declared_variables_are_still_ports(filename, label, identifier):
    mol = read_cdx_molecule(fixture(filename))
    atom = next(a for a in mol.GetAtoms() if a.HasProp('_CDX_ATOM_ID')
                and a.GetProp('_CDX_ATOM_ID') == identifier)
    assert atom.GetAtomicNum() == 0
    assert atom.GetProp('_CDX_ATOM_LABEL') == label
    assert not atom.HasProp('_CDX_UNPARSED_LABEL')
