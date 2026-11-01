"""Original CDX crossed stems: three patents, two full products and old controls."""
import hashlib
from io import BytesIO
import json
from pathlib import Path

from lxml import etree
from pycdxml import cdxml_converter
import pytest
from rdkit import Chem

from cdx_reader import read_cdx_molecule
from cdx_reader._attachments import normalize_wavy_attachment_points
from .crossed_reference import REVIEWED

DATA = Path(__file__).parent / 'data/crossed_attachments'


def test_fixtures_are_original_binary_sources():
    for entry in json.loads((DATA / 'manifest.json').read_text()):
        assert hashlib.sha256((DATA / entry['file']).read_bytes()).hexdigest() == entry['sha256']
    for entry in json.loads((DATA / 'artwork_manifest.json').read_text()):
        assert hashlib.sha256((DATA / entry['file']).read_bytes()).hexdigest() == entry['sha256']


@pytest.mark.parametrize('identifier,smiles,centre,tails,neighbor', [
    ('US20150152119A1_C00066', '*c1cccnc1', '45', ['47', '48', '46'], '41'),
    ('US20140227302A1_C00035', '*CCn1cccn1', '45', ['49', '51', '55'], '43'),
    ('US20160060265A1_C00068', '*c1ccc(C(=O)C2CC2)cc1', '29', ['33', '35', '31'], '21'),
])
@pytest.mark.parametrize('mode', ['bytes', 'path', 'objects'])
def test_crossed_stem_is_a_source_addressed_anonymous_port(identifier,smiles,centre,tails,neighbor,mode):
    path = DATA / (identifier + '.cdx')
    mol = read_cdx_molecule(path.read_bytes() if mode == 'bytes' else path, as_single=mode != 'objects')
    if mode == 'objects':
        assert len(mol) == 1
        mol = mol[0]
    assert Chem.MolToSmiles(mol) == smiles
    atoms = {a.GetProp('_CDX_ATOM_ID'): a for a in mol.GetAtoms() if a.HasProp('_CDX_ATOM_ID')}
    port = atoms[centre]
    assert port.GetAtomicNum() == 0 and port.GetDegree() == 1
    assert port.GetNeighbors()[0].GetProp('_CDX_ATOM_ID') == neighbor
    assert port.GetProp('_CDX_ATTACHMENT_REPAIR') == 'wavy_crossed_stem_bar'
    assert json.loads(port.GetProp('_CDX_ATTACHMENT_TAIL_IDS')) == tails
    assert not port.HasProp('dummyLabel') and not port.HasProp('_MolFileRLabel') and not port.GetAtomMapNum()
    assert all(t not in atoms for t in tails)
    assert not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    raw = cdxml_converter.read_cdx(BytesIO(path.read_bytes())).cdxml
    normalize_wavy_attachment_points(raw)
    once = etree.tostring(raw)
    assert normalize_wavy_attachment_points(raw) == {}
    assert etree.tostring(raw) == once


def test_a_bar_with_two_bare_stem_ends_is_ambiguous_and_untouched():
    # This real drawing has four separate crossed markers intersecting ring bonds.
    # Neither ordinary end belongs to the molecular graph: do not pick a side.
    path = DATA.parent / 'wavy_attachments/US20100119475A1_C00008.cdx'
    raw = cdxml_converter.read_cdx(BytesIO(path.read_bytes())).cdxml
    before = etree.tostring(raw)
    assert normalize_wavy_attachment_points(raw) == {}
    assert etree.tostring(raw) == before


@pytest.mark.parametrize('identifier', sorted(REVIEWED))
def test_source_reviewed_old_references_get_the_actual_corrected_graph(identifier):
    review = REVIEWED[identifier]
    path = DATA.parent / 'patent_regressions' / (identifier + '.cdx')
    assert hashlib.sha256(path.read_bytes()).hexdigest() == review['sha256']
    mol = read_cdx_molecule(path)
    assert Chem.MolToSmiles(mol) == review['expected_repaired_smiles']
    ports = [a for a in mol.GetAtoms() if a.HasProp('_CDX_ATTACHMENT_REPAIR')
             and a.GetProp('_CDX_ATTACHMENT_REPAIR') == 'wavy_crossed_stem_bar']
    assert {a.GetProp('_CDX_ATOM_ID') for a in ports} == set(review['bars'])
