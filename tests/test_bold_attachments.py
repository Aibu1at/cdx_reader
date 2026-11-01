"""Bold attachment bars from original USPTO CDX, and untouched stereo drawings."""
from hashlib import sha256
from io import BytesIO
import json
from pathlib import Path

from lxml import etree
import pytest
from pycdxml import cdxml_converter
from rdkit import Chem

from cdx_reader import read_cdx_molecule
from cdx_reader._attachments import normalize_bold_attachment_points

DATA = Path(__file__).parent / 'data/bold_attachments'


def test_sources_are_unmodified_deposits():
    for row in json.loads((DATA/'manifest.json').read_text()):
        assert sha256((DATA/row['file']).read_bytes()).hexdigest() == row['sha256']


@pytest.mark.parametrize('cid,expected', [
    ('C00463', '*CCCCC(N*)(C(=O)OC)C(F)(F)F'),
    ('C00466', '*CC/C=C/C(=O)N*'),
])
@pytest.mark.parametrize('objects', [False, True])
def test_bold_bars_are_source_addressed_ports(cid, expected, objects):
    raw = (DATA/f'US20120041190A1_{cid}.cdx').read_bytes()
    mol = read_cdx_molecule(raw, as_single=not objects)
    if objects:
        assert len(mol) == 1
        mol = mol[0]
    assert Chem.MolToSmiles(mol) == Chem.CanonSmiles(expected)
    atoms = {a.GetProp('_CDX_ATOM_ID'): a for a in mol.GetAtoms() if a.HasProp('_CDX_ATOM_ID')}
    for aid, tails, neighbour in [('11', ['30', '32'], '13'), ('38', ['40', '42'], '27')]:
        a = atoms[aid]
        assert a.GetAtomicNum() == 0 and a.GetDegree() == 1
        assert a.GetProp('_CDX_ATTACHMENT_REPAIR') == 'bold_two_tail_bar'
        assert json.loads(a.GetProp('_CDX_ATTACHMENT_TAIL_IDS')) == tails
        assert a.GetNeighbors()[0].GetProp('_CDX_ATOM_ID') == neighbour
        assert not any(t in atoms for t in tails)
    xml = cdxml_converter.read_cdx(BytesIO(raw)).cdxml
    assert set(normalize_bold_attachment_points(xml)) == {'11', '38'}
    assert normalize_bold_attachment_points(xml) == {}


@pytest.mark.parametrize('cid', ['C00452', 'C00453', 'C00460'])
def test_real_stereo_and_core_drawings_are_untouched(cid):
    xml = cdxml_converter.read_cdx(BytesIO((DATA/f'US20120041190A1_{cid}.cdx').read_bytes())).cdxml
    before = etree.tostring(xml)
    assert normalize_bold_attachment_points(xml) == {}
    assert etree.tostring(xml) == before
