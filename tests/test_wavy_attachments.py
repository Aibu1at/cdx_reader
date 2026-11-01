"""Original binary patent drawings: attachment bars versus other wavy bonds."""
import hashlib
import io
import json
from pathlib import Path

import pytest
from lxml import etree
from pycdxml import cdxml_converter
from rdkit import Chem

from cdx_reader import read_cdx_molecule
from cdx_reader._attachments import normalize_wavy_attachment_points
from .r7_reference import REVIEWED
from .crossed_reference import prior_comparison_copy

DATA = Path(__file__).parent / 'data/wavy_attachments'


def test_fixtures_are_unmodified_source_files():
    for record in json.loads((DATA/'manifest.json').read_text()):
        assert record['member'].endswith('.CDX')
        assert hashlib.sha256((DATA/record['file']).read_bytes()).hexdigest() == record['sha256']


@pytest.mark.parametrize('ref,expected,centre,tails,neighbor', [
    ('C00005','*Cc1ccccc1B(O)O','31',['33','35'],'29'),
    ('C00006','*Cc1cc(OC(F)(F)F)ccc1B(O)O','31',['33','35'],'29'),
    ('C00013','*CCCCB(O)O','11',['17','19'],'13'),
    ('C00014','*Cc1ccccc1','31',['33','35'],'29'),
    ('C00106','*c1ccc([N+](=O)[O-])cc1','29',['31','33'],'11'),
])
@pytest.mark.parametrize('mode',['bytes','path','objects'])
def test_bar_becomes_an_anonymous_source_addressed_port(ref,expected,centre,tails,neighbor,mode):
    path = DATA/('US20140206647A1_'+ref+'.cdx')
    source = path.read_bytes() if mode=='bytes' else path
    mol = read_cdx_molecule(source,as_single=mode!='objects')
    if mode=='objects':
        assert len(mol)==1
        mol=mol[0]
    assert Chem.MolToSmiles(mol)==expected
    atoms={a.GetProp('_CDX_ATOM_ID'):a for a in mol.GetAtoms() if a.HasProp('_CDX_ATOM_ID')}
    port=atoms[centre]
    assert port.GetAtomicNum()==0 and port.GetDegree()==1
    assert port.GetProp('_CDX_ATTACHMENT_REPAIR')=='wavy_two_tail_bar'
    assert json.loads(port.GetProp('_CDX_ATTACHMENT_TAIL_IDS'))==tails
    assert port.GetNeighbors()[0].GetProp('_CDX_ATOM_ID')==neighbor
    assert not port.GetAtomMapNum() and not port.HasProp('_MolFileRLabel')
    assert not port.HasProp('dummyLabel') and not any(t in atoms for t in tails)
    assert not sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())


@pytest.mark.parametrize('ref',['C00008','C00029'])
def test_one_tail_and_internal_wavy_bonds_are_untouched(ref):
    doc=cdxml_converter.read_cdx(io.BytesIO((DATA/('US20100119475A1_'+ref+'.cdx')).read_bytes()))
    before=etree.tostring(doc.cdxml)
    assert any(b.get('Display')=='Wavy' for b in doc.cdxml.getroot().iter('b'))
    assert normalize_wavy_attachment_points(doc.cdxml)=={}
    assert etree.tostring(doc.cdxml)==before


def test_mixed_wavy_drawing_keeps_other_bonds_and_is_idempotent():
    path=DATA/'US20120035134A1_C00006.cdx'
    doc=cdxml_converter.read_cdx(io.BytesIO(path.read_bytes()))
    bonds={b.get('id'):dict(b.attrib) for b in doc.cdxml.getroot().iter('b')}
    repairs=normalize_wavy_attachment_points(doc.cdxml)
    assert repairs
    tails={t for ts in repairs.values() for t in ts}
    after={b.get('id'):dict(b.attrib) for b in doc.cdxml.getroot().iter('b')}
    assert after=={i:b for i,b in bonds.items() if b.get('B') not in tails and b.get('E') not in tails}
    # The original TIF explicitly equates these four bars with -NR2/-CH2NR2/
    # -OR/-CH2OR. Two have crossed stems; the four formerly retained waves are
    # their marker tails, not bonds conveying molecular stereochemistry.
    assert repairs == {'11': ['26', '28'], '56': ['58', '60', '62'],
                       '77': ['81', '82'], '107': ['108', '109', '110']}
    assert sum(b.get('Display')=='Wavy' for b in after.values())==0
    once=etree.tostring(doc.cdxml)
    assert normalize_wavy_attachment_points(doc.cdxml)=={}
    assert etree.tostring(doc.cdxml)==once


@pytest.mark.parametrize('identifier', sorted(REVIEWED))
def test_source_reviewed_legacy_bars_are_not_complete_molecules(identifier):
    review=REVIEWED[identifier]
    path=DATA.parent/'patent_regressions'/(identifier+'.cdx')
    assert hashlib.sha256(path.read_bytes()).hexdigest()==review['sha256']
    mol=read_cdx_molecule(path)
    # Check any newly reviewed crossed bar first; compare the earlier R7 edit
    # in its historical representation without overwriting that reference.
    mol=prior_comparison_copy(mol, {'id': identifier, 'sha256': review['sha256']})
    assert Chem.MolToSmiles(mol)==review['expected_repaired_smiles']
    ports=[a for a in mol.GetAtoms() if a.HasProp('_CDX_ATTACHMENT_REPAIR')]
    assert {a.GetProp('_CDX_ATOM_ID') for a in ports}==set(review['bars'])
    assert all(a.GetAtomicNum()==0 and a.GetDegree()==1 for a in ports)
