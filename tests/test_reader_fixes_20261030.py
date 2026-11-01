"""Graph conservation and explicit stereo from unchanged, source-reviewed patent CDX."""
import copy
import hashlib
import importlib
import json
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from cdx_reader import read_cdx_molecule
from cdx_reader._io import normalize_integer_widths, read_document
from cdx_reader._normalization import normalize_hollow_wedges, repair_invalid_nitrogen_hydrogens

DATA = Path(__file__).with_name('data') / 'reader_fixes_20261030'
MANIFEST = {r['file']: r for r in json.loads((DATA / 'manifest.json').read_text())}


def fixture(name):
    path = DATA / name
    assert hashlib.sha256(path.read_bytes()).hexdigest() == MANIFEST[name]['sha256']
    return path


def objects(name):
    return {m.GetProp('CDXML_FRAG_ID'): m for m in read_cdx_molecule(fixture(name), as_single=False)
            if m is not None}


def source_atom(mol, identifier):
    return next(a for a in mol.GetAtoms()
                if a.HasProp('_CDX_ATOM_ID') and a.GetProp('_CDX_ATOM_ID') == identifier)


@pytest.mark.parametrize('mode', ['path', 'bytes', 'objects'])
def test_child_warning_restores_the_drawn_urea_hydrogens(mode):
    path = fixture('US20080125424A1_C00020.cdx')
    result = read_cdx_molecule(path.read_bytes() if mode == 'bytes' else path,
                               as_single=mode != 'objects')
    mol = result[0] if mode == 'objects' else result
    assert Chem.MolToInchiKey(mol) == 'OBNCRHQEPUHSQY-UHFFFAOYSA-N'
    assert rdMolDescriptors.CalcMolFormula(mol) == 'C28H26N2O'
    assert all(source_atom(mol, identifier).GetTotalNumHs() == 1 for identifier in ['54', '60'])
    assert not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())


def test_child_warning_allows_the_lactam_to_aromatize_without_changing_its_graph():
    mol = objects('US20070123494A1_C00047.cdx')['14']
    assert Chem.MolToSmiles(mol) == 'Nc1ccc2ncc(=O)[nH]c2n1'
    assert Chem.MolToInchiKey(mol) == 'ZBAKGJJWHNVTTA-UHFFFAOYSA-N'
    assert not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())


@pytest.mark.parametrize('change', ['radical', 'charge', 'query', 'declared_H', 'double_bond', 'no_warning'])
def test_hydrogen_rule_keeps_its_guards_on_a_real_source_atom(change):
    # Controlled mutations of the actual reviewed source, not invented molecules.
    doc = read_document(fixture('US20080125424A1_C00020.cdx'))
    node = doc.cdxml.getroot().xpath('.//n[@id="54"]')[0]
    if change == 'radical': node.set('Radical', 'Doublet')
    elif change == 'charge': node.set('Charge', '1')
    elif change == 'query': node.set('ElementList', '7 8')
    elif change == 'declared_H': node.set('NumHydrogens', '1')
    elif change == 'double_bond':
        bond = next(b for b in node.getparent().findall('b') if '54' in (b.get('B'), b.get('E')))
        bond.set('Order', '2')
    elif change == 'no_warning':
        for item in [node, *node.findall('t')]: item.attrib.pop('Warning', None)
    before = dict(node.attrib)
    assert not repair_invalid_nitrogen_hydrogens(node)
    assert dict(node.attrib) == before


@pytest.mark.parametrize('filename,fragment,port', [
    ('US20210052514A1_C00011.cdx', '72', '73'),
    ('US20210052514A1_C00012.cdx', '40', '41'),
    ('US20210052514A1_C00015.cdx', '159', '161'),
])
def test_methoxy_keeps_its_attachment_port_and_source_identity(filename, fragment, port):
    mol = objects(filename)[fragment]
    assert Chem.MolToSmiles(mol) == '*OC'
    atom = source_atom(mol, port)
    assert atom.GetAtomicNum() == 0 and atom.GetDegree() == 1
    assert atom.GetAtomMapNum() == 0


def test_methoxy_is_independent_of_xml_node_order(monkeypatch):
    path = fixture('US20210052514A1_C00011.cdx')
    expected = {k: Chem.MolToSmiles(m) for k, m in objects(path.name).items()}
    doc = read_document(path)
    for fragment in doc.cdxml.getroot().iter('fragment'):
        children = list(fragment)
        slots = [i for i, child in enumerate(children) if child.tag == 'n']
        for i, node in zip(slots, reversed(fragment.findall('n'))): children[i] = node
        fragment[:] = children
    module = importlib.import_module('cdx_reader.cdx_reader')
    monkeypatch.setattr(module, 'read_document', lambda _: copy.deepcopy(doc))
    assert {k: Chem.MolToSmiles(m) for k, m in objects(path.name).items()} == expected


def test_expansion_preserves_the_potassium_component_of_the_reagent():
    mol = objects('US20110306591A1_C00049.cdx')['36']
    assert Chem.MolToSmiles(mol) == 'CCOC(=O)[S-].[K+]'
    assert Chem.GetFormalCharge(mol) == 0
    assert len(Chem.GetMolFrags(mol)) == 2


def test_expansion_preserves_hydrochloride_and_existing_stereo():
    mol = objects('US20140057938A1_C00038.cdx')['1258']
    assert Chem.MolToInchiKey(mol) == 'KJPNGCJAYJKJOJ-FVGYRXGTSA-N'
    assert any(Chem.MolToSmiles(f) == 'Cl' for f in Chem.GetMolFrags(mol, asMols=True))
    assert len(Chem.FindMolChiralCenters(mol)) == 1


def test_disconnected_nested_group_does_not_hide_the_amine_hydrochloride():
    mol = objects('US20130116228A1_C00023.cdx')['40']
    assert Chem.MolToSmiles(mol) == 'Cl.FC1(C2CC2)CNC1'
    assert Chem.MolToInchiKey(mol) == 'NVPLJZDJMOUWKD-UHFFFAOYSA-N'
    assert mol.HasProp('_CDX_NESTED_EXPANSIONS')


def test_nested_shared_anchor_retains_the_bromide_and_source_atoms():
    mol = objects('US20130178505A1_C00217.cdx')['189']
    assert Chem.MolToInchiKey(mol) == 'RNJUGFMPFUUEGM-NDUABGMUSA-N'
    bromide = source_atom(mol, '196')
    assert (bromide.GetAtomicNum(), bromide.GetFormalCharge(), bromide.GetDegree()) == (35, -1, 0)
    assert mol.HasProp('_CDX_NESTED_EXPANSIONS')


def test_source_co2_is_not_silently_discarded_as_an_acetoxy_repair():
    mol = objects('US20110311444A1_C00025.cdx')['14']
    parts = Chem.GetMolFrags(mol, asMols=True)
    assert len(parts) == 2
    assert 'O=C=O' in [Chem.MolToSmiles(f) for f in parts]


@pytest.mark.parametrize('filename,property_name', [
    ('US20070123494A1_C00048.cdx', 'BracketType'),
    ('US20080161371A1_C00151.cdx', 'SymbolType'),
    ('US20100197697A1_C00155.cdx', 'RectangleType'),
])
def test_graphic_integer_widths_are_lossless_and_idempotent(filename, property_name):
    raw = fixture(filename).read_bytes()
    normalized, repairs = normalize_integer_widths(raw)
    relevant = [r for r in repairs if r['property'] == property_name]
    assert relevant
    for r in relevant:
        start = r['offset']
        assert int.from_bytes(raw[start+2:start+4], 'little') == r['old_width']
        assert int.from_bytes(raw[start+4:start+4+r['old_width']], 'little') == r['value']
        assert r['width'] == 2
    assert normalize_integer_widths(normalized) == (normalized, [])
    doc = read_document(raw)
    assert doc._cdx_binary_repairs == repairs
    assert objects(filename)


@pytest.mark.parametrize('filename,fragment,key', [
    ('US20130178505A1_C00240.cdx', '14', 'LSXVNKZFJLYKDP-AWEZNQCLSA-N'),
    ('US20130317001A1_C00741.cdx', '113', 'OVHJDMLYVHRGRM-TXEJJXNPSA-N'),
])
def test_hollow_wedge_preserves_the_drawn_stereo_and_records_its_style(filename, fragment, key):
    mol = objects(filename)[fragment]
    assert Chem.MolToInchiKey(mol) == key
    assert Chem.MolToInchiKey(Chem.MolFromSmiles(Chem.MolToSmiles(mol))) == key
    repairs = json.loads(mol.GetProp('_CDX_STEREO_REPAIRS'))
    assert len(repairs) == 1
    assert repairs[0]['encoded'] == 'HollowWedgeBegin'
    assert repairs[0]['interpreted'] == 'WedgeBegin'
    assert repairs[0]['fragment_id'] == fragment


def test_reversing_the_hollow_bond_end_keeps_the_same_configuration(monkeypatch):
    filename = 'US20130178505A1_C00240.cdx'
    doc = read_document(fixture(filename))
    bond = next(b for b in doc.cdxml.getroot().iter('b') if b.get('Display') == 'HollowWedgeBegin')
    begin, end = bond.get('B'), bond.get('E')
    bond.set('B', end); bond.set('E', begin); bond.set('Display', 'HollowWedgeEnd')
    module = importlib.import_module('cdx_reader.cdx_reader')
    monkeypatch.setattr(module, 'read_document', lambda _: copy.deepcopy(doc))
    mol = objects(filename)['14']
    assert Chem.MolToInchiKey(mol) == 'LSXVNKZFJLYKDP-AWEZNQCLSA-N'
    repair = json.loads(mol.GetProp('_CDX_STEREO_REPAIRS'))[0]
    assert (repair['begin'], repair['end'], repair['encoded']) == (end, begin, 'HollowWedgeEnd')


def test_hollow_normalization_leaves_other_bond_styles_unchanged():
    doc = read_document(fixture('US20130317001A1_C00741.cdx'))
    others = [b for b in doc.cdxml.getroot().iter('b') if b.get('Display') != 'HollowWedgeBegin']
    before = [dict(b.attrib) for b in others]
    assert len(normalize_hollow_wedges(doc.cdxml)) == 1
    assert [dict(b.attrib) for b in others] == before
    assert normalize_hollow_wedges(doc.cdxml) == []


@pytest.mark.parametrize('case,count', [('US20210079003A1_C00690', 5), ('US20210017189A1_C00004', 14)])
def test_combined_api_retains_every_object_of_the_reaction_scheme(case, count):
    path = DATA.parent / 'patent_regressions' / (case + '.cdx')
    mols = read_cdx_molecule(path, as_single=False)
    assert len(mols) == count and all(m is not None for m in mols)
    combined = read_cdx_molecule(path)
    expected = sorted(Chem.MolToSmiles(f) for m in mols for f in Chem.GetMolFrags(m, asMols=True))
    assert sorted(Chem.MolToSmiles(f) for f in Chem.GetMolFrags(combined, asMols=True)) == expected
