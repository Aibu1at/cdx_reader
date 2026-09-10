"""Variable metadata and abbreviation attachment regressions."""
from types import SimpleNamespace


from lxml import etree


import pytest


from rdkit import Chem


from cdx_reader import parse_cdxml, read_cdx_molecule


from cdx_reader import cdx_reader as reader




def tree(body):
    return etree.ElementTree(etree.fromstring(
        ('<CDXML BondLength="1.5"><page id="1"><fragment id="2">' + body +
         '</fragment></page></CDXML>').encode()))


def read_tree(monkeypatch, xml, as_single=True):
    doc = SimpleNamespace(cdxml=xml, to_cdxml=lambda: etree.tostring(xml).decode())
    monkeypatch.setattr(reader.cdxml_converter, 'read_cdx', lambda source: doc)
    return read_cdx_molecule(b'synthetic CDXML adapter', as_single=as_single)


@pytest.mark.parametrize('node_type', ['Unspecified', 'GenericNickname'])
@pytest.mark.parametrize('suffix', ['', '1', '12', 'a', 'g'])
def test_r_labels_from_styled_segments(monkeypatch, node_type, suffix):
    xml = tree(f'''
      <n id="10" p="0 0" NodeType="{node_type}" GenericNickname="R">
        <t><s face="96">R</s><s face="64">{suffix}</s></t></n>
      <n id="20" p="1 0"/><b id="30" B="10" E="20" Order="1"/>
    ''')
    mol = read_tree(monkeypatch, xml)
    atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
    assert len(atoms) == 1
    atom = atoms[0]
    assert atom.GetProp('_CDX_ATOM_ID') == '10'
    assert atom.GetProp('_CDX_ATOM_LABEL') == 'R' + suffix
    assert atom.GetProp('dummyLabel') == 'R' + suffix
    assert atom.GetAtomMapNum() == 0
    if suffix.isdigit():
        assert atom.GetIntProp('_MolFileRLabel') == int(suffix)
    else:
        assert not atom.HasProp('_MolFileRLabel')


@pytest.mark.parametrize('as_single', [True, False])
def test_repeated_r_numbers_are_not_unique_attachment_ids(monkeypatch, as_single):
    xml = tree('''
      <n id="10" p="0 0" NodeType="Unspecified"><t><s>R1</s></t></n>
      <n id="20" p="1 0"/>
      <n id="30" p="2 0" NodeType="Unspecified"><t><s>R1</s></t></n>
      <b id="40" B="10" E="20"/><b id="50" B="20" E="30"/>
    ''')
    result = read_tree(monkeypatch, xml, as_single=as_single)
    mol = result if as_single else result[0]
    assert Chem.MolToSmiles(mol) == '*C*'
    labels = {a.GetProp('_CDX_ATOM_ID'): a.GetIntProp('_MolFileRLabel')
              for a in mol.GetAtoms() if a.HasProp('_MolFileRLabel')}
    assert labels == {'10': 1, '30': 1}


def test_r_label_survives_adjacent_abbreviation_expansion(monkeypatch):
    xml = tree('''
      <n id="10" p="0 0" NodeType="Unspecified"><t><s>R2</s></t></n>
      <n id="20" p="1 0"/>
      <n id="30" p="2 0" NodeType="Unspecified"><t><s>MeO</s></t></n>
      <b id="40" B="10" E="20"/><b id="50" B="20" E="30"/>
    ''')
    mol = read_tree(monkeypatch, xml)
    assert Chem.MolToSmiles(mol) == '*COC'
    atom = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 0)
    assert atom.GetProp('_CDX_ATOM_ID') == '10'
    assert atom.GetIntProp('_MolFileRLabel') == 2
    assert atom.GetAtomMapNum() == 0


def test_fix_false_does_not_mutate_tree():
    xml = tree('<n id="10" NodeType="Unspecified"><t><s>R1</s></t></n>')
    before = etree.tostring(xml)
    parse_cdxml(xml, fix=False)
    assert etree.tostring(xml) == before


@pytest.mark.parametrize('node_type,label', [
    ('Unspecified', 'X1'), ('Unspecified', 'Y2'), ('Unspecified', 'L'),
    ('Unspecified', 'L12'), ('Unspecified', 'A'), ('Unspecified', 'B1'),
    ('GenericNickname', 'B'), ('GenericNickname', 'Het'), ('GenericNickname', 'Ar'),
])
def test_alternate_variable_labels_preserve_text_and_source_node(monkeypatch, node_type, label):
    xml = tree(f'<n id="10" p="0 0" NodeType="{node_type}"><t><s>{label[0]}</s>'
               f'<s>{label[1:]}</s></t></n><n id="20" p="1 0"/><b id="30" B="10" E="20"/>')
    mol = read_tree(monkeypatch, xml)
    atom = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 0)
    assert atom.GetProp('_CDX_ATOM_ID') == '10'
    assert atom.GetProp('_CDX_ATOM_LABEL') == label
    assert atom.GetProp('_CDX_NODE_TYPE') == node_type
    assert atom.GetProp('dummyLabel') == label
    assert not atom.HasProp('_MolFileRLabel')
    assert atom.GetAtomMapNum() == 0


def test_cut_at_atom_zero_does_not_consume_source_r_isotope(monkeypatch):
    xml = tree('''
      <n id="10" p="0 0" NodeType="Unspecified"><t><s>MeO</s></t></n>
      <n id="20" p="1 0"/><n id="30" p="2 0"/>
      <n id="40" p="3 0" NodeType="Unspecified" Isotope="2"><t><s>R2</s></t></n>
      <b id="50" B="10" E="20"/><b id="60" B="20" E="30"/>
      <b id="70" B="30" E="40"/>
    ''')
    mol = read_tree(monkeypatch, xml)
    assert Chem.MolToSmiles(mol) == '[2*]CCOC'
    r = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 0)
    assert r.GetIsotope() == 2
    assert r.GetProp('_CDX_ATOM_ID') == '40'
    assert r.GetIntProp('_MolFileRLabel') == 2
    assert r.GetAtomMapNum() == 0
