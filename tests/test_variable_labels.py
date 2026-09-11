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
@pytest.mark.parametrize('suffix', ['', '1', '12', 'a', 'g', '1a', '1b', '27z'])
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
    assert reader.get_fragment_smiles('R' + suffix) == (None, None, True)
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


@pytest.mark.parametrize('label', ['R2', 'R1a', 'R1b'])
def test_r_label_survives_adjacent_abbreviation_expansion(monkeypatch, label):
    xml = tree(f'''
      <n id="10" p="0 0" NodeType="Unspecified"><t><s>{label}</s></t></n>
      <n id="20" p="1 0"/>
      <n id="30" p="2 0" NodeType="Unspecified"><t><s>MeO</s></t></n>
      <b id="40" B="10" E="20"/><b id="50" B="20" E="30"/>
    ''')
    mol = read_tree(monkeypatch, xml)
    assert Chem.MolToSmiles(mol) == '*COC'
    atom = next(a for a in mol.GetAtoms() if a.GetAtomicNum() == 0)
    assert atom.GetProp('_CDX_ATOM_ID') == '10'
    assert atom.GetProp('dummyLabel') == label
    if label == 'R2':
        assert atom.GetIntProp('_MolFileRLabel') == 2
    else:
        assert not atom.HasProp('_MolFileRLabel')
    assert atom.GetAtomMapNum() == 0


@pytest.mark.parametrize('as_single', [True, False])
def test_indexed_letter_labels_remain_distinct_and_may_repeat(monkeypatch, as_single):
    xml = tree('''
      <n id="10" p="0 0" NodeType="Unspecified"><t><s>R</s><s>1</s><s>a</s></t></n>
      <n id="20" p="1 0"/>
      <n id="30" p="2 0" NodeType="Unspecified"><t><s>R1b</s></t></n>
      <n id="40" p="1 1" NodeType="Unspecified"><t><s>R1a</s></t></n>
      <b id="50" B="10" E="20"/><b id="60" B="20" E="30"/>
      <b id="70" B="20" E="40"/>
    ''')
    result = read_tree(monkeypatch, xml, as_single=as_single)
    mol = result if as_single else result[0]
    assert Chem.MolToSmiles(mol) == '*C(*)*'
    atoms = {a.GetProp('_CDX_ATOM_ID'): a for a in mol.GetAtoms() if a.GetAtomicNum() == 0}
    assert {key: atom.GetProp('dummyLabel') for key, atom in atoms.items()} == {
        '10': 'R1a', '30': 'R1b', '40': 'R1a'}
    assert all(not atom.HasProp('_MolFileRLabel') and atom.GetAtomMapNum() == 0
               for atom in atoms.values())


def test_fix_false_does_not_mutate_tree():
    xml = tree('<n id="10" NodeType="Unspecified"><t><s>R1</s></t></n>')
    before = etree.tostring(xml)
    parse_cdxml(xml, fix=False)
    assert etree.tostring(xml) == before


@pytest.mark.parametrize('label', ['R8²', 'R²'])
def test_literal_superscript_digits_keep_legacy_placeholder(monkeypatch, label):
    xml = tree(f'''
      <n id="10" p="0 0" NodeType="Unspecified"><t><s>{label}</s></t></n>
      <n id="20" p="1 0"/><b id="30" B="10" E="20"/>
    ''')
    mol = read_tree(monkeypatch, xml)
    atom = next(a for a in mol.GetAtoms() if a.GetProp('_CDX_ATOM_ID') == '10')
    assert atom.GetAtomicNum() == 0
    assert atom.GetProp('_CDX_ATOM_LABEL') == label
    assert atom.GetProp('dummyLabel') == label
    assert not atom.HasProp('_MolFileRLabel')
    assert reader.get_fragment_smiles(label) == (None, None, True)


@pytest.mark.parametrize('label', ['R1', 'R1a'])
def test_variable_does_not_block_existing_valence_repair(monkeypatch, label):
    xml = tree(f'''
      <n id="10" p="0 0"/><n id="20" p="1 0" Element="7" NumHydrogens="2"/>
      <n id="30" p="2 0"/>
      <b id="60" B="10" E="20"/><b id="70" B="20" E="30"/>
      <n id="40" p="0 2" NodeType="Unspecified"><t><s>{label}</s></t></n>
      <n id="50" p="1 2"/><b id="80" B="40" E="50"/>
    ''')
    mol = read_tree(monkeypatch, xml)
    assert mol is not None
    assert Chem.MolToSmiles(mol) == '*C.CNC'
    atom = next(a for a in mol.GetAtoms() if a.GetProp('_CDX_ATOM_ID') == '40')
    assert atom.GetProp('dummyLabel') == label
    assert atom.GetNumExplicitHs() == 0


def test_variable_does_not_enable_unreviewed_hydrogen_repairs(monkeypatch):
    xml = tree('''
      <n id="10" p="0 0"/><n id="20" p="1 0" Element="8" NumHydrogens="2"/>
      <n id="30" p="2 0"/>
      <b id="60" B="10" E="20"/><b id="70" B="20" E="30"/>
      <n id="40" p="0 2" NodeType="Unspecified"><t><s>R1a</s></t></n>
      <n id="50" p="1 2"/><b id="80" B="40" E="50"/>
    ''')
    assert read_tree(monkeypatch, xml) is None


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
