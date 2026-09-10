"""Conservative chemistry repairs and source-graph boundaries."""
from types import SimpleNamespace


from lxml import etree


import pytest


from rdkit import Chem


from cdx_reader import parse_cdxml, read_cdx_molecule


from cdx_reader import cdx_reader as reader


from cdx_reader._normalization import repair_invalid_nitrogen_hydrogens


def tree(body):
    return etree.ElementTree(etree.fromstring(
        ('<CDXML BondLength="1.5"><page id="1"><fragment id="2">' + body +
         '</fragment></page></CDXML>').encode()))


def read_tree(monkeypatch, xml, as_single=True):
    doc = SimpleNamespace(cdxml=xml, to_cdxml=lambda: etree.tostring(xml).decode())
    monkeypatch.setattr(reader.cdxml_converter, 'read_cdx', lambda source: doc)
    return read_cdx_molecule(b'synthetic CDXML adapter', as_single=as_single)


@pytest.mark.parametrize('element,label,charge,isotope', [
    ('11', 'Na+', '1', '23'), ('14', 'Si', '0', '29'),
    ('16', 'S', '0', '34'), ('7', 'NH2', '0', '15'),
])
def test_explicit_elements_charge_and_isotope_are_preserved(element, label, charge, isotope):
    xml = tree(f'<n id="10" Element="{element}" Charge="{charge}" Isotope="{isotope}">'
               f'<t><s>{label}</s></t></n>')
    parse_cdxml(xml)
    atom = xml.find('.//n')
    assert atom.get('Element') == element
    assert atom.get('Charge') == charge
    assert atom.get('Isotope') == isotope


def test_hn_repair_without_explicit_element():
    xml = tree('<n id="10" NodeType="Unspecified"><t><s>H&#13;N</s></t></n>')
    parse_cdxml(xml)
    assert xml.find('.//n').get('Element') == '7'


@pytest.mark.parametrize('extra,order,warning,expected', [
    ('', '1', 'An atom in this label has an invalid valence.', True),
    ('Radical="Doublet"', '1', 'invalid valence', False),
    ('Charge="1"', '1', 'invalid valence', False),
    ('NodeType="GenericNickname"', '1', 'invalid valence', False),
    ('ElementList="6 7"', '1', 'invalid valence', False),
    ('', '2', 'invalid valence', False),
    ('', '1', '', False),
])
def test_hydrogen_repair_is_narrow(extra, order, warning, expected):
    xml = tree(f'''
      <n id="10" Element="7" NumHydrogens="0" Warning="{warning}" {extra}/>
      <n id="20"/><n id="30"/>
      <b id="40" B="10" E="20" Order="{order}"/>
      <b id="50" B="10" E="30" Order="1"/>
    ''')
    atom = xml.find('.//n')
    assert repair_invalid_nitrogen_hydrogens(atom) is expected
    assert (atom.get('NumHydrogens') is None) is expected


@pytest.mark.parametrize('external_points', [1, 2])
def test_nested_graph_requires_unambiguous_attachment(external_points):
    extra = ('<n id="15" NodeType="ExternalConnectionPoint"/>'
             '<b id="16" B="15" E="12"/>' if external_points == 2 else '')
    xml = tree(f'''
      <n id="10" NodeType="Unspecified"><t><s>O2N</s></t>
        <fragment id="11"><n id="12" Element="7"/>
          <n id="13" NodeType="ExternalConnectionPoint"/>
          <b id="14" B="12" E="13"/>{extra}</fragment></n>
      <n id="20"/><b id="30" B="10" E="20"/>
    ''')
    parse_cdxml(xml)
    atom = xml.find('.//n')
    assert atom.get('NodeType') == ('Fragment' if external_points == 1 else 'Unspecified')


@pytest.mark.parametrize('label,element', [('U', '92'), ('V', '23'), ('K', '19')])
def test_ambiguous_neutral_variable_letters_keep_legacy_behavior(label, element):
    xml = tree(f'<n id="10" Element="{element}"><t><s>{label}</s></t></n>')
    parse_cdxml(xml)
    assert xml.find('.//n').get('Element') == '6'


def test_explicit_potassium_ion_is_preserved():
    xml = tree('<n id="10" Element="19" Charge="1"><t><s>K</s></t></n>')
    parse_cdxml(xml)
    assert xml.find('.//n').get('Element') == '19'
    assert xml.find('.//n').get('Charge') == '1'


@pytest.mark.parametrize('suffix,isotope', [('D', '2'), ('T', '3')])
@pytest.mark.parametrize('bonds,order,warning,expected', [
    (1, '1', '', '1'), (2, '1', '', '6'), (1, '2', '', '6'),
    (0, '1', 'invalid valence', '6'),
])
def test_isotopic_hydrogen_requires_a_valid_attachment(suffix, isotope, bonds, order, warning, expected):
    attached = ''.join(f'<n id="{20+i}"/><b id="{30+i}" B="10" E="{20+i}" Order="{order}"/>'
                       for i in range(bonds))
    xml = tree(f'<n id="10" Element="1" Isotope="{isotope}" Warning="{warning}">'
               f'<t><s>{suffix}</s></t></n>{attached}')
    parse_cdxml(xml)
    assert xml.find('.//n').get('Element') == expected
    assert xml.find('.//n').get('Isotope') == isotope


@pytest.mark.parametrize('label,element', [('B', '5'), ('Br', '35'), ('Al', '13'), ('Mg', '12')])
def test_alternate_variables_do_not_reclassify_explicit_elements(label, element):
    xml = tree(f'<n id="10" NodeType="Unspecified" Element="{element}"><t><s>{label}</s></t></n>')
    parse_cdxml(xml)
    assert xml.find('.//n').get('Element') == element


def test_external_connection_is_consumed_but_nested_atoms_and_bonds_survive(monkeypatch):
    xml = tree('''
      <n id="10" p="0 0" NodeType="Fragment"><t><s>Et</s></t>
        <fragment id="11"><n id="12" p="0 0"/><n id="13" p="1 0"/>
          <n id="14" p="2 0" NodeType="ExternalConnectionPoint"/>
          <b id="15" B="12" E="13"/><b id="16" B="13" E="14"/>
        </fragment></n><n id="20" p="3 0" Element="8"/>
      <b id="30" B="10" E="20"/>
    ''')
    mol = read_tree(monkeypatch, xml)
    assert Chem.MolToSmiles(mol) == 'CCO'
    atoms = {a.GetProp('_CDX_ATOM_ID'): a.GetIdx() for a in mol.GetAtoms()}
    assert set(atoms) == {'12', '13', '20'}
    assert mol.GetBondBetweenAtoms(atoms['13'], atoms['20']) is not None
