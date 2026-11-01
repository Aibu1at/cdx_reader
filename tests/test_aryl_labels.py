"""Source-confirmed composite labels; no automatic aryl selection."""
from pathlib import Path
import pytest
from lxml import etree
from rdkit import Chem
from cdx_reader import read_cdx_molecule
from cdx_reader._normalization import variable_label
from tests.test_variable_labels import tree, read_tree


def test_original_ara_arb_table_core():
    raw = (Path(__file__).parent/'data/US20210122713A1-C00143.CDX').read_bytes()
    mol = read_cdx_molecule(raw)
    atoms = {a.GetProp('_CDX_ATOM_ID'):a for a in mol.GetAtoms() if a.HasProp('_CDX_ATOM_ID')}
    assert set(atoms) == {'11','13','15','17','19','21','29','31','33','37','39','41','43','45'}
    for id, label, neighbor in [('43','ArA','33'), ('45','ArB','39'), ('29','R3','15')]:
        a=atoms[id]
        assert a.GetAtomicNum()==0
        assert a.GetProp('_CDX_ATOM_LABEL')==a.GetProp('dummyLabel')==label
        assert a.GetNeighbors()[0].GetProp('_CDX_ATOM_ID')==neighbor
        assert a.GetAtomMapNum()==0


@pytest.mark.parametrize('label', ['ArA','ArB'])
def test_repeated_composite_label_and_neighbor_expansion(monkeypatch,label):
    xml=tree(f'''<n id="10" p="0 0" NodeType="Unspecified"><t><s>Ar</s><s>{label[-1]}</s></t></n>
      <n id="20" p="1 0"/><n id="30" p="2 0" NodeType="Unspecified"><t><s>{label}</s></t></n>
      <n id="40" p="1 1" NodeType="Unspecified"><t><s>MeO</s></t></n>
      <b id="50" B="10" E="20"/><b id="60" B="20" E="30"/><b id="70" B="20" E="40"/>''')
    m=read_tree(monkeypatch,xml)
    assert Chem.MolToSmiles(m)=='*C(*)OC'
    assert sorted(a.GetProp('dummyLabel') for a in m.GetAtoms() if a.GetAtomicNum()==0)==[label,label]


@pytest.mark.parametrize('attributes,label', [
    ('Element="18"','Ar'), ('Element="18" NodeType="Unspecified"','ArA'),
    ('Element="6" NodeType="Unspecified"','ArB'), ('NodeType="Element"','ArA'),
    ('NodeType="Unspecified" Charge="1"','ArA'), ('NodeType="Unspecified" Isotope="13"','ArA'),
    ('NodeType="Unspecified" Radical="Doublet"','ArB'), ('NodeType="Unspecified" ElementList="6 7"','ArA'),
    ('NodeType="Unspecified"','ArC'), ('NodeType="Unspecified"','MeO'),
    ('NodeType="Unspecified"','Na'), ('NodeType="Unspecified"','Boc'),
])
def test_composite_rule_does_not_reclassify_chemistry(attributes,label):
    n=etree.fromstring(f'<n id="10" {attributes}><t><s>{label}</s></t></n>')
    assert variable_label(n, {'MeO','Boc'}) is None


def test_complete_nested_graph_precedes_composite_text():
    xml=tree('''<n id="10" NodeType="Unspecified"><t><s>ArA</s></t>
      <fragment id="11"><n id="12" Element="6"/><n id="13" NodeType="ExternalConnectionPoint"/>
      <b id="14" B="12" E="13"/></fragment></n><n id="20"/><b id="30" B="10" E="20"/>''')
    assert variable_label(xml.find('.//n'),set()) is None
