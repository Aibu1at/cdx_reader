"""Source-reviewed real CDX: functional labels, nested graphs and binary widths."""
import hashlib
import json
from pathlib import Path

import pytest
from rdkit import Chem

from cdx_reader import read_cdx_molecule
from cdx_reader._io import read_document, normalize_integer_widths

DATA=Path(__file__).with_name('data')/'recovery_20261007'
MANIFEST={r['file']:r for r in json.loads((DATA/'manifest.json').read_text())}


def fixture(name):
    path=DATA/name
    assert hashlib.sha256(path.read_bytes()).hexdigest()==MANIFEST[name]['sha256']
    return path


def source_atom(mol,identifier):
    return next(a for a in mol.GetAtoms()if a.HasProp('_CDX_ATOM_ID')
                and a.GetProp('_CDX_ATOM_ID')==identifier)


@pytest.mark.parametrize('mode',['path','bytes','objects'])
def test_nested_sulfone_keeps_the_original_r_center(mode):
    p=fixture('US20140296222A1_C00094.cdx')
    m=read_cdx_molecule(p.read_bytes()if mode=='bytes'else p,as_single=mode!='objects')
    if mode=='objects':
        assert len(m)==1
        assert m[0].GetProp('CDXML_FRAG_ID')=='14'
        m=m[0]
    # Original TIF: N above, CH2-S to the right, pyridyl to the left, methyl back.
    assert Chem.MolToInchiKey(m)=='PRCCLLOFFCXGKV-JTQLQIEISA-N'
    assert Chem.FindMolChiralCenters(m)==[(6,'R')]
    assert source_atom(m,'86').GetAtomicNum()==16
    assert source_atom(m,'86').GetProp('_CDX_GROUP_SOURCE_ATOM')=='43'
    assert len([a for a in m.GetAtoms()if a.GetAtomicNum()==0])==0


def test_recovered_173_retains_the_amide_n_missing_from_its_mol_export():
    m=read_cdx_molecule(fixture('US20140296222A1_C00676.cdx'))
    # Read directly from original TIF, not the MOL export's lost amide N.
    assert m.HasSubstructMatch(Chem.MolFromSmarts('cNC(=O)c1ccc(F)cn1'))
    assert Chem.MolToInchiKey(m)=='VXLAYHUGSSZHSG-KRWDZBQOSA-N'
    assert [cip for _,cip in Chem.FindMolChiralCenters(m)]==['R']


@pytest.mark.parametrize('name,pattern',[
 ('US20120238577A1_C00289.cdx','[CH](c1ccccc1)(NC(=O))C(=O)[OH]'),
 ('US20120238577A1_C00283.cdx','[CH](c1ccccc1)(NC(=O))C(=O)O[CH3]'),
 ('US20110311444A1_C00011.cdx','cO[CH3]'),
 ('US20110311444A1_C00051.cdx','c[NH][CH3]'),
 ('US20180022705A1_C00106.cdx','N1CCN(C(=O)OC(C)(C)C)CC1'),
 ('US20120040968A1_C00212.cdx','N([CH3])[CH3]'),
 ('US20120289572A1_C00043.cdx','[CH](F)F'),
 ('US20120041190A1_C00690.cdx','S(=O)(=O)'),
 ('US20120289572A1_C00024.cdx','C(=O)[NH2]'),
])
def test_condensed_family_has_the_drawn_connectivity(name,pattern):
    ms=read_cdx_molecule(fixture(name),as_single=False)
    assert any(m is not None and m.HasSubstructMatch(Chem.MolFromSmarts(pattern))for m in ms)


@pytest.mark.parametrize('name,anchor,element,labels',[
 ('US20130237515A1_C00029.cdx','31',7,{'R1','R2'}),
 ('US20180022705A1_C00010.cdx','5',7,{'R7','R8'}),
 ('US20180022705A1_C00017.cdx','789',8,{'R15'}),
 ('US20180022705A1_C00015.cdx','747',7,{'R14'}),
 ('US20180022705A1_C00002.cdx','117',6,{'R12'}),
])
def test_composite_groups_preserve_variables_and_honest_source_addresses(name,anchor,element,labels):
    m=read_cdx_molecule(fixture(name))
    assert source_atom(m,anchor).GetAtomicNum()==element
    ports=[a for a in m.GetAtoms()if a.HasProp('_CDX_GROUP_SOURCE_ATOM')
           and a.GetProp('_CDX_GROUP_SOURCE_ATOM')==anchor and a.GetAtomicNum()==0]
    assert {a.GetProp('dummyLabel')for a in ports}==labels
    assert all(not a.HasProp('_CDX_ATOM_ID')for a in ports)
    assert all(a.GetAtomMapNum()==0 for a in ports)
    assert all(not a.HasProp('_CDX_UNPARSED_LABEL')for a in ports)


@pytest.mark.parametrize('name,identifier,raw,normalized,element',[
 ('US20120238577A1_C00269.cdx','85','Oh','OH',8),
 ('US20140128344A1_C00001.cdx','34','OH.','OH',8),
 ('US20140296222A1_C00054.cdx','30','R1,','R1',0),
 ('US20210052514A1_C00024.cdx','240','R¢','R′',0),
])
def test_semantic_normalization_retains_the_original_caption(name,identifier,raw,normalized,element):
    ms=read_cdx_molecule(fixture(name),as_single=False)
    atom=next(a for m in ms if m is not None for a in m.GetAtoms()
              if a.HasProp('_CDX_ATOM_ID')and a.GetProp('_CDX_ATOM_ID')==identifier)
    assert atom.GetAtomicNum()==element
    assert atom.GetProp('_CDX_ATOM_LABEL')==raw
    assert atom.GetProp('_CDX_NORMALIZED_LABEL')==normalized


def test_nested_isocyanide_under_an_nc_label_changes_only_by_the_recorded_nitrile_rule():
    # The text grammar does not replace a nested graph. Since 2026-10-28 the NC label (C next to
    # the bond, nested N(+)#C(-) bonded through N) is read as the drawn nitrile, and says so; the
    # reference of this compound (61) is the nitrile too.
    m=read_cdx_molecule(fixture('US20180022705A1_C00108.cdx'))
    assert not m.HasSubstructMatch(Chem.MolFromSmarts('[C-]#[N+]'))
    assert m.HasSubstructMatch(Chem.MolFromSmarts('[CH2]C#N'))
    assert [r['label'] for r in json.loads(m.GetProp('_CDX_LABEL_REPAIRS'))]==['NC']


def test_valid_nested_variable_keeps_real_source_ids_and_its_stereochemical_neighbor():
    ms=read_cdx_molecule(fixture('US20110135650A1_C00151.cdx'),as_single=False)
    m=ms[0]
    oxygen=source_atom(m,'154')
    variable=source_atom(m,'155')
    center=source_atom(m,'17')
    assert oxygen.GetAtomicNum()==8
    assert variable.GetAtomicNum()==0 and variable.GetProp('dummyLabel')=='R2'
    assert m.GetBondBetweenAtoms(oxygen.GetIdx(),variable.GetIdx()).GetBondType()==Chem.BondType.SINGLE
    assert m.GetBondBetweenAtoms(oxygen.GetIdx(),center.GetIdx()).GetBondType()==Chem.BondType.SINGLE
    assert center.GetChiralTag()==Chem.ChiralType.CHI_TETRAHEDRAL_CCW
    assert not oxygen.HasProp('_CDX_GROUP_SOURCE_ATOM')


@pytest.mark.parametrize('name,count',[
 ('US20100256227A1_C00007.cdx',3),('US20080119443A1_C00004.cdx',7),
])
def test_graphic_enum_widths_recover_all_original_objects_without_renumbering(name,count):
    blob=fixture(name).read_bytes()
    normalized,repairs=normalize_integer_widths(blob)
    assert repairs and all(r['property']in ('GraphicType','LineType','ArrowType')for r in repairs)
    assert normalize_integer_widths(normalized)==(normalized,[])
    doc=read_document(blob)
    ids=[f.get('id')for f in doc.cdxml.getroot().iter('fragment')
         if not any(p.tag=='fragment'for p in f.iterancestors())]
    mols=read_cdx_molecule(blob,as_single=False)
    assert len(ids)==len(mols)==count
    assert [m.GetProp('CDXML_FRAG_ID')for m in mols]==ids


@pytest.mark.parametrize('name', ['US20140171435A1_C00046.cdx','US20210052514A1_C00456.cdx'])
def test_invalid_source_elements_are_not_reinterpreted_as_labels(name):
    # The original TIF/CDX really encodes bonded Ar and triple-bonded F.
    assert read_cdx_molecule(fixture(name)) is None
