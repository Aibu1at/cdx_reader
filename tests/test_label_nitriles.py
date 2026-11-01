"""A cyano label NC that ChemDraw expanded as an isocyanide is read as the nitrile it shows.

US20210017156A1 C00098 is Example 38, named ...pyrimidine-5-carbonitrile; its CDX (and MOL) encode
N(+)#C(-) bonded through N under the label NC, with the C next to the bond. US20200237771A1 C00162
draws reagent 2, ethyl-2-isocyanoacetate (p-0477), with its own N+ and C- atoms: a real isocyanide.
"""
import hashlib
import json
from pathlib import Path

from rdkit import Chem

from cdx_reader import read_cdx_molecule

DATA = Path(__file__).with_name('data') / 'label_nitriles'
MANIFEST = json.loads((DATA / 'manifest.json').read_text())
ISOCYANIDE = Chem.MolFromSmarts('[C-]#[N+]')


def fixture(name):
    path = DATA / name
    assert hashlib.sha256(path.read_bytes()).hexdigest() == MANIFEST[name]['sha256']
    return path


def test_an_nc_label_expanded_as_an_isocyanide_is_read_as_the_drawn_nitrile():
    (mol,) = read_cdx_molecule(fixture('US20210017156A1_C00098.cdx'), as_single=False)
    assert Chem.MolToSmiles(mol) == 'Cn1cnc(S(=O)(=O)N2CCC(Nc3ncc(C#N)c(-c4cnn(CC(F)(F)F)c4)n3)CC2)c1'
    assert not mol.HasSubstructMatch(ISOCYANIDE)
    (repair,) = json.loads(mol.GetProp('_CDX_LABEL_REPAIRS'))
    assert repair['label'] == 'NC' and repair['fragment_id'] == mol.GetProp('CDXML_FRAG_ID')
    assert repair['encoded'].startswith('isocyanide') and repair['read'].startswith('nitrile')


def test_an_isocyanide_drawn_with_its_own_atoms_is_kept():
    mols = read_cdx_molecule(fixture('US20200237771A1_C00162.cdx'), as_single=False)
    found = [m for m in mols if m is not None and m.HasSubstructMatch(ISOCYANIDE)]
    assert [Chem.MolToSmiles(m) for m in found] == ['[C-]#[N+]CC(=O)OCC']
    assert not any(m is not None and m.HasProp('_CDX_LABEL_REPAIRS') for m in mols)
