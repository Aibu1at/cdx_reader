import json
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from cdx_reader import (
    read_cdx_molecule,
    get_fragment_smiles,
    add_first_attachment_point_number,
    join_molecule_fragments,
    FRAGMENT_LIBRARY,
)

DATA_DIR = Path(__file__).parent / "data" / "US20140045811A1-20140213"
REFERENCES_FILE = Path(__file__).parent / "data" / "references.json"
KNOWN_DIFF_FILE = Path(__file__).parent / "data" / "known_differences.txt"

REFERENCES = json.loads(REFERENCES_FILE.read_text())

_known_diffs = set()
if KNOWN_DIFF_FILE.exists():
    _known_diffs = {
        line.strip() for line in KNOWN_DIFF_FILE.read_text().splitlines()
        if line.strip()
    }


def _canonical(smiles_or_mol):
    """Canonical SMILES without stereo for robust comparison."""
    if isinstance(smiles_or_mol, str):
        mol = Chem.MolFromSmiles(smiles_or_mol)
    else:
        mol = smiles_or_mol
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)


@pytest.mark.parametrize("stem", sorted(REFERENCES.keys()))
def test_read_cdx_molecule_matches_reference(stem, request):
    """Every CDX file should parse into a molecule whose canonical SMILES matches
    the curated reference from cdx_comparison_diff.json. Known-hard cases are
    marked xfail via a sidecar file so only genuine regressions break CI."""
    if stem in _known_diffs:
        request.node.add_marker(
            pytest.mark.xfail(
                reason="known CDX parser limitation (listed in known_differences.txt)",
                strict=False,
            )
        )

    cdx_path = DATA_DIR / f"{stem}.CDX"
    mol_from_cdx = read_cdx_molecule(cdx_path, as_single=True)
    assert mol_from_cdx is not None, f"Failed to parse {cdx_path.name}"

    raw_ref = REFERENCES[stem]
    raw_got = Chem.MolToSmiles(mol_from_cdx)
    got = _canonical(mol_from_cdx)
    ref = _canonical(raw_ref)
    assert ref is not None, (
        f"Invalid reference SMILES for {stem}\n"
        f"  file:  {cdx_path}\n"
        f"  raw:   {raw_ref!r}"
    )

    if got != ref:
        got_mol = Chem.MolFromSmiles(got) if got else None
        ref_mol = Chem.MolFromSmiles(ref)
        got_heavy = got_mol.GetNumHeavyAtoms() if got_mol else -1
        ref_heavy = ref_mol.GetNumHeavyAtoms() if ref_mol else -1
        got_formula = rdMolDescriptors.CalcMolFormula(got_mol) if got_mol else "?"
        ref_formula = rdMolDescriptors.CalcMolFormula(ref_mol) if ref_mol else "?"

        pytest.fail(
            f"{stem}: CDX parse differs from reference\n"
            f"  file:      {cdx_path}\n"
            f"  expected:  {ref}\n"
            f"  got:       {got}\n"
            f"  raw got:   {raw_got}\n"
            f"  raw ref:   {raw_ref}\n"
            f"  heavy atoms (got/ref): {got_heavy}/{ref_heavy}"
            + (f"  (diff={got_heavy - ref_heavy:+d})" if got_heavy >= 0 and ref_heavy >= 0 else "")
            + f"\n"
            f"  formula     (got/ref): {got_formula} / {ref_formula}",
            pytrace=False,
        )


def test_read_cdx_molecule_from_bytes():
    """read_cdx_molecule should accept raw bytes, not just a path."""
    cdx_path = next(DATA_DIR.glob("*.CDX"))
    data = cdx_path.read_bytes()
    mol = read_cdx_molecule(data, as_single=True)
    assert mol is not None
    assert mol.GetNumAtoms() > 0


def test_read_cdx_molecule_as_single_false():
    """as_single=False should return a sequence of molecules."""
    cdx_path = next(DATA_DIR.glob("*.CDX"))
    result = read_cdx_molecule(cdx_path, as_single=False)
    assert result is not None
    mols = list(result)
    assert len(mols) >= 1


def test_fragment_library_lookup():
    """FRAGMENT_LIBRARY should resolve common abbreviations."""
    smi, smi_attach, is_scaffold = get_fragment_smiles("Me")
    assert smi == "C"
    assert "[*:1]" in smi_attach
    assert is_scaffold is False

    smi, _, _ = get_fragment_smiles("Ph")
    assert smi  # non-empty

    assert "Boc" in FRAGMENT_LIBRARY
    assert "Fmoc" in FRAGMENT_LIBRARY


def test_fragment_library_unknown_label():
    """Unknown labels should return a (None, None, False) tuple, not raise."""
    result = get_fragment_smiles("NotARealGroup12345")
    assert isinstance(result, tuple)
    assert len(result) == 3
    assert result[0] is None
    assert result[2] is False


def test_fragment_library_scaffold_label():
    """R-group labels should be classified as scaffold."""
    smi, smi_attach, is_scaffold = get_fragment_smiles("R1")
    assert is_scaffold is True
    assert smi is None


def test_add_first_attachment_point_number():
    assert add_first_attachment_point_number("[*]C", 1) == "[*:1]C"
    assert add_first_attachment_point_number("[*:2]C", 3) == "[*:3]C"


def test_join_molecule_fragments():
    scaffold = Chem.MolFromSmiles("[*:1]c1ccccc1")
    decoration = Chem.MolFromSmiles("[*:1]C")
    joined = join_molecule_fragments(scaffold, decoration)
    assert joined is not None
    expected = Chem.MolFromSmiles("Cc1ccccc1")
    assert _canonical(joined) == _canonical(expected)
