import csv
from pathlib import Path

from cdx_reader.cli import convert_folder, extract_tag

DATA_DIR = Path(__file__).parent / "data" / "US20140045811A1-20140213"


def test_extract_tag():
    assert extract_tag("US20140045811A1-20140213-C00001.CDX") == "C00001"
    assert extract_tag("foo-bar-baz.cdx") == "baz"
    assert extract_tag("single.cdx") == "single"


def test_convert_folder(tmp_path):
    out = tmp_path / "out.csv"
    n_ok = convert_folder(DATA_DIR, out)
    assert out.exists()
    assert n_ok > 0

    with out.open() as f:
        rows = list(csv.DictReader(f))

    # Every CDX fixture should appear in the CSV
    cdx_files = sorted(DATA_DIR.glob("*.CDX"))
    assert len(rows) == len(cdx_files)

    # Headers are correct
    assert set(rows[0].keys()) == {"file_name", "tag", "smiles"}

    # At least one row has a parsed SMILES
    assert any(row["smiles"] for row in rows)

    # Tag extraction matches expectation
    for row in rows:
        assert row["tag"] == extract_tag(row["file_name"])
