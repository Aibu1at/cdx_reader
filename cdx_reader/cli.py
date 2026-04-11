"""CLI: convert all CDX files in a folder to a SMILES CSV."""
import argparse
import csv
import logging
from pathlib import Path

from rdkit import Chem

from .cdx_reader import read_cdx_molecule

logger = logging.getLogger(__name__)


def extract_tag(file_name: str) -> str:
    """Tag = last dash-separated component of the file stem.

    Example: 'US20140045811A1-20140213-C00001' -> 'C00001'
    """
    stem = Path(file_name).stem
    return stem.split('-')[-1]


def convert_folder(input_dir: Path, output_csv: Path, pattern: str = '*.CDX',
                   recursive: bool = False) -> int:
    """Walk input_dir for CDX files and write (file_name, tag, smiles) CSV.

    Returns number of rows successfully parsed.
    """
    input_dir = Path(input_dir)
    output_csv = Path(output_csv)
    glob = input_dir.rglob(pattern) if recursive else input_dir.glob(pattern)
    files = sorted(glob)
    n_ok = 0
    n_fail = 0

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['file_name', 'tag', 'smiles'])
        for path in files:
            try:
                mol = read_cdx_molecule(path, as_single=True)
                if mol is None:
                    n_fail += 1
                    smiles = ''
                else:
                    smiles = Chem.MolToSmiles(mol)
                    n_ok += 1
            except Exception as e:
                logger.warning("Failed to parse %s: %s", path.name, e)
                n_fail += 1
                smiles = ''
            writer.writerow([path.name, extract_tag(path.name), smiles])

    print(f"Wrote {output_csv}: {n_ok} parsed, {n_fail} failed, {len(files)} total")
    return n_ok


def main():
    parser = argparse.ArgumentParser(
        description='Convert all CDX files in a folder to a SMILES CSV.')
    parser.add_argument('input_dir', type=Path,
                        help='Directory containing .CDX files')
    parser.add_argument('-o', '--output', type=Path, default=Path('cdx_smiles.csv'),
                        help='Output CSV path (default: cdx_smiles.csv)')
    parser.add_argument('--pattern', default='*.CDX',
                        help='Glob pattern for CDX files (default: *.CDX)')
    parser.add_argument('-r', '--recursive', action='store_true',
                        help='Recurse into subdirectories')
    args = parser.parse_args()

    logging.basicConfig(level=logging.WARNING)
    convert_folder(args.input_dir, args.output, args.pattern, args.recursive)


if __name__ == '__main__':
    main()
