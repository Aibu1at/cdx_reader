# cdx-reader

A Python library (and CLI) for reading ChemDraw CDX/CDXML files and converting them to RDKit molecule objects, with automatic expansion of fragment abbreviations.

## Problem

ChemDraw CDX files frequently contain abbreviated chemical groups (Me, Ph, Boc, CF3, etc.) represented as text labels rather than fully drawn atoms. Standard parsers either ignore these labels or produce molecules with incorrect atom types. `cdx-reader` fixes atom types in the CDXML tree and optionally expands multi-atom abbreviations into proper molecular fragments.

## Features

- **CDX/CDXML parsing** via `pycdxml` with automatic atom type correction
- **Fragment expansion**: 80+ common abbreviation labels (Me, Et, iPr, tBu, Ph, Bn, OMe, CF3, Boc, Fmoc, etc.) expanded into full molecular structure
- **Scaffold detection**: R-group labels (R1, R2, X, Y, ...) recognized and mapped to dummy atoms
- **Valence fixing** post-processing for sanitized RDKit molecules
- **Accepts file paths or raw bytes** as input

## Installation

```bash
pip install cdx-reader
```

### Dependencies

- [RDKit](https://www.rdkit.org/) (>=2022.03)
- [pycdxml](https://github.com/kienerj/pycdxml)

## Quick Start

### Python API

```python
from cdx_reader import read_cdx_molecule
from rdkit import Chem

# From a file path
mol = read_cdx_molecule("structure.cdx")
print(Chem.MolToSmiles(mol))

# From bytes (e.g. read from a database or zip archive)
with open("structure.cdx", "rb") as f:
    mol = read_cdx_molecule(f.read())

# Get individual fragments instead of a combined molecule
mols = read_cdx_molecule("structure.cdx", as_single=False)
```

Multi-atom abbreviations (Me, Ph, Boc, CF3, ...) are always expanded into
their full structures.

### CLI: batch-convert a folder of CDX files

```bash
cdx-to-csv path/to/cdx_folder/ -o smiles.csv
```

Output CSV has three columns: `file_name`, `tag`, `smiles`. The `tag` is
the last dash-separated component of the file name stem — for example
`US20140045811A1-20140213-C00001.CDX` → `C00001`. Files that fail to
parse still appear in the CSV with an empty `smiles` column, so you can
audit which inputs failed.

Options:
- `--pattern '*.cdx'` — case-sensitive glob (default: `*.CDX`)
- `-r` / `--recursive` — recurse into subdirectories
- `-o` / `--output` — output CSV path (default: `cdx_smiles.csv`)


## Supported Fragment Abbreviations

<details>
<summary>Click to expand full list</summary>

**Alkyl groups:** Me, Et, Pr, iPr, nPr, Bu, nBu, iBu, tBu, sBu

**Aromatic groups:** Ph, Bn, Bz

**Oxygen-containing:** MeO, OMe, PhO, EtO, OEt, Ac, AcO, OAc, OCF3

**Nitrogen-containing:** MeN, NMe, NMe2, NEt2, NO2, CN

**Halogenated:** CF3, F3C, CCl3, CHF2, CH2F

**Sulfur-containing:** SMe, SO2Me, MeSO2, SO2NH2, SO3H

**Carbonyl groups:** CHO, CO2H, COOH, CO2Me, CO2Et, CONH2, CONMe2

**Protecting groups:** Boc, Cbz, Fmoc, BocHN, CbzHN, FmocHN, NHBoc, NHCbz, NHAc, NHBn, AcHN, TfaHN

**Sulfonamides:** TsHN, MsHN, NHTs, NHMs

**Other:** B(OH)2, BnO2C, BnO2CHN

</details>

## Limitations

**Abbreviations embedded inside a ring cannot be expanded.** Fragment
expansion works by splitting the molecule at the labeled atom into two
pieces — the scaffold and the abbreviation — and then re-joining them
with the expanded structure. When the labeled atom sits inside a ring,
cutting that single bond does not actually separate the graph into two
components, so the splitter has nothing to replace and the abbreviation
is left as an unexpanded heteroatom (or dummy) in the ring.

Example: [`US20140045811A1-20140213-C00235`](tests/data/US20140045811A1-20140213/US20140045811A1-20140213-C00235.CDX)
contains a ring with an `Me`-labeled atom that participates in two ring
bonds — the parser emits the ring but cannot substitute the methyl in
place, so the result differs from the reference structure.

## License

MIT
