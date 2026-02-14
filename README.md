# cdx-reader

A Python library for reading ChemDraw CDX/CDXML files and converting them to RDKit molecule objects, with automatic expansion of fragment abbreviations.

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

```python
from cdx_reader import read_cdx_molecule

# From a file path
mol = read_cdx_molecule("structure.cdx")

# From bytes (e.g. read from a database or zip archive)
with open("structure.cdx", "rb") as f:
    mol = read_cdx_molecule(f.read())

# Without fragment expansion
mol = read_cdx_molecule("structure.cdx", expand_fragments=False)

# Get individual fragments instead of a combined molecule
mols = read_cdx_molecule("structure.cdx", as_single=False)
```

## API Reference

### `read_cdx_molecule(path_or_bytes, as_single=True, expand_fragments=True)`

Read a CDX file and return an RDKit `Mol` object.

**Parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `path_or_bytes` | `str`, `Path`, or `bytes` | required | Path to a CDX file or raw CDX bytes |
| `as_single` | `bool` | `True` | Combine all fragments into a single molecule |
| `expand_fragments` | `bool` | `True` | Expand multi-atom abbreviations (Me, Ph, Boc, etc.) |

**Returns:** An RDKit `Mol` object (if `as_single=True`) or a tuple of `Mol` objects.

### `get_fragment_smiles(label)`

Look up a text label in the fragment library.

**Parameters:**

| Parameter | Type | Description |
|---|---|---|
| `label` | `str` | Text label from a ChemDraw node (e.g. "Me", "CF3", "R1") |

**Returns:** `(smiles, smiles_with_attachment, is_scaffold)` tuple. Returns `(None, None, True)` for scaffold/R-group labels and `(None, None, False)` for unknown labels.

### `join_molecule_fragments(scaffold, decoration, regenerate_2d_coords=False)`

Join two RDKit molecules at labeled attachment points (`[*:N]`).

**Parameters:**

| Parameter | Type | Default | Description |
|---|---|---|---|
| `scaffold` | `Mol` | required | RDKit molecule with a labeled attachment point |
| `decoration` | `Mol` | required | RDKit molecule with a matching labeled attachment point |
| `regenerate_2d_coords` | `bool` | `False` | Recompute 2D coordinates after joining |

**Returns:** Joined `Mol` object, or `None` on failure.

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

## License

MIT
