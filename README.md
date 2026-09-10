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

- [RDKit](https://www.rdkit.org/) (2025.03.3, pinned in `pyproject.toml`)
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

The reader attempts to expand supported multi-atom abbreviations (Me, Ph, Boc,
CF3, ...). Ambiguous attachments can prevent expansion; see the limitations below.

### R-group metadata and conservative repairs

Styled R labels are read across all text segments. A surviving R dummy carries
`_CDX_ATOM_ID` (the original node), `_CDX_ATOM_LABEL` (original text), and
`dummyLabel` (the normalized display label). Numeric `Rn` labels additionally
carry the integer `_MolFileRLabel`. Repeated R numbers are allowed; these are
chemical labels, not the atom-map numbers used internally to join fragments.
Ordinary SMILES does not serialize all these properties: retain the RDKit object
when a downstream agent needs label-to-atom correspondence.

Variables are not restricted to the letter R. Declared `GenericNickname` nodes
and recognizable uninterpreted variables such as `G`, `J`, `L`, `X1`, `Y2` and
`Z3` retain their full labels as dummies. Known dictionary abbreviations are
expanded first. An ordinary element label is not made into a variable merely
because it is a letter. Every surviving labeled source atom also retains its
original `_CDX_NODE_TYPE` and `_CDX_ATOM_LABEL`, including ambiguous labels whose
chemical interpretation is deferred. Only numeric **R** labels receive
`_MolFileRLabel`; other labels are not assigned invented R numbers.

Source dummy atom maps are saved in `_CDX_ORIGINAL_ATOM_MAP` and removed from the
joining namespace. Only newly created cut atoms identify abbreviation joins;
source R isotopes and repeated numbers never serve that purpose.

Explicit element numbers, charges and isotopes take precedence over text-prefix
guesses such as interpreting Na as nitrogen. Inconsistent variable letters and
warned, chemically impossible metal/isotope drawings retain the legacy behavior.
The repair of uninterpreted HN labels remains separate.

A valid, connected nested graph with one consistent external attachment is used
before abbreviation expansion, including incorrectly marked `Unspecified` nodes.
Query graphs, multiple attachments, impossible valences and covalent graphs that
RDKit rewrites into dative bonds are not promoted. A narrowly targeted zero-H
repair applies to neutral N with two single bonds, an invalid-valence warning,
and no explicit radical/query annotation. Other hydrogen counts are retained.

An explicitly aromatic bond to a variable may have no unique Kekule form. When
that is the sole sanitization obstacle, the original aromatic query graph is
retained without Kekulization. It remains a scaffold, unsuitable for a complete
molecule filter or for assuming a unique InChI. This fallback never applies to
ordinary complete molecules; the returned object sets `_CDX_UNKEKULIZED_QUERY`.
Drawn tetrahedral tags are retained when stereo is recomputed, including centres
whose distinct variable labels appear identical to RDKit as wildcard atoms.

The existing molecule-selection behavior is preserved. In drawings with multiple
structures and unresolved dictionary abbreviations, nested repairs are deferred:
they can otherwise change which component the existing expander selects. An API
that explicitly returns every labeled object for a future agent is separate work.
Fixing an erroneous join through a source R atom can nevertheless change the
components returned by the old expander: an accidental attachment is not retained
as a compatibility requirement. These differences are audited separately.

Regression data and source evidence are documented in
[tests/data/patent_regressions/README.md](tests/data/patent_regressions/README.md).
For reproducible corpus comparisons see [docs/corpus-validation.md](docs/corpus-validation.md).
The completed 212,917-input comparison, acceptance transitions and remaining
limitations are in [docs/validation-report.md](docs/validation-report.md).

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
