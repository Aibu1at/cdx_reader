# Patent regression fixtures

These original binary CDX files and paired MOL files supplement, without replacing,
the existing 46 reference structures. `references.json` records the original patent
and archive member, SHA-256 of both files, expected properties, coverage tags and
the evidence used to choose each expectation. The files are USPTO patent drawings.

The initial diverse controls were selected from independent patents across
2004–2025. For these controls, direct conversion of the unmodified CDXML, the
paired MOL, and the original reader agreed on the complete isomeric graph. They
include tetrahedral and double-bond stereochemistry and pre-expanded groups.
Additional fixtures reproduce the diagnosed defects and cases discovered in the
large local comparison, including isotope and charged salt drawings.

Evidence categories are intentionally distinct:

- `reviewed_diagnostic`: graph checked against explicit CDX properties and the
  named intermediate/product or table entry in the patent; the paired MOL may
  itself be wrong, as explained in that entry.
- `source_and_mol` / independent controls: source graph and paired MOL agree,
  including isomeric SMILES. MOL agreement alone is not used to resolve an
  inconsistent drawing.
- `source_and_patent`: a suspected source typo is checked against the relevant
  compound name and surrounding patent text.
- R-label controls: verify the original node IDs and displayed labels; do not
  claim that the remaining Markush variables form a resolved molecule.
- `deferred_behavior` / `deferred_error`: compatibility expectations for an
  explicitly documented limitation. These are **not chemical gold standards**.
  Known errors are asserted as errors rather than hidden with new `xfail` marks.

The real-file tests verify hashes, formula, full InChIKey and isomeric SMILES,
preserving isotopes. R-label checks also verify attachment to original CDX IDs,
absence of artificial atom maps, existing list mode and MOL R-number round trips.
Small synthetic CDXML tests exercise separately styled R numbers, repeated R
numbers, adjacent abbreviation expansion and the boundaries of chemical repairs.

Full-corpus SQLite databases and before/after JSONL results are intentionally not
shipped here; they contain hundreds of thousands of records. See
[corpus-validation.md](../../../docs/corpus-validation.md) for reproduction.
