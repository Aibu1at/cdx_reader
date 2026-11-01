# Crossed attachment stems, reviewed 2026-09-25

The marker has two opposite wavy tails and two opposite ordinary stem bonds.
Exactly one ordinary end is a bare terminal node; the other belongs to the
substituent. Turn the intersection into an anonymous port and remove the two
wavy tails and the outer stub. Retain the intersection's original CDX ID.
Both-bare, both-connected, non-collinear and labelled/charged marker nodes are
outside this rule. In particular this is not a rule for any wavy or dashed bond.

`manifest.json` hashes the unmodified CDX fixtures. `artwork_manifest.json` hashes
the original review images. Google Patents images were downloaded using the
publication's own image links; the TIFs are unmodified USPTO APPDT members.

## Positive fixtures

- US20150152119A1 C00066 is pyridin-3-yl with a crossed stem, not tert-butylpyridine.
- US20140227302A1 C00035 is a two-carbon chain leading to pyrazole N; the crossed
  stem is outside that chain. The prepared graph is `*CCn1cccn1`.
- US20160060265A1 C00068 is a para-cyclopropylcarbonylphenyl fragment, with the
  crossed stem on the opposite side of the phenyl ring: `*c1ccc(C(=O)C2CC2)cc1`.
- In US20150152119A1, Table 6/p-0373 selects C00069 and C00075 for core C00065.
  The assembled graphs agree with separately drawn products C00030 (25,
  p-0260–0262) and C00032 (35, p-0266–0268), respectively. Full keys:
  `CSHAIEWAKNMFDS-UHFFFAOYSA-N` / `DVKGEUMVNRWJTP-UHFFFAOYSA-N`; formulas
  `C16H17FN4` / `C17H18FN3O`. Ordinary explicit N–H is representational.
  These checks preserve the stereo information present in the read graphs; they
  do not establish an absolute stereochemical assignment absent from those graphs.
  The printed low-resolution ESI values are 285.7 / 300.4, calculated [M+H]+
  285.1510 / 300.1507. They are supporting evidence, not exact-mass gold values.

## Five legacy references

The immutable `patent_regressions/references.json` agrees with the old reader and
MOL export but conflicts with the actual drawn attachment symbols. It remains
unchanged. New source-addressed expectations are in `reviewed_legacy_cases.json`.

| Drawing | Original image | Corrected graph / scope |
|---|---|---|
| US20090018114A1 C00007 | Two crossed bars, after S and at the other ring substituent | `*SC1=C(*)CCOC1` |
| US20090018114A1 C00009 | One three-bond bar and one crossed bar on an unsaturated linker | `*/C=c1/cccc/c1=C/*`; existing double-bond geometry retained |
| US20240051972A1 C00018 | Three methylene vertices between two crossed bars | `*CCC*`, not an alkane with two tert-butyl ends |
| US20250042847A1 C00020 | Three attachment bars on an amine linker | `*CN(C*)CCN*` |
| US20120295881A1 C01043 | Crossed bar at the pyrazole ring | `*c1cnn(C)c1`; the separate, pre-existing interpretation of the misprinted Na ring node is unchanged |

The older C01043 evidence incorrectly called the marker a tert-butyl group. This
review corrects that claim without rewriting the frozen reference. The drawing
really does print Na in the ring; this change does not independently validate or
broaden the legacy Na repair.

`crossed_reference.prior_comparison_copy` first asserts the actual corrected graph,
port IDs, removed nodes, degree, and repair provenance. Only on a copy does it
restore the former carbon representation for the previous regression checks.
`r7_reference.legacy_comparison_copy` then checks the earlier three-bond repairs.
This checks the untouched remainder against frozen answers; the reconstructed
copy is neither an assembly input nor a chemical oracle.

The original US20120035134A1 C00006 explicitly equates -NR2/-CH2NR2/-OR/-CH2OR with
four marked fragments. Two bars have crossed stems. The four wavy bonds left by
the earlier R7 rule are those tails, not stereochemical bonds. Its test now names
all four source-addressed repairs. Real internal wavy bonds remain separately
tested using US20100119475A1 C00029; four detached ambiguous crossed markers in
C00008 remain unchanged.

The local before/after audit covers 9,527 unique CDX inputs from 48 patent archives,
identical RDKit 2025.03.3, separate worker processes and a 30-second file limit.
146 changed inputs match exactly the recorded crossed-bar edits; 9,381 are byte-
identical in the serialized chemical results. All retained raw CDXML attributes,
atoms, labels and stereo outside the edit match; no new exceptions or timeouts.
This is absence of observed degradation on this corpus, not a universal claim.
