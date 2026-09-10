# CDX repairs: validation report

Completed 2026-09-11. Baseline: `452dd4fcf60ce934974625db60fca14c3d6ffefd`.
The validated candidate package fingerprint is
`4c887cdd8909e9694ee429a77d688d95493b6ad245a5048e5cdf88e738a49b98`.
Both versions used Python 3.12.3, RDKit 2025.03.3 and pycdxml 0.1.6b0 at
`0cb554faaa25494b374473852d2e839c91691c73`, with a 30-second per-file limit.

**No degradation was detected on the checked corpus.** This is a bounded
comparison against checked structures and explicit source evidence, not a claim
that every legacy result is chemically correct or that the selected drawing is
the compound intended by a patent paragraph. The existing default selection and
public signatures remain in place. HARVEST's installed package was not changed.

## Changes and representative evidence

| Case | Result and evidence |
| --- | --- |
| US20200323858A1 C00301 | R2/R3/R4/R5 survive on source nodes 51/62/64/66. Styled text, repeated numbers and adjacent abbreviations are checked. |
| Alternate variables | Generic and recognizable uninterpreted labels including G, J, L, X1/Y2/Z3 survive as labeled dummies. Numeric R properties are not invented for other letters. All surviving source text also remains available through `_CDX_ATOM_LABEL`. |
| Source R atoms versus join atoms | Only newly cut dummy atoms connect abbreviations. Source maps and R isotopes cannot be used as temporary connection numbers. Traces cover 98 affected inputs, including 31 old joins consuming source variables and 17 old map collisions; neither occurs in the candidate traces. |
| US20140303156A1 C00009 | The explicit Na+ remains sodium, giving `C10H12N3NaO4`, rather than nitrogen inferred from the label prefix. |
| US20180291030A1 C00033 | The valid nested nitro graph takes precedence over the contradictory outer node type. Formula `C21H19Cl2N3O5`, full InChIKey `KCSNIJFTSPKGEI-XJRBWHPUSA-N`; source stereochemistry is retained. |
| US20140303156A1 C00012 | Only the warned neutral N with two single bonds and contradictory zero H is repaired. Formula `C16H17FN4O3`, no radical electrons; the patent's Example 2 mass data support the explicit graph. |
| US20130029938A1 C00060 | An abbreviation at atom index zero now expands correctly to dimethyl 3-aminopentanedioate: `COC(=O)CC(N)CC(=O)OC`. The intermediate title and six-proton ester signal support the change; the paired MOL also had the wrong ketone. Fmoc and Boc examples exercise the same boundary. |
| US20120225858A1 C00010/C00011 | Stereo cleanup must not erase drawn tetrahedral orientation just because R3 and R4 appear identical as RDKit wildcard atoms. Source-neighbor parity checks protect these centres. |
| US20140005221A1 C00039/C00040 | Explicitly aromatic edges to variable atoms are retained even when no unique Kekule form or InChI exists. The returned scaffold carries `_CDX_UNKEKULIZED_QUERY`; it does not pass the complete-molecule filter. |
| US20100130737A1 C00039 | Added a control for the reported MOL amide error: CDX yields `*NC(C)=O` with R3, while MOL yields `[3*]CC(C)=O`. The existing CDX nitrogen interpretation is preserved. |

In the checked baseline, Python's `parse_cdxml` already concatenated styled
segments. The R-number problem was therefore not solely a first-`<s>` bug in
that function: full labels were not preserved on the returned RDKit atoms, and
source R numbers could interfere with abbreviation joins. Metadata now follows
`_CDX_ATOM_ID` through expansion, rather than an atom index.

Each portable fixture has CDX/MOL SHA-256, original archive member, expected
properties and an evidence explanation in
[references.json](../tests/data/patent_regressions/references.json).
[Variable-label evidence](variable-label-evidence.json) also distinguishes actual
variables from compound identifiers encoded as atom labels. Neither is a carbon
atom merely because the previous reader interpreted it that way.

## Corpus and before/after results

Inventory covered 970 local ZIP archives spanning 2001–2025: 240,093 CDX member
origins, 212,917 unique byte sequences and 838 patent identifiers with CDX.
Deduplication used content SHA-256 and retained every archive/member origin.
Different archive versions were not merged by patent identifier. Both versions
actually reread every unique CDX; historical HARVEST exports were not baselines.

| Outcome | Before | After |
| --- | ---: | ---: |
| RDKit molecule returned | 200,794 | 200,914 |
| Returned `None` | 2,102 | 1,982 |
| Exceptions | 10,021 | 10,021 |
| Timeouts or process crashes | 0 | 0 |
| Accepted by HARVEST's complete-molecule filter | 170,187 | 168,613 |

There are 120 recovered graphs and no previously returned graph lost to an
exception, timeout or `None`. Of all inputs, 13,230 are unchanged, 186,921 change
only metadata, and 12,766 change chemical descriptors or parser status. Every
one of the 199,687 changed inputs was reopened for source inspection.

The audit checked 2,568,138 retained source labels, including 82,634 R labels,
and compared the orientation of 15,494 previously specified tetrahedral centres
in chemically changed results. It also checked retained nested atoms, charges,
isotopes, internal and external bonds, prior bond stereo, formula, full InChIKey,
isomeric SMILES, components, radicals and HARVEST acceptance. No source-label,
nested-graph or retained-stereo mismatch remained.

The final classifications are 186,921 metadata changes and 12,766 source
corrections. No review flags remain unresolved. These classifications concern
the changed source feature; an unchanged ambiguity elsewhere in the same drawing
does not thereby become a chemically validated structure.

### Acceptance transitions and all flagged difference classes

| Flag or transition | Inputs | Disposition |
| --- | ---: | --- |
| Accepted → rejected | 1,583 | Every result contains a retained, source-identified placeholder explaining rejection. These include real variables and compound-reference labels previously read as ordinary atoms. |
| Rejected → accepted | 9 | Successful abbreviation joins restore complete selected structures. Source attachment traces support the changes; this does not establish which compound a multi-structure drawing should represent. |
| Previous MOL agreement lost | 342 | Explicit source graphs/labels support the change. MOL agreement alone is not an override. First-atom expansion cases additionally require a successful join to a generated cut atom. |
| Component count changed | 37 | Traces isolate corrected abbreviation joins or source-map collisions. The selection algorithm was not rewritten. |
| Total radical count increased | 5 | All additional radicals occur in additionally retained source components. The baseline assigns exactly the same radicals before its erroneous abbreviation join discards those components; no radical count increased on previously retained atoms. |

Flags overlap. The 1,643 flagged inputs are all accounted for by these classes.
The five radical cases are US20140088079A1 C00022, US20200062798A1
C00033/C00034/C00039 and US20210047311A1 C00073. Their exact atom IDs, baseline
pre-expansion radicals and join traces are in [r-join-evidence.json](r-join-evidence.json).
They remain scaffold/legacy-chemistry limitations, not examples of repaired
radical chemistry.

The nine acceptance gains are US20240336629A1
C00552/C00562/C00588/C00601/C00606/C00607, US20130029938A1 C00061,
US20140024636A1 C00224 and US20120208833A1 C00010. They cover CF3 at the first
atom, ester expansion and Cbz expansion. A correct local expansion does not
change the separate responsibility of choosing the relevant compound. Their
final SMILES, full InChIKeys and input hashes are in
[accepted-gains.json](accepted-gains.json).

## Tests and integration

- Original baseline: 54 passing tests and one existing expected failure.
- Final suite: **318 passed, one existing xfail**. No new xfails and no changes
  to the original 46 chemical references. The additional portable corpus has
  **97 distinct binary CDX fixtures from 56 patents**, with paired MOL files,
  plus synthetic boundary cases.
- Public path, bytes and list input modes, metadata-to-atom correspondence,
  repeated R numbers, R-number MOL round trips, isotope/stereo cases and CLI are
  exercised. Plain SMILES cannot serialize all label metadata.
- CLI on all 97 new fixtures: 94 parsed, 3 known failures (two dependency errors
  and one deferred `None` case).
- Temporary HARVEST import: available and resolved to the candidate checkout.
  Its relevant 58 unit tests passed; `parse_cdx_bytes` accepted 79 of the 97
  fixtures and rejected 18. No package installation or production switch occurred.

## Deferred issues and limits

The same 10,021 exception payloads remain. They were separately reread with the
final descriptor handling and compared exactly: zero differences.

| Existing exception class | Inputs |
| --- | ---: |
| Legacy CDX document header | 5,412 |
| GraphicType decoder requires two bytes | 2,211 |
| No RDKit molecules; empty `reduce` | 1,296 |
| INT8 decoder size error | 807 |
| Unknown object tag 32772, including unsupported GenericList | 196 |
| ExternalConnectionType decoder requires two bytes | 72 |
| Unknown object tag 32775 | 16 |
| Invalid XML text | 6 |
| RDKit duplicate-bond error | 4 |
| Atom-map number outside supported range | 1 |

[`ExternalConnectionType` is specified as INT8](https://chemapps.stolaf.edu/iupac/cdx/sdk/Node.htm), but the pinned pycdxml decoder
requires two bytes. This is an upstream binary-decoding problem. No properties
were discarded and no alternative reader was added as a production dependency.

The separate claim that every missing `ExternalConnectionPoint` is a lost
substituent attachment was not confirmed. In two checked archive versions,
465 decoded points are inside groups whose parents have one, two or three outer
bonds; none has an unattached parent. Such points normally disappear when the
internal and outer graph are connected. Source atoms and resulting bonds, rather
than dummy counts alone, establish whether a group survived. Archive fingerprints,
counts and seven-point examples are in
[external-connection-evidence.json](external-connection-evidence.json). Decoding
failures and ambiguous nested groups still limit this assessment.

Other deferred work includes ambiguous U/V/K and D/T labels, low-valent Sn,
warned K coordination, invalid halogen/metal attachments, multi-attachment or
recursive nested groups, mixed drawings where nested promotion would disturb
selection, and unrelated valence errors exposed by valid nested expansion. The
latter retry the legacy nested treatment instead of resetting all hydrogens.
The original in-ring abbreviation limitation remains. Explicit extraction of
all labeled objects/alternative groups belongs to a future API.

## Artifacts and reproduction

[corpus-validation.md](corpus-validation.md) gives the inventory, isolated reader,
comparison, source audit, join tracing and review commands. Compact summaries and
environment fingerprints are in [validation-summary.json](validation-summary.json).
Portable test inputs occupy about 1.3 MiB; full databases and JSONL are local only.

The local artifact root is the HARVEST checkout's
`results/diagnostics/cdx_regression_20260910/`:

- `baseline.sqlite`, engine `baseline2`, and `candidate_final.sqlite`, engine
  `candidate_labels6`: complete results and every input origin;
- `baseline-errors-only.sqlite`, engine `baseline_descriptor_check`: the
  separate repeat of all original exceptions;
- `final-comparison/`: complete before/after differences and environment;
- `final-source-audit/`: per-input source evidence and initial review flags;
- `final-review/`: final classifications and all 1,643 reviewed cases;
- `cli_final.csv` and `harvest_integration_final.json`: CLI/integration outputs.

Intermediate engines were diagnostic runs and are not the final comparison.
