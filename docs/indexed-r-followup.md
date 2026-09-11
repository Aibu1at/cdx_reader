# Indexed letter R labels: follow-up

Completed 2026-09-11, against the previously validated local revision
`4427b4b8533bc3f69c1554b150f0e457fd093626`. The candidate package fingerprint is
`62e244a2b2bee3dc930a489f73d5938fffc8dd792678a7648dddf0cfffa43086`.
The preceding [full-corpus report](validation-report.md) describes the earlier
package, not a full rerun of this follow-up.

## Behavior

Uninterpreted `R1a`, `R1b`, `R2x`, `R27z` and equivalent generic nicknames now
become distinct dummy atoms. Their original styled text is concatenated and
retained in `_CDX_ATOM_LABEL`; `dummyLabel` is the display label. The source
`_CDX_ATOM_ID` and bonds identify the attachment. Indexed letter labels do not
receive `_MolFileRLabel`: `R1a` and `R1b` must not collapse into numeric R1.
Repeated labels and expansion of neighboring abbreviations remain supported.
The public fragment lookup retains its previous literal superscript-digit
coverage, including `R8²`, alongside the newly supported suffixes.

Adding a dummy initially blocked the existing hydrogen repair in three related
files: US20120201816A1 C00154, US20170050964A1 C00151 and US20200140436A1 C00150.
Their neutral nitrogen has two single bonds and two explicit H atoms. The
previous reader already reduced this to one H, but trying RDKit's undefined
default valence (-1) on the new dummy aborted that correction. In this setting,
we now skip dummy atoms and permit only that isolated `NH2 -> NH` correction.
Charged, radical and query atoms and other hydrogen adjustments are excluded
from this additional path. The existing Kekule/query fallback is preserved.

Four real CDX/MOL pairs were added, with source IDs, attachment neighbors and
hashes in [references.json](../tests/data/patent_regressions/references.json):
US20120220583A1 C00002, US20170362217A1 C00007, US20120077856A1 C00007 and
US20200140436A1 C00150. One existing deferred ClZn fixture also contains eight
previously misread R10b/R11b labels. Its old chemical descriptors are retained
as compatibility checks on a comparison copy with only these eight variables
mapped back to their former carbon interpretation. The actual output is checked
separately for dummy identity, complete labels and source neighbors. The original
46 chemical references are unchanged.

## Validation and limits

Selection inspected all 212,917 results from the preceding full run and selected
2,365 with retained `R + digits + lowercase letter` source text. We added all
101 portable controls and all 1,982 previous `None` results. After deduplication:
**4,442 CDX inputs, 4,919 archive/member origins**. Both versions actually reread
these bytes in isolated workers, using the same Python 3.12.3, RDKit 2025.03.3,
pycdxml revision and 30-second limit. Reproducing selection yielded identical
input and origin sets.

| Outcome | Before | After |
| --- | ---: | ---: |
| Molecule returned | 2,458 | 2,461 |
| `None` | 1,982 | 1,979 |
| Exceptions | 2 | 2 |
| Timeout / crash | 0 | 0 |
| Accepted by the existing HARVEST filter | 436 | 79 |

All 2,458 previously returned graphs remain. In 2,365 inputs, **11,260 source
atoms change only from carbon to labeled variable**. Every other stored atom
property and every bond (endpoints, order, stereo) is identical. Component
counts are unchanged. The source audit additionally checked 51,438 labels,
including 17,902 R labels, and 1,095 retained tetrahedral centers.

All 357 accepted-to-rejected transitions are explained by source variables
previously mistaken for ordinary carbon. There are no rejected-to-accepted
transitions. The three recovered graphs are US20210053989A1 C00782/C00802 and
US20140194413A1 C00146. Replaying both valence helpers on the same imported graph
isolates four `NH2 -> NH` changes, each on neutral N with two single bonds and
an explicit ChemDraw invalid-valence warning. Atom identities, connectivity,
charges, isotopes and drawn tetrahedral orientation are preserved. All three
remain scaffolds. Their other legacy label interpretations are not certified as
correct merely because a graph can now be returned.

The final review accounts for all 2,368 changed inputs; **no unexplained flags
remain**. Compact evidence is in [indexed-r-summary.json](indexed-r-summary.json)
and [indexed-r-delta-evidence.json](indexed-r-delta-evidence.json).

- Parser tests: **343 passed, one pre-existing expected failure**; no new xfails.
- Portable corpus: **101 CDX fixtures from 60 patents**, plus the original 46.
- CLI: 98 parsed, three previously documented failures.
- Temporary HARVEST import: enabled and using the candidate; 58 relevant tests
  pass. `parse_cdx_bytes` accepts 79 of 101 fixtures and rejects 22.

The bounded conclusion is **no detected degradation on the checked corpus**.
The focused selection cannot cover indexed labels absent from the previous
selected molecule or previous exception payloads outside the portable controls.
It is not another full 212,917-input run. Intermediate broader valence variants
were rejected during testing; their extra recovered graphs are not final results.
HARVEST's installed parser and production configuration were not switched.

## Explicitly deferred by the user

- Automatic molecule selection and CDX/MOL export arbitration (review point 2):
  a future agent will choose the molecule through a different mechanism.
- A new return contract / InChI policy for aromatic queries (review point 4).
- `Ar` in US20100081610A1 C00062 is a source-document defect: the nodes encode
  disconnected argon although the patent defines an aryl group. Do not globally
  reinterpret explicit argon or invent missing source bonds.

The [previously documented pycdxml limitations](corpus-validation.md#deferred-dependency-and-interpretation-problems)
remain. Development remains local; nothing was published.

## Reproduction

Use the same environment as in [corpus-validation.md](corpus-validation.md).
Keep the earlier reader at `4427b4b` and the candidate in separate checkouts.
From the candidate checkout:

```bash
python scripts/select_indexed_r.py --source-db /data/candidate_final.sqlite \
  --engine candidate_labels6 --output-db /data/indexed-r.sqlite
python scripts/compare_corpus.py run --db /data/indexed-r.sqlite \
  --repo /checkouts/previous --engine before --workers 6 --timeout 30
python scripts/compare_corpus.py run --db /data/indexed-r.sqlite \
  --repo /checkouts/candidate --engine after --workers 6 --timeout 30
python scripts/compare_corpus.py compare --db /data/indexed-r.sqlite \
  --before before --after after --output /data/comparison
python scripts/audit_corpus.py --db /data/indexed-r.sqlite \
  --before-db /data/indexed-r.sqlite --before before --after after \
  --before-repo /checkouts/previous --after-repo /checkouts/candidate \
  --output /data/source-audit
python scripts/audit_indexed_r_delta.py --db /data/indexed-r.sqlite \
  --before before --after after --before-repo /checkouts/previous \
  --after-repo /checkouts/candidate --output /data/delta-evidence.json
python scripts/trace_joins.py --evidence /data/source-audit/evidence.jsonl \
  --before-repo /checkouts/previous --after-repo /checkouts/candidate \
  --output /data/joins.json
python scripts/review_audit.py --evidence /data/source-audit/evidence.jsonl \
  --traces /data/joins.json --delta-evidence /data/delta-evidence.json \
  --output /data/review
```

Local results live in HARVEST's
`results/diagnostics/cdx_regression_20260910/indexed-r-followup/`.
The final engines are `before_suffix` and `after_suffix_narrow` in
`comparison.sqlite`; final outputs are `final-comparison/`,
`final-source-audit/`, `final-delta-audit.json`, and `final-review/`.
