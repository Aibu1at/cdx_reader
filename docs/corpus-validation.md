# Reproducing the CDX comparison

Use one Python environment for both source revisions. The verified environment is
Python 3.12.3, RDKit 2025.03.3 and pycdxml 0.1.6b0 at commit
`0cb554faaa25494b374473852d2e839c91691c73`. The production dependency declaration
is unchanged; no alternative reader is required or installed by these tools.

Keep two checkouts, one at baseline `452dd4fcf60ce934974625db60fca14c3d6ffefd`
and one at the candidate revision. Run the tools from the candidate checkout:

```bash
python scripts/compare_corpus.py inventory --db /data/cdx-audit.sqlite /data/patent-zips
python scripts/compare_corpus.py run --db /data/cdx-audit.sqlite \
  --repo /checkouts/cdx-baseline --engine baseline --workers 6 --timeout 30
python scripts/compare_corpus.py run --db /data/cdx-audit.sqlite \
  --repo /checkouts/cdx-candidate --engine candidate --workers 6 --timeout 30
python scripts/compare_corpus.py compare --db /data/cdx-audit.sqlite \
  --before baseline --after candidate --output /data/cdx-comparison
python scripts/audit_corpus.py --db /data/cdx-audit.sqlite \
  --before-db /data/cdx-audit.sqlite --before baseline --after candidate \
  --before-repo /checkouts/cdx-baseline --after-repo /checkouts/cdx-candidate \
  --output /data/cdx-source-audit
python scripts/trace_joins.py --evidence /data/cdx-source-audit/evidence.jsonl \
  --before-repo /checkouts/cdx-baseline --after-repo /checkouts/cdx-candidate \
  --output /data/cdx-joins.json
python scripts/review_audit.py --evidence /data/cdx-source-audit/evidence.jsonl \
  --traces /data/cdx-joins.json --output /data/cdx-review
```

Inventory accepts multiple directories of ZIP files, reads archive members without
bulk extraction, hashes CDX content and retains every ZIP/member origin. Identical
bytes share one parse result. Different archive versions remain separate origins;
different CDX content is never merged by patent number. Changed inputs during
re-inventory or parsing are errors. Inventory failures must be resolved first.

Each parse runs in a worker process, with the same per-file timeout on both sides.
Crashes and timeouts are recorded, and the worker is replaced. Source snapshots
prevent edits during a run from changing later results. Repeating a run resumes
missing files; a changed source fingerprint, dependency version or timeout requires
a new engine name. Completed rows are committed in small batches. Do not run two
copies of the same engine concurrently. `--limit` supports a smoke test followed
by an unrestricted resume.

The comparison includes isomeric and non-isomeric SMILES, full InChIKey, formula,
charge, isotopes, radicals, components, R metadata and retained source atom IDs.
Atom and bond details are stored in each result for source-level assessment.
`harvest_accepted` mirrors the current HARVEST filter: a successful result without
atomic numbers 0, 99 or 100. It is an acceptance test, not proof of chemical identity.
Different dependency versions or incomplete paired runs fail comparison explicitly.
When a retained aromatic query has no unique InChI, the result records
`inchikey_error` separately and keeps its graph, formula and parser status. The
RDKit molecule marks this case with `_CDX_UNKEKULIZED_QUERY`. Descriptor failure
does not mean that CDX decoding failed. Measurements preserve the same SMILES
then InChI order on both revisions.

The source audit reopens every changed CDX, verifies all retained text against original IDs,
checks retained atoms and bonds of repaired nested groups, and records the rule
responsible for each changed node. It flags lost resolution, changed component
counts, added radicals, changed pre-existing bond stereo, disagreement with a
previously matching MOL and unclassified changes for review. A local source-graph
check proves fidelity to the explicit drawing; patent text is still needed where
the drawing itself is inconsistent. A full comparison has to be completed before
interpreting this audit as a complete assessment.

Review `changes.jsonl`, `evidence.jsonl` and `review.jsonl`, not only success counts.
The final review gate rejects unexplained losses, stereo changes and source
mismatches. Join traces distinguish consuming a source variable from consuming a
temporary cut atom. A known source correction does not excuse an unrelated flag.
No result from a previous HARVEST run substitutes for rerunning the CDX reader.
The final claim is limited to **no detected degradation on the checked corpus**.

The completed run is described in [validation-report.md](validation-report.md).
Its baseline errors were also reread with the final descriptor implementation:
all 10,021 error payloads were identical. Keeping query InChI failures separate
from parser failures introduced no change to those baseline results.

For a temporary HARVEST integration check, run from its checkout without
installing the candidate into its environment:

```bash
PYTHONPATH=/checkouts/cdx-candidate .venv/bin/python -c \
  'import cdx_reader; from patent_processor.cdx_extractor import CDX_PARSING_AVAILABLE; print(cdx_reader.__file__, CDX_PARSING_AVAILABLE)'
PYTHONPATH=/checkouts/cdx-candidate .venv/bin/python -m pytest \
  tests/unit/patent_processor/test_zip_source_cdx.py tests/unit/bindingdb/test_helpers.py -q
```

Verify that the printed module path belongs to the candidate and availability is
`True`; a silently disabled CDX import is not a passing integration result.

## Deferred dependency and interpretation problems

- pycdxml does not implement `CDXGenericList` (property `0x0441`); its outer
  exception may misleadingly report unknown object tag 32772.
- Some ChemDraw versions wrote a one-byte `GraphicType`. The pinned pycdxml
  decoder requires two bytes. No properties are dropped to bypass this error.
- The pinned pycdxml `CDXExternalConnectionType.from_bytes` requires two bytes,
  whereas ChemDraw specifies an INT8 property. This causes separate decoding
  failures before any molecule can be returned; it is not lost attachment
  metadata in successfully decoded molecules.
- U/V/K and multi-bond D/T can be patent variables despite misleading element
  numbers. Their ambiguous elemental encoding retains the previous behavior;
  their original labels are preserved for a future agent. Explicit generic and
  recognizable uninterpreted variables, including letters other than R, are
  preserved as placeholders.
- Invalid nested attachments, including halogen-attached ClZn/BrZn, and warned
  K+ coordination drawings need separate chemistry handling.
- Nested repairs in mixed multi-structure drawings with unresolved abbreviations
  are deferred to preserve current selection. Returning all objects and labels
  is a future API, not a parser repair in this change.
- A valid nested group can expose an unrelated surrounding valence error that
  the existing general valence fixer cannot handle. The reader retries without
  nested promotion in that case. Low-valent Sn and warned K coordination also
  remain deferred rather than introducing RDKit radical assignments as fixes.
- Some Markush rings contain explicitly aromatic bonds to variable atoms. Their
  graph can be retained without a unique Kekule form; it is not a complete
  molecule and must not be judged by a complete-molecule SMILES round trip.

`ExternalConnectionPoint` normally connects a nested nickname or alternative
group to its parent. Its disappearance after connecting the parent and internal
atoms is expected, as specified in the
[ChemDraw NodeType documentation](https://chemapps.stolaf.edu/iupac/cdx/sdk/properties/Node_Type.htm).
A count of missing dummy atoms alone does not establish lost substituent
attachment points. The regression suite checks the resulting bond and source
atom IDs. Standalone alternative-group extraction remains separate work.

The existing in-ring abbreviation limitation and its original expected-failure
test remain unchanged.
