# Original-source review of six historical references

Reviewed 2026-09-24. These six cases in `../patent_regressions/references.json`
were previously called complete molecules because the unmodified CDX/RDKit result
agreed with the MOL export. Original USPTO TIF artwork contradicts that claim:
the apparent isopropyl groups are wavy attachment bars.

The original reference JSON, CDX and MOL files remain byte-for-byte unchanged.
`reviewed_legacy_cases.json` is a separate, source-addressed amendment, with original
TIFs and verbatim XML context excerpts. Downloaded CDX bytes were checked against
the hashes of the existing fixtures before associating the TIFs with them.

| Drawing | Source evidence | Scope of the new check |
|---|---|---|
| US20040001895A1 C00017 | paragraph P-0032 defines a two-ended Y linker; TIF shows two bars on adjacent cyclohexane positions | two anonymous ports; cyclohexane-1,2-diyl connectivity |
| US20040001895A1 C00038 | same Y alternative in the claims; original artwork shows the same two bars | same two ports, separate original CDX retained |
| US20090018114A1 C00009 | p-0065/0066 defines an A fragment; the lower drawn bar has the strict R7 topology | **partial repair only**; the other junction has more than one ordinary continuation and is outside the rule; the remainder is not a verified chemical gold structure |
| US20130005755A1 C00025 | Table A, original drawing has a bar at the end of the ethyl chain | `*CCc1ccc(-c2cc[nH]c(=O)c2)cc1` |
| US20160002213A1 C00040 | p-0504/0505 defines R1 option (ii-1); TIF places the bar on the heterocycle's imine carbon | `*C1=NOCC(=O)N1` |
| US20160002213A1 C00042 | p-0508/0509 defines R1 option (iii-1); TIF places the bar on the heteroaromatic carbon | `*c1noc(=O)[nH]1` |

The old graphs still provide a regression check for atoms and bonds outside R7.
`tests/r7_reference.py` first asserts the actual new graph and the exact source IDs
of the ports/deleted tails. On a copy only, it then restores those source carbon
nodes to compare the untouched remainder with the frozen historical graph, formula
and InChIKey. It does not accept either output indiscriminately, skip cases, update
the old answers or mark a new expected failure. `test_wavy_attachments.py` separately
asserts the real output including the ports. The historical pre-R7 suite also ran
unchanged on the preserved reader snapshot: 359 passed, one existing xfail.

The new `US20140206647A1_*.cdx` fixtures cover attachment bars, boronic acids,
nested nitro expansion and separately drawn final products. Source-selected
assembly is checked in HARVEST against those products and BindingDB's full keys
for US8981106 compounds 10 and 12. `US20100119475A1` fixtures are negative cases
(one tail and an internal wavy bond); `US20120035134A1_C00006` mixes attachment bars
with four other wavy bonds that must remain untouched.
