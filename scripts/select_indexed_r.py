"""Select indexed-letter labels and previous failures for a focused reread.

The preceding run supplies selection only. Both readers must reread the selected
CDX bytes with compare_corpus.py; its old outputs are not comparison baselines.
Labels absent from the preceding selected molecule are outside this selection.
"""
import argparse
import json
from pathlib import Path
import re
import sqlite3

from compare_corpus import connect


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-db', required=True)
    parser.add_argument('--engine', required=True)
    parser.add_argument('--output-db', required=True)
    parser.add_argument('--references', default=str(
        Path(__file__).resolve().parents[1] / 'tests/data/patent_regressions/references.json'))
    args = parser.parse_args()
    if Path(args.output_db).exists():
        raise SystemExit('Choose a new output database; existing runs are preserved.')
    source = sqlite3.connect(f'file:{Path(args.source_db).resolve()}?mode=ro', uri=True)
    total = source.execute('SELECT COUNT(*) FROM files').fetchone()[0]
    count = source.execute('SELECT COUNT(*) FROM results WHERE engine=?', (args.engine,)).fetchone()[0]
    if total != count:
        raise SystemExit('Source selection run is incomplete.')
    references = json.loads(Path(args.references).read_text())
    controls = {r['sha256'] for r in references}
    matched, failed = set(), set()
    for sha, payload in source.execute('SELECT sha,payload FROM results WHERE engine=?', (args.engine,)):
        result = json.loads(payload)
        if result['status'] == 'none':
            failed.add(sha)
        if any(re.fullmatch(r'R\d+[a-z]', re.sub(r'\s+', '', str(label.get('_CDX_ATOM_LABEL', ''))))
               for label in result.get('labels', [])):
            matched.add(sha)
    selected = matched | failed | controls
    rows = [source.execute('SELECT * FROM files WHERE sha=?', (sha,)).fetchone() for sha in sorted(selected)]
    if any(row is None for row in rows):
        raise SystemExit('A portable control is missing from the source inventory.')
    output = connect(args.output_db)
    output.executemany('INSERT INTO files VALUES (?,?,?,?)', rows)
    for sha in sorted(selected):
        output.executemany('INSERT INTO origins VALUES (?,?,?,?)',
                           source.execute('SELECT * FROM origins WHERE sha=?', (sha,)))
    summary = {'scanned_results': total, 'matching_results': len(matched),
               'portable_controls': len(controls), 'previous_none_controls': len(failed),
               'selected_total': len(selected), 'source_engine': args.engine,
               'scope': 'Retained R + digits + lowercase letter labels, portable controls and previous None results; labels absent from the preceding selected result are not covered.'}
    output.execute('INSERT INTO metadata VALUES (?,?)', ('selection', json.dumps(summary)))
    output.commit()
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
