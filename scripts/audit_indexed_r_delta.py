"""Strict follow-up gate: only indexed R substitutions or isolated NH2 repairs.

Use together with audit_corpus.py, which rereads source labels and stereo. This
checks every stored atom/bond, including atoms generated during abbreviation
expansion, and traces the few recovered graphs through both valence helpers.
"""
import argparse
from collections import Counter
import hashlib
from io import BytesIO
import json
from pathlib import Path
import re
import sqlite3
from zipfile import ZipFile

from audit_corpus import load_reader, tetrahedral_parities
from compare_corpus import revision_info


def trace_recovery(data, old, new):
    from rdkit import Chem
    from pycdxml import cdxml_converter
    document = cdxml_converter.read_cdx(BytesIO(data))
    nodes = {n.get('id'): n for n in document.cdxml.getroot().iter('n')}
    calls = []
    original = new.fix_valence

    def traced(mol):
        before = Chem.Mol(mol)
        result = original(mol)
        if result is None:
            return None
        changes = []
        assert old.fix_valence(Chem.Mol(before)) is None
        assert before.GetNumAtoms() == result.GetNumAtoms()
        assert before.GetNumBonds() == result.GetNumBonds()
        assert tetrahedral_parities(before) == tetrahedral_parities(result)
        for left, right in zip(before.GetBonds(), result.GetBonds()):
            assert (left.GetBeginAtomIdx(), left.GetEndAtomIdx(), left.GetBondDir(), left.GetStereo()) == (
                right.GetBeginAtomIdx(), right.GetEndAtomIdx(), right.GetBondDir(), right.GetStereo())
            if left.GetBondType() != right.GetBondType():
                # Sanitization may aromatize alternating single/double bonds.
                assert left.IsInRing() and right.IsInRing()
                assert left.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.AROMATIC)
                assert right.GetBondType() in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE, Chem.BondType.AROMATIC)
        for left, right in zip(before.GetAtoms(), result.GetAtoms()):
            assert (left.GetAtomicNum(), left.GetFormalCharge(), left.GetIsotope(), left.GetNoImplicit()) == (
                right.GetAtomicNum(), right.GetFormalCharge(), right.GetIsotope(), right.GetNoImplicit())
            assert not right.GetNumRadicalElectrons()
            if left.GetNumExplicitHs() == right.GetNumExplicitHs():
                continue
            identifier = left.GetProp('_CDX_ATOM_ID')
            source = nodes[identifier]
            assert left.GetAtomicNum() == 7 and not left.GetFormalCharge()
            assert not left.HasQuery() and not left.GetNumRadicalElectrons()
            assert left.GetNumExplicitHs() == 2 and right.GetNumExplicitHs() == 1
            assert left.GetDegree() == 2
            assert all(b.GetBondType() == Chem.BondType.SINGLE for b in left.GetBonds())
            assert source.get('Element') == '7' and source.get('NumHydrogens') == '2'
            assert source.get('Charge', '0') == '0' and not source.get('Radical')
            changes.append({'id': identifier, 'source': dict(source.attrib), 'explicit_h': [2, 1]})
        assert changes and any(a.GetAtomicNum() == 0 for a in before.GetAtoms())
        calls.append({'old_helper_returns_none': True, 'nitrogen_repairs': changes})
        return result

    new.fix_valence = traced
    try:
        mol = new.read_cdx_molecule(data)
    finally:
        new.fix_valence = original
    assert mol is not None and calls
    return calls


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', required=True)
    parser.add_argument('--before', required=True)
    parser.add_argument('--after', required=True)
    parser.add_argument('--before-repo', required=True)
    parser.add_argument('--after-repo', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    db = sqlite3.connect(args.db)
    for engine, repo in [(args.before, args.before_repo), (args.after, args.after_repo)]:
        recorded = json.loads(db.execute('SELECT value FROM metadata WHERE key=?', ('run:' + engine,)).fetchone()[0])
        assert recorded['reader_source_sha256'] == revision_info(repo)['reader_source_sha256']
    old, new = load_reader('delta_before', args.before_repo), load_reader('delta_after', args.after_repo)
    counts, recoveries, unresolved = Counter(), [], []
    for sha, path, member, aa, bb in db.execute('''
        SELECT a.sha,f.zip,f.member,a.payload,b.payload FROM results a JOIN results b USING(sha)
        JOIN files f USING(sha) WHERE a.engine=? AND b.engine=?''', (args.before, args.after)):
        a, b = json.loads(aa), json.loads(bb)
        counts['compared'] += 1
        try:
            if a['status'] == 'none' and b['status'] == 'ok':
                assert not b['harvest_accepted'] and not b['radicals']
                with ZipFile(path) as archive:
                    data = archive.read(member)
                assert hashlib.sha256(data).hexdigest() == sha
                calls = trace_recovery(data, old, new)
                recoveries.append({'sha': sha, 'member': member, 'calls': calls})
                counts['isolated_nh2_recoveries'] += 1
                continue
            assert a['status'] == b['status'], 'parser status lost'
            if a['status'] != 'ok':
                assert a == b, 'failure payload changed'
                counts['unchanged_failures'] += 1
                continue
            assert a['bonds'] == b['bonds'], 'bond topology, order or stereo changed'
            assert len(a['atoms']) == len(b['atoms']), 'atom set changed'
            assert len(a['components']) == len(b['components']), 'component count changed'
            before_labels = {x['index']: x for x in a['labels']}
            after_labels = {x['index']: x for x in b['labels']}
            changed = 0
            for left, right in zip(a['atoms'], b['atoms']):
                if left == right:
                    continue
                old_label = before_labels[left['index']]
                new_label = after_labels[right['index']]
                label = re.sub(r'\s+', '', old_label['_CDX_ATOM_LABEL'])
                assert re.fullmatch(r'R\d+[a-z]', label), 'unrelated atom changed'
                assert left['element'] == 6 and right['element'] == 0
                assert {**left, 'element': 0} == right, 'other atom properties changed'
                assert new_label['_CDX_ATOM_LABEL'] == old_label['_CDX_ATOM_LABEL']
                assert new_label['dummyLabel'] == label and '_MolFileRLabel' not in new_label
                changed += 1
            counts['indexed_variable_atoms'] += changed
            counts['indexed_variable_inputs' if changed else 'unchanged_graphs'] += 1
        except Exception as exc:
            unresolved.append({'sha': sha, 'member': member, 'error': repr(exc)})
    assert counts['compared'] == db.execute('SELECT COUNT(*) FROM files').fetchone()[0]
    output = {'counts': dict(counts), 'recovered_graphs': recoveries, 'unresolved': unresolved}
    Path(args.output).write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps({'counts': dict(counts), 'unresolved': len(unresolved)}))
    if unresolved:
        raise SystemExit('Unexpected changes remain; inspect the output.')


if __name__ == '__main__':
    main()
