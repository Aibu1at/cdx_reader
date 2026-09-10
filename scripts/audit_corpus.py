"""Check changed CDX results against source atoms, attachments and styled labels.

This supplements compare_corpus.py. A passed local check establishes fidelity to
the explicit CDX graph, not the identity of a compound mentioned in patent prose.
Review flags must be assessed before describing a candidate as regression free.
"""
import argparse
import collections
import copy
import hashlib
import importlib.util
import io
import itertools
import multiprocessing as mp
import re
import json
import logging
from pathlib import Path
import sqlite3
import sys
from zipfile import ZipFile


def load_reader(name, repo):
    package = Path(repo) / 'cdx_reader'
    spec = importlib.util.spec_from_file_location(
        name, package / '__init__.py', submodule_search_locations=[str(package)])
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return sys.modules[name + '.cdx_reader']


def keyed_atoms(result):
    return {str(a['cdx_id']): a for a in result.get('atoms', []) if a['cdx_id'] is not None}


def keyed_bonds(result):
    atoms = {a['index']: str(a['cdx_id']) for a in result.get('atoms', [])}
    return {tuple(sorted((atoms[a], atoms[b]))): (order, stereo)
            for a, b, order, stereo in result.get('bonds', [])
            if atoms[a] != 'None' and atoms[b] != 'None'}


def tetrahedral_parities(mol):
    """Compare tetrahedral orientation independently of RDKit atom ordering."""
    result = {}
    if mol is None:
        return result
    for atom in mol.GetAtoms():
        if not atom.HasProp('_CDX_ATOM_ID'):
            continue
        neighbors = list(atom.GetNeighbors())
        if any(not neighbor.HasProp('_CDX_ATOM_ID') for neighbor in neighbors):
            continue
        ids = [neighbor.GetProp('_CDX_ATOM_ID') for neighbor in neighbors]
        inversions = sum(ids[i] > ids[j] for i in range(len(ids)) for j in range(i + 1, len(ids)))
        tag = str(atom.GetChiralTag())
        parity = {'CHI_TETRAHEDRAL_CW': 1, 'CHI_TETRAHEDRAL_CCW': -1}.get(tag, 0)
        result[atom.GetProp('_CDX_ATOM_ID')] = (tuple(sorted(ids)), parity * (-1) ** inversions)
    return result


def _initialize(args):
    global _old, _new, _normalization, _Chem, _converter, _archives
    _old = load_reader('audit_before', args.before_repo)
    _new = load_reader('audit_after', args.after_repo)
    _normalization = sys.modules['audit_after._normalization']
    from pycdxml import cdxml_converter
    from rdkit import Chem, RDLogger
    _Chem, _converter = Chem, cdxml_converter
    _archives = collections.OrderedDict()
    logging.disable(logging.CRITICAL)
    RDLogger.DisableLog('rdApp.*')


def _audit_case(case):
    sha, path, member, aa, bb = case
    old, new, normalization = _old, _new, _normalization
    Chem, cdxml_converter = _Chem, _converter
    if path not in _archives:
        if len(_archives) >= 8:
            _archives.popitem(last=False)[1].close()
        _archives[path] = ZipFile(path)
    _archives.move_to_end(path)
    zf = _archives[path]
    a, b = json.loads(aa), json.loads(bb)
    record = {'sha': sha, 'zip': path, 'member': member, 'flags': [], 'changes': []}
    flags = record['flags']
    atoms, bonds = keyed_atoms(b), keyed_bonds(b)
    chemical = a.get('smiles') != b.get('smiles') or a['status'] != b['status']
    try:
        data = zf.read(member)
        if hashlib.sha256(data).hexdigest() != sha:
            raise ValueError('Input changed since comparison')
        doc = cdxml_converter.read_cdx(io.BytesIO(data))
        original = {n.get('id'): n for n in doc.cdxml.getroot().iter('n')}
        # Check every retained text label, including non-R variables and labels
        # whose interpretation is deferred to a future agent.
        expected_all = normalization.capture_atom_labels(doc.cdxml)
        expected_r = normalization.capture_r_labels(doc.cdxml)
        actual = {str(item['cdx_id']): item for item in b.get('labels', [])
                  if '_CDX_ATOM_LABEL' in item}
        for identifier, atom in atoms.items():
            if identifier not in expected_all:
                continue
            raw, node_type = expected_all[identifier]
            label = actual.get(identifier, {})
            if str(label.get('_CDX_ATOM_LABEL')) != raw:
                # Older comparison files used RDKit's automatic conversion of
                # string properties (e.g. "3." -> 3.0). Check the actual atom
                # when that lossy serialization, rather than parsing, differs.
                reparsed = new.read_cdx_molecule(data)
                recovered = {atom.GetProp('_CDX_ATOM_ID'): atom.GetProp('_CDX_ATOM_LABEL')
                             for atom in reparsed.GetAtoms() if atom.HasProp('_CDX_ATOM_ID')
                             and atom.HasProp('_CDX_ATOM_LABEL')} if reparsed is not None else {}
                if recovered.get(identifier) != raw:
                    flags.append('source_label_mismatch:' + identifier)
            if atom['element'] == 0 and label.get('dummyLabel') != re.sub(r'\s+', '', raw):
                flags.append('dummy_label_mismatch:' + identifier)
            if atom['element'] == 0 and identifier in expected_r:
                normalized = expected_r[identifier][1]
                if (normalized[1:].isdigit() and int(normalized[1:]) <= 2147483647
                        and label.get('_MolFileRLabel') != int(normalized[1:])):
                    flags.append('r_number_mismatch:' + identifier)
        if any(identifier not in expected_all for identifier in actual):
            flags.append('label_without_source')
        record['retained_source_labels_checked'] = len(actual)
        record['retained_r_labels_checked'] = sum(identifier in expected_r and atoms[identifier]['element'] == 0
                                                 for identifier in actual if identifier in atoms)
        if chemical:
            before_mol = old.read_cdx_molecule(data)
            after_mol = new.read_cdx_molecule(data)
            old_stereo, new_stereo = tetrahedral_parities(before_mol), tetrahedral_parities(after_mol)
            checked = 0
            for identifier in old_stereo.keys() & new_stereo.keys():
                neighbors, parity = old_stereo[identifier]
                if parity and new_stereo[identifier][0] == neighbors:
                    checked += 1
                    if new_stereo[identifier][1] != parity:
                        flags.append('tetrahedral_stereo_changed:' + identifier)
            record['retained_stereocenters_checked'] = checked
            left, right = copy.deepcopy(doc.cdxml), copy.deepcopy(doc.cdxml)
            old.parse_cdxml(left)
            new.parse_cdxml(right)
            old_nodes = {n.get('id'): n for n in left.getroot().iter('n')}
            for node in right.getroot().iter('n'):
                identifier = node.get('id')
                previous, source = old_nodes[identifier], original[identifier]
                delta = {key: [previous.get(key), node.get(key)]
                         for key in sorted(set(previous.attrib) | set(node.attrib))
                         if previous.get(key) != node.get(key)}
                if not delta:
                    continue
                change = {'id': identifier, 'label': new.atom_label(node), 'attributes': delta}
                if (node.get('NodeType') == 'Fragment'
                        and normalization.has_unambiguous_nested_fragment(source)):
                    change['rule'] = 'explicit_nested_graph'
                    fragment = source.find('fragment')
                    ext = next(n for n in fragment.findall('n')
                               if n.get('NodeType') == 'ExternalConnectionPoint')
                    children = {n.get('id'): n for n in fragment.findall('n') if n is not ext}
                    retained = set(children) & set(atoms)
                    if retained and retained != set(children):
                        flags.append('partially_retained_nested_graph:' + identifier)
                    for child_id in retained:
                        child, atom = children[child_id], atoms[child_id]
                        if (atom['element'] != int(child.get('Element', '6'))
                                or atom['charge'] != int(child.get('Charge', '0'))
                                or atom['isotope'] != int(child.get('Isotope', '0'))):
                            if keyed_atoms(a).get(child_id) != atom:
                                # Source Fragment nodes can already contain a
                                # deferred variable encoded as a metal; a
                                # retained legacy interpretation is not a new
                                # nested-graph corruption.
                                previous_atom = keyed_atoms(a).get(child_id, {})
                                if any(previous_atom.get(k) != atom[k]
                                       for k in ('element', 'charge', 'isotope')):
                                    flags.append('nested_atom_mismatch:' + child_id)
                            change.setdefault('legacy_source_atom_mismatches', []).append(child_id)
                    for bond in fragment.findall('b'):
                        begin, end = bond.get('B'), bond.get('E')
                        if ext.get('id') in (begin, end):
                            internal = end if begin == ext.get('id') else begin
                            outer = next(bd for bd in source.getparent().findall('b')
                                         if identifier in (bd.get('B'), bd.get('E')))
                            neighbor = outer.get('E') if outer.get('B') == identifier else outer.get('B')
                            begin, end = internal, neighbor
                            order = outer.get('Order', '1')
                        else:
                            order = bond.get('Order', '1')
                        if begin in atoms and end in atoms:
                            got = bonds.get(tuple(sorted((begin, end))))
                            expected = {'1': 'SINGLE', '2': 'DOUBLE', '3': 'TRIPLE'}[order]
                            if not got or got[0] not in (expected, 'AROMATIC'):
                                flags.append('nested_bond_mismatch:' + identifier)
                    change['retained_atoms_checked'] = len(retained)
                    if not retained:
                        change['rule'] = ('nested_repair_deferred' if identifier in atoms else
                                          'nested_graph_not_in_selected_result')
                elif 'NumHydrogens' in delta:
                    change['rule'] = 'contradictory_nitrogen_hydrogens'
                    if identifier in atoms and atoms[identifier]['radicals']:
                        flags.append('nitrogen_still_radical:' + identifier)
                elif node.get('Element') == '0' and normalization.variable_label(source, new.FRAGMENT_LIBRARY):
                    change['rule'] = 'declared_variable_placeholder'
                    change['source_node_type'] = source.get('NodeType')
                    if identifier in atoms and atoms[identifier]['element'] != 0:
                        flags.append('variable_not_retained_as_dummy:' + identifier)
                elif 'Element' in delta and node.get('Element') == source.get('Element'):
                    change['rule'] = 'explicit_element_preserved'
                    if identifier in atoms:
                        atom = atoms[identifier]
                        if (atom['element'] != int(source.get('Element'))
                                or atom['charge'] != int(source.get('Charge', '0'))
                                or atom['isotope'] != int(source.get('Isotope', '0'))):
                            flags.append('explicit_element_mismatch:' + identifier)
                elif normalization.r_label(source):
                    change['rule'] = 'r_label'
                else:
                    change['rule'] = 'unclassified'
                    flags.append('unclassified_source_change:' + identifier)
                record['changes'].append(change)
            if not record['changes']:
                reparsed = new.read_cdx_molecule(data)
                maps = {a.GetProp('_CDX_ATOM_ID'): a.GetIntProp('_CDX_ORIGINAL_ATOM_MAP')
                        for a in reparsed.GetAtoms() if a.HasProp('_CDX_ORIGINAL_ATOM_MAP')
                        and a.HasProp('_CDX_ATOM_ID')} if reparsed is not None else {}
                if maps:
                    record['changes'].append({'rule': 'variable_maps_detached', 'original_maps': maps})
                elif after_mol is not None and after_mol.HasProp('_CDX_UNKEKULIZED_QUERY'):
                    record['changes'].append({'rule': 'aromatic_variable_graph_retained',
                                              'inchikey_error': b.get('inchikey_error')})
                else:
                    # A repaired boundary in the existing expander permits a
                    # dictionary abbreviation at RDKit atom index zero. Replay
                    # that one helper with the baseline implementation to
                    # isolate its effect, rather than guessing from a SMILES.
                    split = new.split_mol_at_atom
                    new.split_mol_at_atom = old.split_mol_at_atom
                    try:
                        replay = new.read_cdx_molecule(data)
                    finally:
                        new.split_mol_at_atom = split
                    if replay is not None:
                        for atom in replay.GetAtoms():
                            atom.SetAtomMapNum(0)
                    replay_smiles = Chem.MolToSmiles(replay) if replay is not None else None
                    if replay_smiles == a.get('smiles') and a.get('smiles') != b.get('smiles'):
                        record['changes'].append({'rule': 'first_atom_abbreviation_expanded',
                                                  'baseline_split_replay_matches': True})
                    else:
                        flags.append('chemical_change_without_source_change')
            try:
                mol = Chem.MolFromMolBlock(zf.read(member[:-4] + '.MOL').decode(errors='replace'),
                                           sanitize=True, removeHs=False)
                key = Chem.MolToInchiKey(mol) if mol is not None else None
                record['mol_inchikey'] = key
                record['before_matches_mol'] = bool(key) and a.get('inchikey') == key
                record['after_matches_mol'] = bool(key) and b.get('inchikey') == key
                if record['before_matches_mol'] and not record['after_matches_mol']:
                    flags.append('prior_mol_agreement_lost')
            except Exception as exc:
                record['mol_error'] = str(exc)
            if a.get('harvest_accepted') and not b.get('harvest_accepted'):
                flags.append('acceptance_loss')
                placeholders = []
                for identifier, atom in atoms.items():
                    if atom['element'] != 0:
                        continue
                    node = original.get(identifier)
                    if node is None:
                        continue
                    label = normalization.atom_label(node)
                    declared = (normalization.variable_label(node, new.FRAGMENT_LIBRARY)
                                or label in new.SCAFFOLD_LABELS or new.get_fragment_smiles(label)[2])
                    if declared:
                        placeholders.append({'id': identifier, 'label': label,
                                             'source_type': node.get('NodeType'),
                                             'source_degree': len(normalization._bonds_to(node.getparent(), identifier))})
                record['source_placeholders_explaining_rejection'] = placeholders
            if a['status'] == 'ok' and b['status'] != 'ok':
                flags.append('graph_loss')
            if a['status'] == b['status'] == 'ok':
                if len(a['components']) != len(b['components']):
                    flags.append('component_count_changed')
                if b['radicals'] > a['radicals']:
                    flags.append('radicals_increased')
                old_bonds = keyed_bonds(a)
                for pair in old_bonds.keys() & bonds.keys():
                    if old_bonds[pair][1] != 'STEREONONE' and old_bonds[pair][1] != bonds[pair][1]:
                        flags.append('bond_stereo_changed:' + ','.join(pair))
        record['classification'] = ('review' if flags else
                                    'source_graph_verified' if chemical else 'source_metadata_verified')
    except Exception as exc:
        flags.append('audit_error:' + str(exc))
        record['classification'] = 'review'
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--db', required=True)
    parser.add_argument('--before-db', required=True)
    parser.add_argument('--before', default='baseline')
    parser.add_argument('--after', default='candidate')
    parser.add_argument('--before-repo', required=True)
    parser.add_argument('--after-repo', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--workers', type=int, default=6)
    parser.add_argument('--limit', type=int)
    parser.add_argument('--chemistry-only', action='store_true', help='Skip unchanged graphs with only label metadata changes')
    args = parser.parse_args()
    # Pool's feeder thread consumes this read-only cursor; no other thread uses
    # the connection while results are streamed to workers.
    db = sqlite3.connect(args.db, check_same_thread=False)
    db.execute('ATTACH DATABASE ? AS previous', (args.before_db,))
    from compare_corpus import revision_info
    for table, engine, repo in [('previous.metadata', args.before, args.before_repo),
                                ('metadata', args.after, args.after_repo)]:
        metadata = db.execute(f'SELECT value FROM {table} WHERE key=?', ('run:' + engine,)).fetchone()
        if (not metadata or json.loads(metadata[0])['reader_source_sha256']
                != revision_info(repo)['reader_source_sha256']):
            raise SystemExit('Audit source does not match the recorded reader revision: ' + engine)
    cases = db.execute('''SELECT a.sha,f.zip,f.member,a.payload,b.payload
        FROM previous.results a JOIN results b ON a.sha=b.sha JOIN files f ON f.sha=a.sha
        WHERE a.engine=? AND b.engine=? AND (
        json_extract(a.payload,'$.smiles') IS NOT json_extract(b.payload,'$.smiles') OR
        json_extract(a.payload,'$.status') IS NOT json_extract(b.payload,'$.status') OR
        (?=0 AND json_extract(a.payload,'$.labels') IS NOT json_extract(b.payload,'$.labels')))
        ORDER BY f.zip,f.member''', (args.before, args.after, int(args.chemistry_only)))
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    counts, reasons = collections.Counter(), collections.Counter()
    if args.workers <= 0 or args.limit is not None and args.limit <= 0:
        raise SystemExit('workers and limit must be positive')
    if args.limit:
        cases = itertools.islice(cases, args.limit)
    with (output / 'evidence.jsonl').open('w') as evidence, (output / 'review.jsonl').open('w') as review:
        with mp.get_context('spawn').Pool(args.workers, initializer=_initialize, initargs=(args,)) as pool:
            for i, record in enumerate(pool.imap(_audit_case, cases, chunksize=16)):
                counts[record['classification']] += 1
                for change in record['changes']:
                    reasons[change['rule']] += 1
                evidence.write(json.dumps(record) + '\n')
                if record['flags']:
                    review.write(json.dumps(record) + '\n')
                if (i + 1) % 1000 == 0:
                    print('audit', i + 1, flush=True)
    summary = {'changed_inputs_checked': sum(counts.values()), 'classifications': dict(counts), 'node_rules': dict(reasons)}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary), flush=True)


if __name__ == '__main__':
    main()
