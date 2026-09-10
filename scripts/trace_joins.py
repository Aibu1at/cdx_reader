"""Trace abbreviation joins implicated by a completed source audit.

Source variable nodes must never be consumed as temporary cut atoms. This tool
reopens flagged CDX inputs and instruments both readers without modifying them.
"""
import argparse
import collections
import json
import logging
from pathlib import Path
from zipfile import ZipFile

from audit_corpus import load_reader


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', required=True)
    parser.add_argument('--before-repo', required=True)
    parser.add_argument('--after-repo', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    old = load_reader('join_before', args.before_repo)
    new = load_reader('join_after', args.after_repo)
    from rdkit import Chem, RDLogger
    from pycdxml import cdxml_converter
    from io import BytesIO
    import hashlib
    logging.disable(logging.CRITICAL)
    RDLogger.DisableLog('rdApp.*')
    records = []
    for line in Path(args.evidence).open():
        source = json.loads(line)
        if ('component_count_changed' not in source['flags'] and not any(
                c['rule'] in ('variable_maps_detached', 'first_atom_abbreviation_expanded')
                for c in source['changes'])):
            continue
        with ZipFile(source['zip']) as archive:
            data = archive.read(source['member'])
        if hashlib.sha256(data).hexdigest() != source['sha']:
            raise ValueError('Input changed since source audit')
        document = cdxml_converter.read_cdx(BytesIO(data))
        labels = new.capture_atom_labels(document.cdxml)
        row = {'sha': source['sha'], 'member': source['member'],
               'flags': source['flags'], 'traces': {}}
        retained = {}
        for reader, engine in [(old, 'before'), (new, 'after')]:
            original = reader.join_molecule_fragments
            calls = []

            def traced(scaffold, decoration, **kwargs):
                maps = {a.GetAtomMapNum() for a in decoration.GetAtoms()
                        if a.GetAtomicNum() == 0 and a.GetAtomMapNum()}
                used = [a for a in scaffold.GetAtoms()
                        if a.GetAtomicNum() == 0 and a.GetAtomMapNum() in maps]
                call = {'source_attachments': [
                    {'id': a.GetProp('_CDX_ATOM_ID'),
                     'text': labels.get(a.GetProp('_CDX_ATOM_ID')),
                     'isotope': a.GetIsotope()}
                    for a in used if a.HasProp('_CDX_ATOM_ID')],
                    'generated_attachments': sum(not a.HasProp('_CDX_ATOM_ID') for a in used),
                    'decoration_smiles': Chem.MolToSmiles(decoration),
                    'source_neighbors': [neighbor.GetProp('_CDX_ATOM_ID')
                                         for atom in used for neighbor in atom.GetNeighbors()
                                         if neighbor.HasProp('_CDX_ATOM_ID')]}
                result = original(scaffold, decoration, **kwargs)
                call['joined'] = result is not None
                calls.append(call)
                return result

            reader.join_molecule_fragments = traced
            try:
                mol = reader.read_cdx_molecule(data)
            finally:
                reader.join_molecule_fragments = original
            row['traces'][engine] = {
                'joins': calls,
                'components': len(Chem.GetMolFrags(mol)) if mol is not None else None,
                'radicals': {a.GetProp('_CDX_ATOM_ID'): [a.GetAtomicNum(), a.GetNumRadicalElectrons()]
                             for a in mol.GetAtoms() if a.GetNumRadicalElectrons()
                             and a.HasProp('_CDX_ATOM_ID')} if mol is not None else {}}
            retained[engine] = {a.GetProp('_CDX_ATOM_ID') for a in mol.GetAtoms()
                                if a.HasProp('_CDX_ATOM_ID')} if mol is not None else set()
        for engine in ('before', 'after'):
            row[engine + '_consumed_source_variables'] = any(
                c['joined'] and c['source_attachments'] for c in row['traces'][engine]['joins'])
            row[engine + '_source_map_collision'] = any(
                not c['joined'] and c['source_attachments'] and c['generated_attachments']
                for c in row['traces'][engine]['joins'])
        before_radicals = row['traces']['before']['radicals']
        after_radicals = row['traces']['after']['radicals']
        row['increased_radicals_on_retained_atoms'] = {
            identifier: value for identifier, value in after_radicals.items()
            if identifier in retained['before'] and value[1] > before_radicals.get(identifier, [0, 0])[1]}
        row['radicals_in_additional_source_components'] = sorted(set(after_radicals) - retained['before'])
        if 'radicals_increased' in source['flags']:
            # Isolate selection from chemistry: the baseline already assigns
            # these radicals before its erroneous abbreviation join discards
            # their components. Merely finding new source IDs is insufficient.
            expand = old.expand_fragments_in_mol
            old.expand_fragments_in_mol = lambda mol, fragments: mol
            try:
                unexpanded = old.read_cdx_molecule(data)
            finally:
                old.expand_fragments_in_mol = expand
            baseline_radicals = {a.GetProp('_CDX_ATOM_ID'): a.GetNumRadicalElectrons()
                                 for a in unexpanded.GetAtoms() if a.HasProp('_CDX_ATOM_ID')
                                 and a.GetNumRadicalElectrons()} if unexpanded is not None else {}
            row['baseline_before_expansion_radicals'] = baseline_radicals
            row['additional_radicals_match_baseline_before_expansion'] = all(
                baseline_radicals.get(identifier) == after_radicals[identifier][1]
                for identifier in row['radicals_in_additional_source_components'])
        records.append(row)
    Path(args.output).write_text(json.dumps(records, indent=2) + '\n')
    counts = collections.Counter(key for row in records for key in (
        'before_consumed_source_variables', 'after_consumed_source_variables',
        'before_source_map_collision', 'after_source_map_collision') if row[key])
    print(json.dumps({'records': len(records), **counts}))


if __name__ == '__main__':
    main()
