"""Classify source-audited differences using explicit join and label evidence.

This is a gate, not a chemical gold generator. Unknown flags remain unresolved;
MOL agreement alone never overrides an explicit source variable or graph.
"""
import argparse
import collections
import json
from pathlib import Path


def classify(record, trace, valence_trace=None):
    remaining = []
    evidence = []
    rules = {c['rule'] for c in record['changes']}
    isolated_nh2_repair = bool(valence_trace and valence_trace['calls'] and all(
        call['old_helper_returns_none'] and call['nitrogen_repairs'] and all(
            repair['explicit_h'] == [2, 1]
            and repair['source'].get('Element') == '7'
            and repair['source'].get('NumHydrogens') == '2'
            for repair in call['nitrogen_repairs']) for call in valence_trace['calls']))
    if isolated_nh2_repair:
        rules.add('isolated_nh2_repair_with_variables')
    source_rules = {'explicit_nested_graph', 'explicit_element_preserved',
                    'contradictory_nitrogen_hydrogens', 'declared_variable_placeholder',
                    'first_atom_abbreviation_expanded', 'variable_maps_detached',
                    'aromatic_variable_graph_retained'}
    joined_to_cut = bool(trace and any(c['joined'] and c['generated_attachments'] == 1
                                      and not c['source_attachments']
                                      for c in trace['traces']['after']['joins']))
    namespace_fix = bool(trace and not trace['after_consumed_source_variables'] and (
        trace['before_consumed_source_variables'] or trace['before_source_map_collision']))
    for flag in record['flags']:
        if flag == 'acceptance_loss' and record.get('source_placeholders_explaining_rejection'):
            evidence.append('source_placeholders_preclude_complete_molecule')
        elif flag == 'chemical_change_without_source_change' and isolated_nh2_repair:
            evidence.append('source_nh2_repair_verified_by_valence_replay')
        elif flag == 'prior_mol_agreement_lost' and rules & source_rules:
            if 'first_atom_abbreviation_expanded' in rules and not joined_to_cut:
                remaining.append(flag)
            else:
                evidence.append('explicit_source_graph_or_label_overrides_mol_agreement')
        elif flag == 'component_count_changed' and (namespace_fix or joined_to_cut):
            evidence.append('existing_selection_after_corrected_abbreviation_join')
        elif (flag == 'radicals_increased' and trace and (namespace_fix or joined_to_cut)
              and not trace['increased_radicals_on_retained_atoms']
              and trace['radicals_in_additional_source_components']
              and trace.get('additional_radicals_match_baseline_before_expansion')):
            evidence.append('radicals_belong_to_additional_source_components')
        else:
            remaining.append(flag)
    category = ('unresolved' if remaining else 'metadata_only'
                if record['classification'] == 'source_metadata_verified'
                else 'source_correction')
    return {'sha': record['sha'], 'member': record['member'], 'classification': category,
            'source_rules': sorted(rules), 'review_evidence': sorted(set(evidence)),
            'remaining_flags': remaining}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence', required=True)
    parser.add_argument('--traces', required=True)
    parser.add_argument('--delta-evidence', help='Successful audit_indexed_r_delta.py result')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    traces = {r['sha']: r for r in json.loads(Path(args.traces).read_text())}
    valence_traces = {}
    if args.delta_evidence:
        delta = json.loads(Path(args.delta_evidence).read_text())
        if delta['unresolved']:
            raise SystemExit('The supplied delta audit has unresolved changes.')
        valence_traces = {r['sha']: r for r in delta['recovered_graphs']}
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    counts, explanations, flags = collections.Counter(), collections.Counter(), collections.Counter()
    exceptional = []
    labels = r_labels = stereo = 0
    with (output / 'classifications.jsonl').open('w') as classifications:
        for line in Path(args.evidence).open():
            record = json.loads(line)
            result = classify(record, traces.get(record['sha']), valence_traces.get(record['sha']))
            counts[result['classification']] += 1
            explanations.update(result['review_evidence'])
            flags.update(result['remaining_flags'])
            labels += record.get('retained_source_labels_checked', 0)
            r_labels += record.get('retained_r_labels_checked', 0)
            stereo += record.get('retained_stereocenters_checked', 0)
            classifications.write(json.dumps(result) + '\n')
            if record['flags']:
                exceptional.append(result)
    summary = {'classifications': dict(counts), 'review_evidence': dict(explanations),
               'remaining_flags': dict(flags), 'source_labels_checked': labels,
               'r_labels_checked': r_labels, 'retained_stereocenters_checked': stereo}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    (output / 'reviewed_cases.json').write_text(json.dumps(exceptional, indent=2) + '\n')
    print(json.dumps(summary))
    if counts['unresolved']:
        raise SystemExit('Unresolved differences remain; inspect reviewed_cases.json.')


if __name__ == '__main__':
    main()
