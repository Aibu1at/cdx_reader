"""The comparison harness must not silently skip failures or mix revisions."""
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import zipfile

import pytest

SPEC = importlib.util.spec_from_file_location(
    'compare_corpus', Path(__file__).parents[1] / 'scripts' / 'compare_corpus.py')
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)


def test_inventory_keeps_origins_without_merging_distinct_content(tmp_path):
    source = tmp_path / 'source'
    source.mkdir()
    for name, data in [('US1-version1', b'one'), ('US1-version2', b'two'), ('US2', b'one')]:
        with zipfile.ZipFile(source / (name + '.zip'), 'w') as archive:
            archive.writestr('C00001.CDX', data)
    db_path = tmp_path / 'corpus.sqlite'
    runner.inventory(SimpleNamespace(db=db_path, inputs=[source]))
    db = runner.connect(db_path)
    assert db.execute('SELECT COUNT(*) FROM files').fetchone()[0] == 2
    assert db.execute('SELECT COUNT(*) FROM origins').fetchone()[0] == 3
    assert {r[0] for r in db.execute('SELECT patent FROM origins')} == {'US1-version1', 'US1-version2', 'US2'}
    db.close()
    with zipfile.ZipFile(source / 'US1-version1.zip', 'w') as archive:
        archive.writestr('C00001.CDX', b'changed')
    with pytest.raises(SystemExit, match='Inventory incomplete'):
        runner.inventory(SimpleNamespace(db=db_path, inputs=[source]))


def test_resume_and_timeout_record_every_input(tmp_path):
    # Use the harness through its CLI so multiprocessing spawn imports __main__
    # normally, just as it does in a real corpus run.
    import subprocess
    import sys
    repo = tmp_path / 'reader'
    package = repo / 'cdx_reader'
    package.mkdir(parents=True)
    implementation = package / '__init__.py'
    implementation.write_text('''import time
def read_cdx_molecule(data, as_single=True):
    if data == b'slow':
        time.sleep(5)
    if data == b'bad':
        raise ValueError('known decoding failure')
    return None
''')
    archive_path = tmp_path / 'input.zip'
    with zipfile.ZipFile(archive_path, 'w') as archive:
        for name, data in [('A.CDX', b'good'), ('B.CDX', b'bad'), ('C.CDX', b'slow')]:
            archive.writestr(name, data)
    db_path = tmp_path / 'corpus.sqlite'
    runner.inventory(SimpleNamespace(db=db_path, inputs=[tmp_path]))
    command = [sys.executable, str(Path(runner.__file__)), 'run', '--db', str(db_path),
               '--repo', str(repo), '--engine', 'test', '--workers', '1', '--timeout', '1']
    subprocess.run(command, check=True, capture_output=True, text=True, timeout=30)
    db = runner.connect(db_path)
    rows = {sha: json.loads(payload) for sha, payload in db.execute('SELECT sha,payload FROM results')}
    assert rows[hashlib.sha256(b'good').hexdigest()]['status'] == 'none'
    assert rows[hashlib.sha256(b'bad').hexdigest()]['error'] == 'known decoding failure'
    assert rows[hashlib.sha256(b'slow').hexdigest()]['status'] == 'timeout'
    resume = subprocess.run(command, check=True, capture_output=True, text=True, timeout=30)
    assert 'pending 0' in resume.stdout
    assert db.execute('SELECT COUNT(*) FROM results').fetchone()[0] == 3
    implementation.write_text(implementation.read_text() + '\n# different source\n')
    rejected = subprocess.run(command, capture_output=True, text=True, timeout=30)
    assert rejected.returncode != 0
    assert 'Revision/environment changed' in rejected.stderr
    db.close()


@pytest.mark.parametrize('value', ['0', '-1'])
def test_run_budgets_must_be_positive(value):
    import argparse
    with pytest.raises(argparse.ArgumentTypeError):
        runner.positive_integer(value)
    with pytest.raises(argparse.ArgumentTypeError):
        runner.positive_number(value)


def test_stereo_audit_ignores_atom_reordering_but_detects_inversion():
    from rdkit import Chem
    spec = importlib.util.spec_from_file_location(
        'audit_corpus', Path(__file__).parents[1] / 'scripts' / 'audit_corpus.py')
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    mol = Chem.MolFromSmiles('N[C@@H](C)C(=O)O')
    for atom in mol.GetAtoms():
        atom.SetIntProp('_CDX_ATOM_ID', atom.GetIdx() + 100)
    shuffled = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    assert audit.tetrahedral_parities(mol) == audit.tetrahedral_parities(shuffled)
    next(a for a in shuffled.GetAtoms() if a.GetProp('_CDX_ATOM_ID') == '101').InvertChirality()
    assert audit.tetrahedral_parities(mol)['101'][1] == -audit.tetrahedral_parities(shuffled)['101'][1]


def test_missing_query_inchi_is_not_a_parser_failure(monkeypatch):
    from rdkit import Chem
    mol = Chem.MolFromSmiles('*c1ccccc1')

    def undefined_key(mol):
        raise ValueError('query has no unique InChI')

    monkeypatch.setattr(Chem, 'MolToInchiKey', undefined_key)
    result = runner.describe(mol)
    assert result['status'] == 'ok'
    assert result['smiles'] == '*c1ccccc1'
    assert result['formula'] == 'C6H5*'
    assert result['inchikey'] is None
    assert result['inchikey_error']['message'] == 'query has no unique InChI'
    assert not result['harvest_accepted']


@pytest.mark.parametrize('unexplained', ['graph_loss', 'tetrahedral_stereo_changed:7'])
@pytest.mark.parametrize('with_valence_replay', [False, True])
def test_review_gate_does_not_excuse_graph_or_stereo_loss_with_a_variable_label(unexplained, with_valence_replay):
    spec = importlib.util.spec_from_file_location(
        'review_audit', Path(__file__).parents[1] / 'scripts' / 'review_audit.py')
    review = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(review)
    record = {'sha': 'input', 'member': 'C00001.CDX',
              'classification': 'review',
              'changes': [{'rule': 'declared_variable_placeholder'}],
              'source_placeholders_explaining_rejection': [{'id': '8', 'label': 'X1'}],
              'flags': ['acceptance_loss', unexplained]}
    valence_trace = None
    if with_valence_replay:
        record['flags'].append('chemical_change_without_source_change')
        valence_trace = {'calls': [{'old_helper_returns_none': True, 'nitrogen_repairs': [
            {'explicit_h': [2, 1], 'source': {'Element': '7', 'NumHydrogens': '2'}}]}]}
    result = review.classify(record, None, valence_trace)
    assert result['classification'] == 'unresolved'
    assert result['remaining_flags'] == [unexplained]
