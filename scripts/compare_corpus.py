"""Reproducible, resumable before/after CDX comparison (no network calls).

Run with the same Python/dependencies for both revisions. --repo selects the
reader checkout inside isolated worker processes; it does not install it.
The SQLite inventory stores content hashes and every original ZIP/member pair.
"""
import argparse
import collections
import hashlib
import json
import math
import multiprocessing as mp
from multiprocessing.connection import wait
from pathlib import Path
import sqlite3
import shutil
import sys
import time
import zipfile


def connect(path):
    db = sqlite3.connect(str(path), timeout=60)
    db.execute('PRAGMA journal_mode=WAL')
    db.execute('PRAGMA synchronous=NORMAL')
    db.executescript('''
        CREATE TABLE IF NOT EXISTS files (
            sha TEXT PRIMARY KEY, zip TEXT, member TEXT, size INTEGER);
        CREATE TABLE IF NOT EXISTS origins (
            sha TEXT, zip TEXT, member TEXT, patent TEXT,
            PRIMARY KEY(zip, member));
        CREATE TABLE IF NOT EXISTS results (
            engine TEXT, sha TEXT, payload TEXT, seconds REAL,
            PRIMARY KEY(engine, sha));
        CREATE TABLE IF NOT EXISTS metadata (key TEXT PRIMARY KEY, value TEXT);
    ''')
    return db


def inventory(args):
    db = connect(args.db)
    archives = sorted({str(z.resolve()) for directory in args.inputs
                       for z in Path(directory).glob('*.zip')})
    failures = []
    for i, path in enumerate(archives):
        try:
            with zipfile.ZipFile(path) as z:
                for member in sorted(z.namelist()):
                    if not member.lower().endswith('.cdx'):
                        continue
                    data = z.read(member)
                    sha = hashlib.sha256(data).hexdigest()
                    previous = db.execute('SELECT sha FROM origins WHERE zip=? AND member=?',
                                          (path, member)).fetchone()
                    if previous and previous[0] != sha:
                        raise ValueError('Input changed since inventory: ' + member)
                    db.execute('INSERT OR IGNORE INTO files VALUES (?,?,?,?)',
                               (sha, path, member, len(data)))
                    db.execute('INSERT OR IGNORE INTO origins VALUES (?,?,?,?)',
                               (sha, path, member, Path(path).stem))
            db.commit()
        except Exception as exc:
            db.rollback()
            failures.append({'zip': path, 'error': str(exc)})
        if (i + 1) % 50 == 0:
            print('inventory', i + 1, '/', len(archives), flush=True)
    info = {'archives': len(archives), 'failures': failures,
            'unique_files': db.execute('SELECT COUNT(*) FROM files').fetchone()[0],
            'origins': db.execute('SELECT COUNT(*) FROM origins').fetchone()[0],
            'patents_with_cdx': db.execute('SELECT COUNT(DISTINCT patent) FROM origins').fetchone()[0]}
    db.execute('INSERT OR REPLACE INTO metadata VALUES (?,?)', ('inventory', json.dumps(info)))
    db.commit()
    print(json.dumps(info), flush=True)
    if failures:
        raise SystemExit('Inventory incomplete; inspect failures before running.')


def describe(mol):
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors
    if mol is None:
        return {'status': 'none'}
    atoms = []
    labels = []
    for atom in mol.GetAtoms():
        props = atom.GetPropsAsDict(includePrivate=True, includeComputed=False,
                                   autoConvertStrings=False)
        entry = {'index': atom.GetIdx(), 'cdx_id': props.get('_CDX_ATOM_ID'),
                 'element': atom.GetAtomicNum(), 'charge': atom.GetFormalCharge(),
                 'isotope': atom.GetIsotope(), 'radicals': atom.GetNumRadicalElectrons(),
                 'explicit_h': atom.GetNumExplicitHs(), 'no_implicit': atom.GetNoImplicit(),
                 'query': atom.DescribeQuery() if atom.HasQuery() else ''}
        atoms.append(entry)
        selected = {k: props[k] for k in ('_CDX_ATOM_LABEL', '_MolFileRLabel',
                                          'dummyLabel', 'atomLabel', 'molFileAlias') if k in props}
        if selected:
            labels.append({'cdx_id': entry['cdx_id'], 'index': atom.GetIdx(), **selected})
    clean = Chem.Mol(mol)
    # Atom maps are annotations, not isotopes or stereochemistry.
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
    components = sorted(Chem.MolToSmiles(f, isomericSmiles=True)
                        for f in Chem.GetMolFrags(clean, asMols=True, sanitizeFrags=False))
    # Keep the baseline measurement order. SMILES initializes RDKit's stereo
    # caches; invoking InChI first on large patent drawings can be much slower.
    smiles = Chem.MolToSmiles(clean, isomericSmiles=True)
    graph_smiles = Chem.MolToSmiles(clean, isomericSmiles=False)
    # An aromatic query with unspecified variable electron counts can have a
    # valid retained graph and no unique InChI. Do not report a descriptor's
    # exception as an exception from the CDX parser.
    descriptor_error = {}
    try:
        inchikey = Chem.MolToInchiKey(Chem.Mol(clean))
    except Exception as exc:
        inchikey = None
        descriptor_error['inchikey_error'] = {'type': type(exc).__name__, 'message': str(exc)}
    return {'status': 'ok', 'smiles': smiles,
            'graph_smiles': graph_smiles,
            'inchikey': inchikey, **descriptor_error,
            'formula': rdMolDescriptors.CalcMolFormula(clean),
            'charge': Chem.GetFormalCharge(clean), 'components': components,
            'radicals': sum(a['radicals'] for a in atoms),
            'atoms': atoms, 'bonds': [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(),
                                      str(b.GetBondType()), str(b.GetStereo())) for b in mol.GetBonds()],
            'labels': labels,
            # Exact current HARVEST parse_cdx_bytes acceptance after parsing:
            # None/error and atomic numbers 0,99,100 are rejected.
            'harvest_accepted': not any(a['element'] in (0, 99, 100) for a in atoms)}


def worker(pipe, repo):
    sys.path.insert(0, repo)
    from cdx_reader import read_cdx_molecule
    from rdkit import RDLogger
    import logging
    logging.disable(logging.CRITICAL)
    RDLogger.DisableLog('rdApp.*')
    cache = collections.OrderedDict()
    pipe.send(('ready', None))
    while True:
        task = pipe.recv()
        if task is None:
            break
        sha, archive, member = task
        start = time.monotonic()
        try:
            if archive not in cache:
                if len(cache) >= 8:
                    cache.popitem(last=False)[1].close()
                cache[archive] = zipfile.ZipFile(archive)
            cache.move_to_end(archive)
            data = cache[archive].read(member)
            if hashlib.sha256(data).hexdigest() != sha:
                raise ValueError('CDX content changed after inventory')
            result = describe(read_cdx_molecule(data, as_single=True))
        except Exception as exc:
            result = {'status': 'error', 'error_type': type(exc).__name__,
                      'error': str(exc), 'cause': str(exc.__cause__) if exc.__cause__ else None}
        pipe.send(('result', (sha, result, time.monotonic() - start)))
    for z in cache.values():
        z.close()


def revision_info(repo):
    import importlib.metadata
    from rdkit import rdBase
    digest = hashlib.sha256()
    for path in sorted((Path(repo) / 'cdx_reader').glob('*.py')):
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    pycdxml = importlib.metadata.distribution('pycdxml')
    return {'reader_source_sha256': digest.hexdigest(), 'repo': str(Path(repo).resolve()),
            'rdkit': rdBase.rdkitVersion, 'python': sys.version,
            'pycdxml': pycdxml.version,
            'pycdxml_vcs': json.loads(pycdxml.read_text('direct_url.json') or '{}').get('vcs_info')}


def run(args):
    db = connect(args.db)
    info = revision_info(args.repo)
    info['timeout'] = args.timeout
    key = 'run:' + args.engine
    previous = db.execute('SELECT value FROM metadata WHERE key=?', (key,)).fetchone()
    if previous and json.loads(previous[0]) != info:
        raise SystemExit('Revision/environment changed: use a new engine name, preserving old results.')
    db.execute('INSERT OR REPLACE INTO metadata VALUES (?,?)', (key, json.dumps(info)))
    db.commit()
    # Workers always import an immutable snapshot. Editing the working checkout
    # during a long run cannot change the code used by a replacement worker.
    snapshot = Path(args.db).resolve().parent / 'reader-snapshots' / info['reader_source_sha256']
    package = snapshot / 'cdx_reader'
    package.mkdir(parents=True, exist_ok=True)
    for source in (Path(args.repo) / 'cdx_reader').glob('*.py'):
        target = package / source.name
        if not target.exists():
            shutil.copyfile(source, target)
        if target.read_bytes() != source.read_bytes():
            raise SystemExit('Source snapshot has been modified: ' + str(target))
    tasks = db.execute('''SELECT sha,zip,member FROM files WHERE sha NOT IN
                       (SELECT sha FROM results WHERE engine=?) ORDER BY zip,member''', (args.engine,)).fetchall()
    if args.limit:
        tasks = tasks[:args.limit]
    print('run', args.engine, 'pending', len(tasks), json.dumps(info), flush=True)
    if not tasks:
        return
    ctx = mp.get_context('spawn')
    slots = {}
    pending = iter(tasks)
    done = 0
    buffered = []

    def flush():
        if buffered:
            db.executemany('INSERT INTO results VALUES (?,?,?,?)', buffered)
            db.commit()
            buffered.clear()

    def start_slot():
        parent, child = ctx.Pipe()
        process = ctx.Process(target=worker, args=(child, str(snapshot)))
        process.start()
        child.close()
        slots[parent] = {'process': process, 'task': None, 'since': time.monotonic()}

    def finish(sha, result, seconds):
        nonlocal done
        buffered.append((args.engine, sha, json.dumps(result), seconds))
        done += 1
        if done % 100 == 0:
            flush()
        if done % 1000 == 0:
            print(args.engine, done, '/', len(tasks), flush=True)

    def dispatch(pipe):
        try:
            task = next(pending)
        except StopIteration:
            pipe.send(None)
            slot = slots.pop(pipe)
            slot['process'].join(5)
            if slot['process'].is_alive():
                slot['process'].terminate()
                slot['process'].join()
            pipe.close()
            return
        slots[pipe].update(task=task, since=time.monotonic())
        pipe.send(task)

    for _ in range(min(args.workers, len(tasks))):
        start_slot()
    try:
        while slots:
            for pipe in wait(list(slots), timeout=0.25):
                slot = slots[pipe]
                try:
                    event, payload = pipe.recv()
                except EOFError:
                    event, payload = 'crash', None
                if event == 'ready':
                    dispatch(pipe)
                elif event == 'result':
                    sha, result, seconds = payload
                    if sha != slot['task'][0]:
                        raise RuntimeError('Worker result does not match assigned input')
                    finish(sha, result, seconds)
                    dispatch(pipe)
                else:
                    if slot['task'] is None:
                        raise RuntimeError('Worker failed during initialization')
                    finish(slot['task'][0], {'status': 'crash'}, time.monotonic() - slot['since'])
                    slot['process'].join(5)
                    pipe.close()
                    del slots[pipe]
                    start_slot()
            for pipe, slot in list(slots.items()):
                budget = args.timeout if slot['task'] else max(60, args.timeout)
                if time.monotonic() - slot['since'] <= budget:
                    continue
                slot['process'].terminate()
                slot['process'].join(5)
                if slot['process'].is_alive():
                    slot['process'].kill()
                    slot['process'].join()
                if slot['task'] is None:
                    raise RuntimeError('Worker initialization timed out')
                finish(slot['task'][0], {'status': 'timeout'}, time.monotonic() - slot['since'])
                pipe.close()
                del slots[pipe]
                start_slot()
    finally:
        for pipe, slot in slots.items():
            slot['process'].terminate()
            slot['process'].join()
            pipe.close()
        flush()
    print('completed', args.engine, done, flush=True)


CHEMISTRY_FIELDS = ('status', 'smiles', 'graph_smiles', 'inchikey', 'formula',
                    'charge', 'components', 'radicals', 'harvest_accepted',
                    'error_type', 'error', 'cause')


def compare(args):
    db = connect(args.db)
    before_table = 'results'
    if args.before_db:
        db.execute('ATTACH DATABASE ? AS previous', (args.before_db,))
        before_table = 'previous.results'
    before_metadata = 'previous.metadata' if args.before_db else 'metadata'
    before_info = db.execute(f'SELECT value FROM {before_metadata} WHERE key=?',
                             ('run:' + args.before,)).fetchone()
    after_info = db.execute('SELECT value FROM metadata WHERE key=?', ('run:' + args.after,)).fetchone()
    if not before_info or not after_info:
        raise SystemExit('Missing run metadata: both revisions must be recorded.')
    before_info, after_info = json.loads(before_info[0]), json.loads(after_info[0])
    if any(before_info.get(k) != after_info.get(k)
           for k in ('python', 'rdkit', 'pycdxml', 'pycdxml_vcs', 'timeout')):
        raise SystemExit('Cannot compare runs with different dependency versions or time limits.')
    total = db.execute('SELECT COUNT(*) FROM files').fetchone()[0]
    counts = collections.Counter()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    with (output / 'changes.jsonl').open('w') as changed:
        for sha, before, after in db.execute(f'''
            SELECT a.sha,a.payload,b.payload FROM {before_table} a JOIN results b ON a.sha=b.sha
            WHERE a.engine=? AND b.engine=? ORDER BY a.sha''', (args.before, args.after)):
            a, b = json.loads(before), json.loads(after)
            counts['compared'] += 1
            counts['before_' + a['status']] += 1
            counts['after_' + b['status']] += 1
            counts['before_accepted'] += bool(a.get('harvest_accepted'))
            counts['after_accepted'] += bool(b.get('harvest_accepted'))
            chemistry_changed = any(a.get(k) != b.get(k) for k in CHEMISTRY_FIELDS)
            metadata_changed = a.get('labels') != b.get('labels')
            if chemistry_changed:
                counts['chemistry_or_status_changed'] += 1
            elif metadata_changed:
                counts['metadata_only'] += 1
            else:
                counts['unchanged'] += 1
            if a.get('harvest_accepted') and not b.get('harvest_accepted'):
                counts['accepted_to_rejected'] += 1
            if not a.get('harvest_accepted') and b.get('harvest_accepted'):
                counts['rejected_to_accepted'] += 1
            if chemistry_changed or metadata_changed:
                changed.write(json.dumps({'sha': sha, 'chemistry_changed': chemistry_changed,
                                          'before': a, 'after': b}) + '\n')
    counts['inventory'] = total
    counts['missing_pairs'] = total - counts['compared']
    summary = dict(counts)
    (output / 'summary.json').write_text(json.dumps(summary, indent=2))
    (output / 'environment.json').write_text(json.dumps({'before': before_info, 'after': after_info}, indent=2))
    print(json.dumps(summary, indent=2))
    if counts['missing_pairs']:
        raise SystemExit('Comparison incomplete: missing results.')


def positive_number(value):
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError('must be positive')
    return number


def positive_integer(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError('must be positive')
    return number


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('inventory')
    p.add_argument('--db', required=True)
    p.add_argument('inputs', nargs='+')
    p.set_defaults(func=inventory)
    p = sub.add_parser('run')
    p.add_argument('--db', required=True)
    p.add_argument('--repo', required=True)
    p.add_argument('--engine', required=True)
    p.add_argument('--workers', type=positive_integer, default=6)
    p.add_argument('--timeout', type=positive_number, default=30)
    p.add_argument('--limit', type=positive_integer)
    p.set_defaults(func=run)
    p = sub.add_parser('compare')
    p.add_argument('--db', required=True)
    p.add_argument('--before', default='baseline')
    p.add_argument('--before-db', help='Read baseline from a separate inventory/result database')
    p.add_argument('--after', default='candidate')
    p.add_argument('--output', required=True)
    p.set_defaults(func=compare)
    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main()
