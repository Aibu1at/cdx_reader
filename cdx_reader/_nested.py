"""Flatten multi-attachment CDXML groups whose external ends share one atom.

RDKit's CDXML molzip can lose these fragments before returning a molecule. No
orientation is inferred here: every external bond already points at the same
internal atom. Disconnected components (such as a counterion) are retained.
Non-equivalent ends and query graphs stay untouched.
"""
import json
from rdkit import Chem


def flatten_shared_anchor_groups(tree):
    records = {}
    for wrapper in list(tree.getroot().iter('n')):
        children = wrapper.findall('fragment')
        if len(children) != 1 or wrapper.get('NodeType') != 'Fragment':
            continue
        fragment = children[0]
        nodes, bonds = fragment.findall('n'), fragment.findall('b')
        ids = {n.get('id') for n in nodes}
        external = {n.get('id') for n in nodes if n.get('NodeType') == 'ExternalConnectionPoint'}
        if len(external) < 2 or None in ids or len(ids) != len(nodes):
            continue
        if any(n.find('fragment') is not None or n.get('ElementList') is not None
               or n.get('GenericList') is not None or n.get('Radical') is not None
               or n.get('NodeType') not in (None, 'Element', 'ExternalConnectionPoint') for n in nodes):
            continue
        edges = {i: [] for i in ids}
        valid = True
        for b in bonds:
            a, z = b.get('B'), b.get('E')
            if a not in ids or z not in ids or a == z or b.get('Order', '1') not in ('1','2','3'):
                valid = False
                break
            edges[a].append(z); edges[z].append(a)
        if not valid or any(len(edges[e]) != 1 for e in external):
            continue
        anchors = {edges[e][0] for e in external}
        if len(anchors) != 1 or anchors & external:
            continue
        anchor = next(iter(anchors))
        outer = wrapper.getparent()
        outside = [b for b in outer.findall('b') if wrapper.get('id') in (b.get('B'),b.get('E'))]
        external_bonds = [b for b in bonds if b.get('B') in external or b.get('E') in external]
        if (len(outside) != len(external)
                or any(b.get('Order','1') != '1' for b in outside + external_bonds)):
            continue
        # Validate the specified graph, without the legacy valence fixer.
        try:
            probe = Chem.RWMol(); mapped = {}
            for n in nodes:
                a = Chem.Atom(0 if n.get('id') in external else int(n.get('Element','6')))
                a.SetFormalCharge(int(n.get('Charge','0')))
                a.SetIsotope(int(n.get('Isotope','0')))
                if n.get('NumHydrogens') is not None:
                    a.SetNumExplicitHs(int(n.get('NumHydrogens'))); a.SetNoImplicit(True)
                mapped[n.get('id')] = probe.AddAtom(a)
            for b in bonds:
                probe.AddBond(mapped[b.get('B')],mapped[b.get('E')],
                              {'1':Chem.BondType.SINGLE,'2':Chem.BondType.DOUBLE,'3':Chem.BondType.TRIPLE}[b.get('Order','1')])
            Chem.SanitizeMol(probe)
            if any(a.GetNumRadicalElectrons() for a in probe.GetAtoms()):
                continue
        except (ValueError, RuntimeError, OverflowError):
            continue
        source = wrapper.get('id')
        records[source] = dict(anchor=anchor, atoms=sorted(ids-external),
                               external=sorted(external), operation='flatten_shared_anchor')
        for b in outside:
            for end in ('B','E'):
                if b.get(end) == source:
                    b.set(end,anchor)
                    b.attrib.pop('BeginAttach' if end == 'B' else 'EndAttach',None)
        for n in list(fragment):
            if n.get('id') not in external and n not in external_bonds:
                outer.append(n)
        outer.remove(wrapper)
        # RDKit reads CDXML sequentially: endpoints must precede their bonds.
        for b in outer.findall('b'):
            outer.remove(b); outer.append(b)
    return records


def restore_nested_provenance(mol, records):
    if mol is None or not records:
        return mol
    present = {a.GetProp('_CDX_ATOM_ID'):a for a in mol.GetAtoms() if a.HasProp('_CDX_ATOM_ID')}
    applied = {}
    for source, record in records.items():
        if record['anchor'] not in present:
            continue
        applied[source] = record
        for identifier in record['atoms']:
            if identifier in present:
                present[identifier].SetProp('_CDX_GROUP_SOURCE_ATOM', source)
    if applied:
        mol.SetProp('_CDX_NESTED_EXPANSIONS', json.dumps(applied, sort_keys=True))
    return mol
