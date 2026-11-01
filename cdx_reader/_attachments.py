"""Recognize a drawn two-ended wavy attachment bar before CDXML loses display data.

This is deliberately not a general interpretation of wavy bonds. A bar has two
bare terminal carbon nodes opposite each other and an ordinary stem. A stem may
continue beyond the bar to a unique bare terminal stub opposite the substituent.
"""
import json
import math

from ._normalization import atom_label


def _bare_carbon(node):
    return (node is not None and node.get('Element') in (None, '6')
            and node.get('NodeType') in (None, 'Element')
            and not atom_label(node) and node.find('fragment') is None
            and node.get('Charge') in (None, '0')
            and node.get('Isotope') in (None, '0')
            and node.get('Radical') in (None, '0', 'None')
            and node.get('ElementList') is None and node.get('GenericList') is None)


def _opposite(centre, first, second):
    try:
        c, a, b = [tuple(map(float, n.get('p', '').split())) for n in (centre, first, second)]
        if any(len(p) != 2 or not all(map(math.isfinite, p)) for p in (c, a, b)):
            return False
        u, v = (a[0] - c[0], a[1] - c[1]), (b[0] - c[0], b[1] - c[1])
        length = math.hypot(*u) * math.hypot(*v)
        # At most about eight degrees from antiparallel. No stem-angle condition:
        # the source drawings include oblique bars.
        return length > 0 and sum(x * y for x, y in zip(u, v)) / length <= -0.99
    except (TypeError, ValueError):
        return False


def normalize_wavy_attachment_points(tree):
    """Modify only confirmed bars; return centre CDX ID -> removed tail IDs.

    Work within each fragment so a nested abbreviation cannot donate nodes or
    bonds to its parent. IDs are source IDs, never newly invented attachment maps.
    """
    return _normalize_attachment_points(tree, 'Wavy', crossed=True)


def normalize_bold_attachment_points(tree):
    """The same two-tail bar drawn in bold, with one ordinary stem.

    Bold bonds inside a skeleton are untouched. Unlike the reviewed wavy
    crossed-stem convention, no four-neighbour bold case is inferred here.
    """
    return _normalize_attachment_points(tree, 'Bold', crossed=False)


def _normalize_attachment_points(tree, display, crossed):
    repairs = {}
    for fragment in tree.getroot().iter('fragment'):
        nodes = {n.get('id'): n for n in fragment.findall('n')}
        if None in nodes or len(nodes) != len(fragment.findall('n')):
            continue
        edges = {identifier: [] for identifier in nodes}
        for bond in fragment.findall('b'):
            begin, end = bond.get('B'), bond.get('E')
            if begin in edges:
                edges[begin].append((bond, end))
            if end in edges:
                edges[end].append((bond, begin))
        selected = []
        for identifier, node in nodes.items():
            links = edges[identifier]
            if not _bare_carbon(node) or len(links) not in ((3, 4) if crossed else (3,)):
                continue
            waves = [(b, other) for b, other in links if b.get('Display') == display]
            stems = [(b, other) for b, other in links if b.get('Display') in (None, 'Solid')]
            if (len(waves) != 2 or len(stems) != len(links) - 2
                    or any(other not in nodes for _, other in stems)):
                continue
            if any(b.get('Order', '1') != '1' or b.get('Display2') not in (None, 'Solid')
                   for b, _ in links):
                continue
            tails = [other for _, other in waves]
            if len(set(tails)) != 2 or any(not _bare_carbon(nodes.get(t))
                                           or len(edges.get(t, [])) != 1 for t in tails):
                continue
            if not _opposite(node, nodes[tails[0]], nodes[tails[1]]):
                continue
            removed = list(waves)
            if len(stems) == 2:
                # A crossed stem has exactly one dispensable end. If both ordinary
                # ends are bare terminals, or both lead into graphs, do not guess
                # which side of the bar belongs to the substituent.
                stubs = [(b, other) for b, other in stems
                         if _bare_carbon(nodes[other]) and len(edges[other]) == 1]
                if (len(stubs) != 1 or len({other for _, other in links}) != 4
                        or not _opposite(node, nodes[stems[0][1]], nodes[stems[1][1]])):
                    continue
                removed += stubs
            selected.append((identifier, node, removed))
        for identifier, node, removed in selected:
            node.set('Element', '0')
            node.set('NumHydrogens', '0')
            for bond, tail in removed:
                fragment.remove(bond)
                fragment.remove(nodes[tail])
            repairs[identifier] = [tail for _, tail in removed]
    return repairs


def restore_attachment_provenance(mol, repairs, style='wavy'):
    if mol is not None:
        for atom in mol.GetAtoms():
            if atom.HasProp('_CDX_ATOM_ID'):
                identifier = atom.GetProp('_CDX_ATOM_ID')
                if identifier in repairs and atom.GetAtomicNum() == 0:
                    operation = (f'{style}_crossed_stem_bar' if len(repairs[identifier]) == 3
                                 else f'{style}_two_tail_bar')
                    atom.SetProp('_CDX_ATTACHMENT_REPAIR', operation)
                    atom.SetProp('_CDX_ATTACHMENT_TAIL_IDS', json.dumps(repairs[identifier]))
    return mol
