"""Read a cyano label "NC" the way it is drawn, when ChemDraw expanded it the other way round.

An author writing a cyano group to the left of a ring types "NC", so that the C stands next to
the bond. ChemDraw expands the label by its letters, N first, and stores an isocyanide bonded
through N (external-N(+)#C(-)); its MOL export says the same. The picture and the patent's name
(...carbonitrile) say nitrile. Only that exact case is changed: the label text NC displayed to
the right of its box, one single bond to a neighbour on the right, and a nested graph of exactly
the external point, N(+) and C(-). An isocyanide drawn with its own N+ and C- atoms, or a label
whose N is next to the bond, is left as encoded. Each change is returned, so the reading can be
shown with what the file encodes.
"""


def _position(node):
    try:
        x, y = map(float, node.get('p', '').split())
        return x, y
    except ValueError:
        return None


def restore_label_nitriles(tree):
    """Turn each such NC label's nested N(+)#C(-) into a neutral C#N bonded through C; return
    the changes as {'cdx_atom_id', 'fragment_id', 'label', 'encoded', 'read'}."""
    root = tree.getroot()
    repairs = []
    for node in list(root.iter('n')):
        if ''.join(node.xpath('./t//s/text()')) != 'NC' or node.get('LabelDisplay') != 'Right':
            continue
        nested = node.findall('fragment')
        if len(nested) != 1:
            continue
        atoms = {a.get('id'): a for a in nested[0].findall('n')}
        bonds = nested[0].findall('b')
        external = [i for i, a in atoms.items() if a.get('NodeType') == 'ExternalConnectionPoint']
        nitrogen = [i for i, a in atoms.items() if a.get('Element') == '7' and a.get('Charge') == '1']
        carbon = [i for i, a in atoms.items() if a.get('Element', '6') == '6' and a.get('Charge') == '-1'
                  and a.get('NodeType') != 'ExternalConnectionPoint']
        if len(atoms) != 3 or len(bonds) != 2 or not len(external) == len(nitrogen) == len(carbon) == 1:
            continue
        if any(a.get('Isotope') is not None or a.get('Radical') is not None for a in atoms.values()):
            continue
        edges = {frozenset((b.get('B'), b.get('E'))): b.get('Order', '1') for b in bonds}
        if edges != {frozenset((external[0], nitrogen[0])): '1', frozenset((nitrogen[0], carbon[0])): '3'}:
            continue
        parent = node.getparent()
        outer = [b for b in parent.findall('b') if node.get('id') in (b.get('B'), b.get('E'))]
        if len(outer) != 1 or outer[0].get('Order', '1') != '1':
            continue
        other_id = outer[0].get('E') if outer[0].get('B') == node.get('id') else outer[0].get('B')
        other = next((a for a in parent.findall('n') if a.get('id') == other_id), None)
        here, there = _position(node), _position(other) if other is not None else None
        if here is None or there is None or there[0] <= here[0]:
            continue
        # The external point now bonds to C, and N is the terminal atom: swap the two elements.
        for identifier, element, text in ((nitrogen[0], '6', 'C'), (carbon[0], '7', 'N')):
            atom = atoms[identifier]
            atom.set('Element', element)
            atom.set('Charge', '0')
            for at, piece in enumerate(atom.findall('./t/s')):
                piece.text = text if at == 0 else ''
        outer_fragment = next((a for a in node.iterancestors() if a.tag == 'fragment'), None)
        repairs.append({'cdx_atom_id': node.get('id'),
                        'fragment_id': outer_fragment.get('id') if outer_fragment is not None else None,
                        'label': 'NC', 'encoded': 'isocyanide, N(+)#C(-) bonded through N',
                        'read': 'nitrile, C#N bonded through C'})
    return repairs
