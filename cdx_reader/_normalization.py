"""Conservative CDXML repairs and preservation of existing R-group labels."""
import re
from rdkit import Chem


R_LABEL = re.compile(r"R(?:\d+[a-z]?|[a-g]|[′']+)?\Z")
ELEMENT_SYMBOLS = {Chem.GetPeriodicTable().GetElementSymbol(number) for number in range(1, 119)}


def atom_label(node):
    """Read all styled segments, including separately styled R-group numbers."""
    text = node.find('t')
    if text is None:
        return ''
    return ''.join(''.join(segment.itertext()) for segment in text.iter('s')).strip()


def r_label(node):
    label = re.sub(r'\s+', '', atom_label(node))
    return label if R_LABEL.fullmatch(label) else None


def variable_label(node, abbreviations):
    """Recognize declared placeholders without guessing from an element letter."""
    label = re.sub(r'\s+', '', atom_label(node))
    if r_label(node):
        return label
    if not label or label in abbreviations:
        return None
    if node.get('NodeType') == 'GenericNickname':
        return label
    if node.get('NodeType') == 'Unspecified':
        # US20210122713A1 Table 3 writes Ar with a separately styled A/B.
        # ChemDraw also cached an invalid argon graph under these text nodes.
        # Only these source-confirmed labels are added: preserve any explicitly
        # encoded element or usable nested chemical graph.
        if (label in ('ArA', 'ArB') and node.get('Element') in (None, '0')
                and node.get('ElementList') is None
                and node.get('GenericList') is None
                and node.get('Charge') in (None, '0')
                and node.get('Isotope') is None
                and node.get('Radical') is None
                and not has_unambiguous_nested_fragment(node)):
            return label
        # Indexed non-element letters are frequent alternate variable names.
        # Ordinary elemental labels and H-count formulas retain their chemistry.
        if (re.fullmatch(r'[A-Z]\d+', label)
                and node.get('Element') in (None, '0', '99', '100')):
            return label
        if label not in ELEMENT_SYMBOLS and re.fullmatch(r'[ABEGJKLMQUVWXYZ](?:[a-z])?', label):
            return label
    return None


def _bonds_to(parent, node_id):
    return [bond for bond in parent.findall('b')
            if node_id in (bond.get('B'), bond.get('E'))]


def has_unambiguous_nested_fragment(node):
    """Recognize a connected nested graph with exactly one attachment point.

    Multi-attachment groups and query graphs need their own interpretation and
    are intentionally not repaired by this rule.
    """
    fragments = node.findall('fragment')
    if len(fragments) != 1:
        return False
    fragment = fragments[0]
    nodes = fragment.findall('n')
    ids = {child.get('id') for child in nodes}
    if None in ids or len(ids) != len(nodes) or len(ids) < 2:
        return False
    external = [child for child in nodes if child.get('NodeType') == 'ExternalConnectionPoint']
    if len(external) != 1:
        return False
    allowed = {None, 'Element', 'Unspecified', 'Fragment', 'ExternalConnectionPoint'}
    if any(child.get('NodeType') not in allowed or child.get('ElementList') is not None
           or child.get('GenericList') is not None or child.get('Radical') is not None
           or variable_label(child, set()) is not None
           for child in nodes):
        return False
    edges = {identifier: set() for identifier in ids}
    for bond in fragment.findall('b'):
        begin, end = bond.get('B'), bond.get('E')
        if begin not in ids or end not in ids or begin == end:
            return False
        edges[begin].add(end)
        edges[end].add(begin)
    external_id = external[0].get('id')
    if len(edges[external_id]) != 1:
        return False
    visited, pending = set(), [external_id]
    while pending:
        current = pending.pop()
        if current not in visited:
            visited.add(current)
            pending.extend(edges[current] - visited)
    parent = node.getparent()
    if visited != ids or parent is None or len(_bonds_to(parent, node.get('id'))) != 1:
        return False
    internal_attachment = _bonds_to(fragment, external_id)
    outer_attachment = _bonds_to(parent, node.get('id'))
    if (len(internal_attachment) != 1
            or internal_attachment[0].get('Order', '1') != outer_attachment[0].get('Order', '1')):
        return False
    # A connected drawing can still contain impossible chemistry (e.g. oxygen
    # bonded to three fluorines, or CH3 at two attachment bonds). Do not promote
    # those ambiguous drawings merely because an internal fragment exists.
    probe = Chem.RWMol()
    atom_indices = {}
    try:
        for child in nodes:
            if child.find('fragment') is not None or child.get('NodeType') == 'Fragment':
                return False  # recursively nested attachments need a separate rule
            if (child.get('NodeType') == 'Unspecified' and child.get('Element') is None
                    and atom_label(child) not in ('', 'C', 'CH', 'CH2', 'CH3')):
                return False  # a text placeholder is not a fully specified graph
            atomic_number = 0 if child is external[0] else int(child.get('Element', '6'))
            atom = Chem.Atom(atomic_number)
            atom.SetFormalCharge(int(child.get('Charge', '0')))
            atom.SetIsotope(int(child.get('Isotope', '0')))
            if child.get('NumHydrogens') is not None:
                atom.SetNumExplicitHs(int(child.get('NumHydrogens')))
                atom.SetNoImplicit(True)
            atom_indices[child.get('id')] = probe.AddAtom(atom)
        orders = {'1': Chem.BondType.SINGLE, '2': Chem.BondType.DOUBLE,
                  '3': Chem.BondType.TRIPLE}
        for bond in fragment.findall('b'):
            if bond.get('Order', '1') not in orders:
                return False
            probe.AddBond(atom_indices[bond.get('B')], atom_indices[bond.get('E')],
                          orders[bond.get('Order', '1')])
        Chem.SanitizeMol(probe)
        if any(a.GetNumRadicalElectrons() for a in probe.GetAtoms()):
            return False
        if any(b.GetBondType() not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE,
                                      Chem.BondType.TRIPLE, Chem.BondType.AROMATIC)
               for b in probe.GetBonds()):
            # Sanitization may silently turn an overbonded halogen next to a
            # metal into a dative bond. Such a rewritten graph is not evidence
            # that the original nested attachment was chemically valid.
            return False
    except (ValueError, RuntimeError, OverflowError):
        return False
    return True


def explicit_element_is_usable(node):
    """Preserve explicit atoms without reinterpreting ambiguous variables.

    ChemDraw also assigns element numbers to variable letters such as U/V/K,
    and encodes T variables as tritium. Their general interpretation is outside
    these repairs; charge/isotope annotations distinguish explicit metal atoms.
    """
    label = atom_label(node)
    if (label == 'Sn' and node.get('Element') == '50' and node.get('Charge', '0') == '0'
            and node.get('NumHydrogens') == '0' and node.get('Radical') is None):
        parent = node.getparent()
        bonds = _bonds_to(parent, node.get('id')) if parent is not None else []
        if (bonds and all(b.get('Order', '1') in ('1', '2', '3') for b in bonds)
                and sum(int(b.get('Order', '1')) for b in bonds) < 4):
            # Low-valent Sn needs electron/coordination handling: RDKit's
            # default radical assignment does not represent Sn(II) correctly.
            return False
    if label == 'K' and 'invalid valence' in node.get('Warning', '').lower():
        # RDKit 2025.03.3 turns a warned covalent K+ salt drawing into a dative
        # complex with seven radical electrons on K. Leave that case unresolved
        # under the existing parser until salt/coordination handling is repaired.
        return False
    if (node.get('Element') == '11' and node.get('Charge', '0') == '0'
            and 'invalid valence' in node.get('Warning', '').lower()):
        parent = node.getparent()
        bonds = _bonds_to(parent, node.get('id')) if parent is not None else []
        # A documented typo places "Na" at an aromatic N position, with a
        # single and a double bond. Preserve the existing N interpretation of
        # that inconsistent drawing, while retaining genuine sodium atoms.
        if len(bonds) > 1 or any(b.get('Order', '1') in ('2', '3') for b in bonds):
            return False
    if (label in ('U', 'V', 'K') and node.get('Charge', '0') == '0'
            and node.get('Isotope', '0') == '0'):
        return False
    if node.get('Element') == '1' and label in ('D', 'T'):
        parent = node.getparent()
        if parent is None or node.get('Warning'):
            return False
        bonds = _bonds_to(parent, node.get('id'))
        return len(bonds) == 1 and bonds[0].get('Order', '1') == '1'
    return True


def has_ambiguous_multimolecule_expansion(root, abbreviations):
    """Do not alter selection when repairs interact with unexpanded groups.

    The existing abbreviation expander may choose one component from a drawing.
    Promoting a nested group in another component can change that choice. Defer
    nested repairs in these mixed drawings until selection has an explicit API.
    """
    drawings = [fragment for fragment in root.iter('fragment')
                if not any(parent.tag == 'n' for parent in fragment.iterancestors())]
    return len(drawings) > 1 and any(
        atom_label(node).replace('\n', '').replace('\r', '') in abbreviations
        and node.find('fragment') is None for node in root.iter('n'))


def normalize_nested_fragment(node):
    if node.get('NodeType') not in (None, 'Unspecified', 'Fragment'):
        return False
    if not has_unambiguous_nested_fragment(node):
        return False
    node.set('NodeType', 'Fragment')
    return True


def repair_invalid_nitrogen_hydrogens(node):
    """Discard a contradictory zero-H hint, not explicit radical chemistry."""
    if (node.get('Element') != '7' or node.get('NumHydrogens') != '0'
            or node.get('Charge', '0') != '0'
            or node.get('NodeType') not in (None, 'Element')
            or not any('invalid valence' in item.get('Warning', '').lower()
                       for item in [node, *node.findall('t')])
            or node.get('Radical') is not None
            or node.get('ElementList') is not None
            or node.get('GenericList') is not None
            or node.find('fragment') is not None):
        return False
    parent = node.getparent()
    if parent is None:
        return False
    bonds = _bonds_to(parent, node.get('id'))
    if len(bonds) != 2 or any(bond.get('Order', '1') != '1' for bond in bonds):
        return False
    del node.attrib['NumHydrogens']
    return True


def normalize_hollow_wedges(tree):
    """Translate explicit hollow single-bond wedges to RDKit-supported styles.

    Preserve the begin/end direction and retain the encoded style as provenance.
    Bold, Hash and Wavy have other meanings and are deliberately not converted.
    """
    styles = {'HollowWedgeBegin': 'WedgeBegin', 'HollowWedgeEnd': 'WedgeEnd'}
    repairs = []
    for bond in tree.getroot().iter('b'):
        encoded = bond.get('Display')
        if encoded not in styles or bond.get('Order', '1') != '1':
            continue
        fragments = [p for p in bond.iterancestors() if p.tag == 'fragment']
        if not fragments:
            continue
        repairs.append(dict(fragment_id=fragments[-1].get('id'), bond_id=bond.get('id'),
                            begin=bond.get('B'), end=bond.get('E'),
                            encoded=encoded, interpreted=styles[encoded]))
        bond.set('Display', styles[encoded])
    return repairs


def capture_r_labels(tree):
    return {node.get('id'): (atom_label(node), r_label(node))
            for node in tree.getroot().iter('n')
            if node.get('id') and r_label(node)}


def capture_atom_labels(tree):
    return {node.get('id'): (atom_label(node), node.get('NodeType', ''))
            for node in tree.getroot().iter('n') if node.get('id') and atom_label(node)}


def detach_variable_atom_maps(mol, labels):
    """Keep source dummy maps out of the fragment-joining namespace."""
    for atom in mol.GetAtoms():
        if (atom.GetAtomicNum() == 0 and atom.HasProp('_CDX_ATOM_ID')
                and atom.GetProp('_CDX_ATOM_ID') in labels and atom.GetAtomMapNum()):
            atom.SetIntProp('_CDX_ORIGINAL_ATOM_MAP', atom.GetAtomMapNum())
            atom.SetAtomMapNum(0)


def restore_atom_labels(mol, labels):
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        if not atom.HasProp('_CDX_ATOM_ID'):
            continue
        identifier = atom.GetProp('_CDX_ATOM_ID')
        if identifier not in labels:
            continue
        raw, node_type = labels[identifier]
        atom.SetProp('_CDX_ATOM_LABEL', raw)
        atom.SetProp('_CDX_NODE_TYPE', node_type)
        if atom.GetAtomicNum() == 0:
            atom.SetProp('dummyLabel', re.sub(r'\s+', '', raw))
    return mol


def restore_r_labels(mol, labels):
    """Attach metadata to the surviving original atom, never to an atom index."""
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 0 or not atom.HasProp('_CDX_ATOM_ID'):
            continue
        identifier = atom.GetProp('_CDX_ATOM_ID')
        if identifier not in labels:
            continue
        original, normalized = labels[identifier]
        atom.SetProp('_CDX_ATOM_LABEL', original)
        atom.SetProp('dummyLabel', normalized)
        if normalized[1:].isdigit():
            number = int(normalized[1:])
            # RDKit integer properties use signed 32-bit storage. Preserve the
            # textual label even when a malformed number exceeds that range.
            if number <= 2147483647:
                atom.SetIntProp('_MolFileRLabel', number)
    return mol


def sanitize_aromatic_variable_graph(mol, labels):
    """Retain an explicitly aromatic query edge without guessing its Kekule form.

    A ring variable can stand for carbon or nitrogen, so its electron donation
    is not defined. RDKit cannot always Kekulize such a graph. Work on the
    pristine input: failed Kekulization can already have changed bond orders.
    """
    if not any(bond.GetBondType() == Chem.BondType.AROMATIC and any(
            atom.GetAtomicNum() == 0 and atom.HasProp('_CDX_ATOM_ID')
            and atom.GetProp('_CDX_ATOM_ID') in labels
            for atom in (bond.GetBeginAtom(), bond.GetEndAtom()))
            for bond in mol.GetBonds()):
        return None
    probe = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(probe)
        return None  # This fallback is only for an actual Kekulization failure.
    except Chem.KekulizeException:
        pass
    except (ValueError, RuntimeError):
        return None
    probe = Chem.Mol(mol)
    try:
        Chem.SanitizeMol(probe, sanitizeOps=Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
    except (ValueError, RuntimeError):
        return None
    probe.SetBoolProp('_CDX_UNKEKULIZED_QUERY', True)
    return probe
