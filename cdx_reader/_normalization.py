"""Preserve source variable labels without using them as joining IDs."""
import re

from rdkit import Chem

R_LABEL = re.compile(r'R(?:\d+|[a-g])?\Z')

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
        # Indexed non-element letters are frequent alternate variable names.
        # Ordinary elemental labels and H-count formulas retain their chemistry.
        if (re.fullmatch(r'[A-Z]\d+', label)
                and node.get('Element') in (None, '0', '99', '100')):
            return label
        if label not in ELEMENT_SYMBOLS and re.fullmatch(r'[ABEGJKLMQUVWXYZ](?:[a-z])?', label):
            return label
    return None

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
