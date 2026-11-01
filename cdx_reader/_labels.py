"""Preserve undecoded text as uncertainty, never as ordinary carbon or a port."""
import json
import re
from rdkit import Chem
from ._normalization import atom_label, _bonds_to


UNPARSED = '_CDX_UNPARSED_LABEL'
ANNOTATIONS = '_CDX_ANNOTATIONS'
_ANNOTATION = re.compile(r'(?:\((?:[RSEZ](?:[,/][RSEZ])*|[IVX]+)\)|and|or|\+|'
                         r'enantiomer\s*\d+|diastereomer\s*\d+|racemic)', re.I)


def simple_element_label(label):
    """An element with optional explicitly written H and charge, no prefixes."""
    label = re.sub(r'\s+', '', label).replace('−', '-')
    match = re.fullmatch(r'(?:(H[1-9]?)?([A-Z][a-z]?)(H[1-9]?)?)([+-][1-9]?|[1-9][+-])?', label)
    if not match:
        return None
    before, symbol, after, charge = match.groups()
    if before and after:
        return None
    if symbol not in {'H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I'}:
        return None
    hydrogen = before or after
    h_count = (int(hydrogen[1:] or '1') if hydrogen else None)
    if charge:
        sign = -1 if '-' in charge else 1
        charge = sign * int(charge.replace('+', '').replace('-', '') or '1')
    return symbol, h_count, charge


def mark_unparsed(node, label):
    node.set('_CDXUnparsedLabel', label)
    node.set('Element', '0')
    node.set('NodeType', 'Unspecified')


def restore_unparsed(mol, labels, declared_variables=()):
    if mol is not None:
        for atom in mol.GetAtoms():
            if atom.HasProp('_CDX_ATOM_ID'):
                text = labels.get(atom.GetProp('_CDX_ATOM_ID'))
                identifier = atom.GetProp('_CDX_ATOM_ID')
                if (text is None and atom.GetAtomicNum() == 0
                        and identifier not in declared_variables
                        and atom.HasProp('_CDX_ATOM_LABEL')
                        and atom.GetProp('_CDX_ATOM_LABEL').strip() != '*'
                        and not atom.HasProp('_CDX_ATTACHMENT_REPAIR')):
                    # RDKit can leave an opaque nested Fragment/Nickname as a
                    # dummy. Its text is not thereby a declared variable.
                    text = atom.GetProp('_CDX_ATOM_LABEL')
                if text is not None:
                    atom.SetProp(UNPARSED, text)
                    atom.SetProp('_CDX_LABEL_STATUS', 'unparsed')
    return mol


def detach_annotations(tree):
    """Detach only recognized, disconnected text nodes. Keep original object IDs."""
    annotations, empty_objects = [], []
    root = tree.getroot()
    for fragment in list(root.iter('fragment')):
        if any(a.tag in ('fragment', 'n') for a in fragment.iterancestors()):
            continue
        local = []
        for node in list(fragment.findall('n')):
            label = atom_label(node)
            if (node.get('NodeType') not in ('Unspecified', 'Fragment') or node.get('Element') not in (None, '0')
                    or node.get('Charge', '0') != '0' or node.get('Isotope') is not None
                    or node.get('Radical') is not None or _bonds_to(fragment, node.get('id'))
                    or not _ANNOTATION.fullmatch(label)):
                continue
            entry = {'text': label, 'cdx_atom_id': node.get('id'), 'fragment_id': fragment.get('id')}
            local.append(entry)
            annotations.append(entry)
            fragment.remove(node)
        if local and not fragment.findall('n') and not fragment.findall('b'):
            empty_objects.append((fragment.get('id'), local))
            fragment.getparent().remove(fragment)
    return annotations, empty_objects


def annotation_object(identifier, entries):
    mol = Chem.Mol()
    mol.SetProp('CDXML_FRAG_ID', identifier)
    mol.SetProp('_CDX_OBJECT_KIND', 'annotation')
    mol.SetProp(ANNOTATIONS, json.dumps(entries, ensure_ascii=False))
    return mol
