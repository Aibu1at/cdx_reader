"""Small structural grammar for reversible condensed functional-group labels.

This is not a general formula-to-structure guesser. Recognized families have a
single central attachment atom: alkoxy/amino, carbonyl, sulfonyl and halomethyl.
The original bond count validates the expansion. Valid embedded CDX graphs take
precedence, even when their chemistry differs from the visible abbreviation.
"""
import json
import re
from copy import deepcopy

from lxml import etree
from rdkit import Chem
from rdkit.Chem import rdDepictor

from ._normalization import atom_label, has_unambiguous_nested_fragment, variable_label

VARIABLE = r'R(?:\d+[a-z]?|[a-g])?'
UNIT = re.compile(r'(R\d+[a-z]?|R[a-g]?|Boc|Cbz|Fmoc|CH3|H3C|Me|Et|Ph|Bn|H)(\d*)')


def _leaf(text, library):
    if re.fullmatch(VARIABLE,text):
        return '*', text
    if text in ('CH3','H3C'):
        return 'C', None
    if text in ('Me','Et','Ph','Bn','Boc','Cbz','Fmoc') and text in library:
        mol = Chem.MolFromSmiles(library[text][1])
        ports = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        if len(ports) != 1 or ports[0].GetDegree() != 1:
            return None
        anchor = ports[0].GetNeighbors()[0].GetIdx()
        rw = Chem.RWMol(mol); removed = ports[0].GetIdx(); rw.RemoveAtom(removed)
        fragment = rw.GetMol(); fragment.UpdatePropertyCache(strict=False)
        return Chem.MolToSmiles(fragment,rootedAtAtom=anchor-(anchor>removed)), None
    return None


def _build(element, branches, degree, hydrogens=None):
    rw = Chem.RWMol(); center = Chem.Atom(element)
    if hydrogens is not None:
        center.SetNumExplicitHs(hydrogens); center.SetNoImplicit(True)
    rw.AddAtom(center)
    for smiles, label, order in branches:
        leaf = Chem.MolFromSmiles(smiles)
        if leaf is None:
            return None
        if label:
            ports = [a for a in leaf.GetAtoms() if a.GetAtomicNum() == 0]
            if len(ports) != 1:
                return None
            ports[0].SetProp('dummyLabel',label)
        start = rw.GetNumAtoms(); rw.InsertMol(leaf)
        rw.AddBond(0,start,order)
    for _ in range(degree):
        port = Chem.Atom(0); port.SetBoolProp('_external',True)
        rw.AddBond(0,rw.AddAtom(port),Chem.BondType.SINGLE)
    mol = rw.GetMol()
    try:
        Chem.SanitizeMol(mol)
        if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
            return None
        return mol
    except (ValueError, RuntimeError):
        return None


def parse_condensed(text, degree, library):
    """Return a graph with designated external ports, or None for an ambiguity."""
    single, double = Chem.BondType.SINGLE, Chem.BondType.DOUBLE
    if degree not in (1,2):
        return None
    # R-O-C(=O)- and H-O-C(=O)- in either condensed writing direction.
    for body in ('CO2','O2C'):
        tail = text[len(body):] if text.startswith(body) else text[:-len(body)] if text.endswith(body) else ''
        if degree == 1 and tail:
            leaf = ('',None) if tail == 'H' else _leaf(tail,library)
            if leaf is not None:
                return _build(6,[('O',None,double),('O'+leaf[0],leaf[1],single)],degree)
    # Carbonyl amide: NH2/H2N is a terminal group; CO/OC fixes the carbonyl.
    if degree == 1 and re.fullmatch(r'(?:(?:NH2|H2N)(?:CO|OC)|(?:CO|OC)(?:NH2|H2N))',text):
        return _build(6,[('O',None,double),('N',None,single)],degree)
    match = re.fullmatch(r'(.+)\(O\)C|C\(O\)(.+)',text)
    if degree == 1 and match:
        leaf = _leaf(match.group(1) or match.group(2),library)
        if leaf:
            return _build(6,[('O',None,double),(leaf[0],leaf[1],single)],degree)
    match = re.fullmatch(r'O([12])S|SO([12]?)',text)
    if degree == 2 and match:
        count = int(match.group(1) or match.group(2) or '1')
        return _build(16,[('O',None,double)]*count,degree)
    # Alkoxy/aryloxy/variable-O: counts such as O2 are not silently dropped.
    if degree == 1 and (text.startswith('O') or text.endswith('O')):
        tail = text[1:] if text.startswith('O') else text[:-1]
        leaf = _leaf(tail,library)
        if leaf:
            return _build(8,[(leaf[0],leaf[1],single)],degree)
    # Amino center and explicitly counted substituents/hydrogens, in either
    # direction. R12 is one label, never twelve copies of R.
    if text.startswith('N') or text.endswith('N'):
        tail = text[1:] if text.startswith('N') else text[:-1]
        pos, branches, hs, explicit_h = 0, [], 0, False
        while pos < len(tail):
            match = UNIT.match(tail,pos)
            if not match:
                return None
            unit, amount = match.groups(); amount = int(amount or 1)
            if amount > 3 or amount < 1:
                return None
            pos = match.end()
            if unit == 'H':
                hs += amount; explicit_h = True
            else:
                leaf = _leaf(unit,library)
                if leaf is None:
                    return None
                branches.extend([(leaf[0],leaf[1],single)]*amount)
        if branches and degree + len(branches) + hs <= 3:
            return _build(7,branches,degree,hs if explicit_h else None)
    # Counted halomethyl: the H count belongs to the central carbon.
    match = re.fullmatch(r'(F|Cl|Br)([123])C(H\d?)?|C(H\d?)?(F|Cl|Br)([123])',text)
    if degree == 1 and match:
        element = match.group(1) or match.group(5)
        count = int(match.group(2) or match.group(6))
        hydrogen = match.group(3) or match.group(4)
        hs = int(hydrogen[1:] or 1) if hydrogen else 0
        if count + hs + degree == 4:
            return _build(6,[(element,None,single)]*count,degree,hs)
    return None


def _write_label(node, value):
    for text in node.findall('t'):
        node.remove(text)
    if value:
        text = etree.SubElement(node,'t')
        etree.SubElement(text,'s').text = value


def _usable_nested_graph(node):
    if has_unambiguous_nested_fragment(node):
        return True
    # An already declared fragment may contain genuine variable atoms, e.g.
    # O-R2. Validate its topology with those variables as dummies, then keep the
    # original graph and its real atom IDs instead of rebuilding the same group.
    if node.get('NodeType') != 'Fragment' or len(node.findall('fragment')) != 1:
        return False
    parent = etree.Element('fragment')
    probe = deepcopy(node)
    parent.append(probe)
    for bond in node.getparent().findall('b'):
        if node.get('id') in (bond.get('B'), bond.get('E')):
            parent.append(deepcopy(bond))
    changed = False
    for child in probe.find('fragment').findall('n'):
        if variable_label(child, set()) is not None:
            child.set('NodeType', 'Element')
            child.set('Element', '0')
            child.set('NumHydrogens', '0')
            child.attrib.pop('GenericNickname', None)
            _write_label(child, '')
            changed = True
    return changed and has_unambiguous_nested_fragment(probe)


def expand_condensed_labels(tree, library):
    """Expand only undecoded text, preserving raw labels and generated-atom origins."""
    root = tree.getroot()
    maximum = max((int(n.get('id')) for n in root.iter() if (n.get('id') or '').isdigit()),default=0)
    generated, expansions, normalized = {}, {}, {}
    fonts = {f.get('id'): f.get('name','').lower() for f in root.iter('font')}
    for node in list(root.iter('n')):
        if node.get('NodeType') not in (None,'Unspecified','Fragment') or node.get('Element') not in (None,'0'):
            continue
        if node.get('ElementList') is not None or node.get('GenericList') is not None:
            continue
        if node.get('Charge','0') != '0' or node.get('Isotope','0') != '0' or node.get('Radical') is not None:
            continue
        raw = atom_label(node); text = re.sub(r'\s+','',raw)
        if not text:
            continue
        # Decode prime only in a Symbol-font run, not a literal cent sign.
        parts = []
        for t in node.findall('t'):
            for s in t.iter('s'):
                chunk = ''.join(s.itertext())
                if 'symbol' in fonts.get(s.get('font'),''):
                    chunk = chunk.replace('¢','′')
                parts.append(chunk)
        decoded = re.sub(r'\s+','', ''.join(parts))
        cleaned = decoded.rstrip(',;.')
        # These narrowly typed H-case errors were measured separately. An
        # explicitly encoded element (including Nh) never enters this path.
        if re.fullmatch(r'[NOS]h[1-3]?',cleaned) and node.get('Warning'):
            cleaned = cleaned[0]+'H'+cleaned[2:]
        parent = node.getparent()
        outside = [b for b in parent.findall('b') if node.get('id') in (b.get('B'),b.get('E'))]
        if not outside or any(b.get('Order','1') != '1' for b in outside):
            continue
        if _usable_nested_graph(node):
            continue
        if (re.fullmatch(r'R(?:\d+[a-z]?|[a-g]|[′\']+)?|[ABEGJKLMQUVWXYZ]\d*',cleaned)
                or re.fullmatch(r'(?:NH[0-3]?|OH[0-2]?|SH[0-2]?)',cleaned)):
            if cleaned != text and node.find('fragment') is None:
                _write_label(node,cleaned)
                normalized[node.get('id')] = cleaned
            continue
        # Existing dictionary groups without a nested graph keep their legacy
        # path except for multi-bond centers it cannot expand (BocN, MeN...).
        if text in library and len(outside) == 1:
            continue
        graph = parse_condensed(cleaned,len(outside),library)
        if graph is None:
            continue
        if node.get('Geometry') == 'Tetrahedral' or node.get('AS') in ('R','S','r','s'):
            continue
        source = node.get('id')
        try:
            # The CDXML writer stores a bond order as an integer: an aromatic 1.5 would become 1.
            Chem.Kekulize(graph, clearAromaticFlags=True)
        except Chem.KekulizeException:
            pass  # as before: written with the aromatic flags
        rdDepictor.Compute2DCoords(graph)
        conf = graph.GetConformer(); center_pos = conf.GetAtomPosition(0)
        try:
            origin = tuple(map(float,node.get('p','').split()))
            if len(origin) != 2:
                continue
        except ValueError:
            continue
        # All external ports in this grammar share atom 0. Existing outside
        # edges retain endpoints, order, wedge/hash information and IDs.
        mapped = {0:source}
        for f in node.findall('fragment'):
            node.remove(f)
        for k in ('NodeType','BondOrdering','GenericNickname','Warning'):
            node.attrib.pop(k,None)
        _write_label(node,'')
        created, variables = [], []
        for atom in graph.GetAtoms():
            if atom.HasProp('_external'):
                continue
            idx = atom.GetIdx()
            if idx == 0:
                target = node
            else:
                maximum += 1; identifier = str(maximum); mapped[idx] = identifier
                pos = conf.GetAtomPosition(idx)
                target = etree.SubElement(parent,'n',id=identifier,
                    p=f'{origin[0]+(pos.x-center_pos.x)*10:.4f} {origin[1]+(pos.y-center_pos.y)*10:.4f}')
                generated[identifier] = source; created.append(identifier)
            target.set('Element',str(atom.GetAtomicNum()))
            target.set('NumHydrogens',str(atom.GetTotalNumHs()))
            if atom.GetAtomicNum() == 0 and atom.HasProp('dummyLabel'):
                label = atom.GetProp('dummyLabel'); _write_label(target,label)
                target.set('NodeType','Unspecified'); variables.append(label)
        for bond in graph.GetBonds():
            a,b = bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
            if a in mapped and b in mapped:
                maximum += 1
                etree.SubElement(parent,'b',id=str(maximum),B=mapped[a],E=mapped[b],
                                 Order=str(int(bond.GetBondTypeAsDouble())))
        for b in parent.findall('b'):
            parent.remove(b); parent.append(b)
        expansions[source] = dict(label=raw,normalized=cleaned,generated_atoms=created,
                                  variables=variables,operation='condensed_label')
    return generated, expansions, normalized


def restore_condensed_provenance(mol, generated, expansions, normalized):
    if mol is None:
        return mol
    present = set()
    for a in mol.GetAtoms():
        if not a.HasProp('_CDX_ATOM_ID'):
            continue
        identifier = a.GetProp('_CDX_ATOM_ID')
        if identifier in generated:
            a.SetProp('_CDX_GROUP_SOURCE_ATOM',generated[identifier])
            a.ClearProp('_CDX_ATOM_ID')
        if identifier in expansions:
            present.add(identifier)
            a.SetProp('_CDX_CONDENSED_LABEL',expansions[identifier]['normalized'])
        if identifier in normalized:
            a.SetProp('_CDX_NORMALIZED_LABEL',normalized[identifier])
            if a.GetAtomicNum() == 0:
                a.SetProp('dummyLabel',normalized[identifier])
    if present:
        mol.SetProp('_CDX_LABEL_EXPANSIONS',json.dumps({s:expansions[s] for s in sorted(present)},sort_keys=True))
    return mol
