import re
import logging
import functools
from pycdxml import cdxml_converter

from rdkit import Chem
from rdkit.Chem import AllChem, CombineMols
from rdkit.Chem.rdchem import BondType as BT, Mol, RWMol, BondType
from rdkit import RDLogger

import sys
from io import StringIO, BytesIO
from pathlib import Path
from ._normalization import (
    R_LABEL, atom_label, r_label, capture_r_labels, restore_r_labels,
    normalize_nested_fragment, repair_invalid_nitrogen_hydrogens,
    explicit_element_is_usable,
    has_ambiguous_multimolecule_expansion,
    variable_label, capture_atom_labels, restore_atom_labels, detach_variable_atom_maps,
    sanitize_aromatic_variable_graph,
)

RDLogger.DisableLog('rdApp.*')
logger = logging.getLogger(__name__)

# Disable logging for  pycdxml (gives  ~12% acceleration)
logging.getLogger('pycdxml').setLevel(logging.CRITICAL)
logging.getLogger('pycdxml.chemdraw_objects').setLevel(logging.CRITICAL)


def fix_valence(mol):
    """Attempt to repair valence issues on an RDKit molecule in-place.

    Returns the sanitized molecule, or None if it cannot be fixed.
    """
    from rdkit.Chem import GetPeriodicTable

    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass

    try:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL)
        # No fixes needed
        return mol
    except Exception as exc:
        correcting_atomic_valence = isinstance(exc, Chem.AtomValenceException)

    repairing_with_variables = correcting_atomic_valence and any(
        atom.GetAtomicNum() == 0 for atom in mol.GetAtoms())
    pt = GetPeriodicTable()
    for atom in mol.GetAtoms():
        if repairing_with_variables and atom.GetAtomicNum() == 0:
            # A variable has no default elemental valence. RDKit reports -1;
            # using it as a hydrogen target can abort an unrelated N repair.
            continue
        # Set all explicit Hs to 0 and let RDKit figure it out
        diff_valence = pt.GetDefaultValence(atom.GetAtomicNum()) - atom.GetTotalValence()
        if diff_valence != 0:
            if repairing_with_variables and not (
                atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 0
                and not atom.HasQuery() and not atom.GetNumRadicalElectrons()
                and atom.GetDegree() == 2 and atom.GetNumExplicitHs() == 2
                and diff_valence == -1
                and all(b.GetBondType() == BT.SINGLE for b in atom.GetBonds())
            ):
                # Permit the already-supported NH2 -> NH correction beside a
                # newly recognized variable. Other formerly failed drawings
                # need separate evidence before enabling their H heuristics.
                return None
            numh = atom.GetNumExplicitHs()
            try:
                atom.SetNumExplicitHs(numh + diff_valence)
            except Exception:
                logger.error("can't fix valence error")
                return None
    try:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SANITIZE_ALL)
        # No fixes needed
        return mol
    except Exception:
        return None


SYMBOL_TO_ATOMIC_NUM = {
    "1": "H", "2": "He", "3": "Li", "4": "Be", "5": "B", "6": "C", "7": "N", "8": "O",
    "9": "F", "10": "Ne", "11": "Na", "12": "Mg", "13": "Al", "14": "Si", "15": "P",
    "16": "S", "17": "Cl", "19": "K", "20": "Ca", "35": "Br", "53": "I",  "0":"*",
    "100": "None"
}

ATOMIC_NUM_TO_SYMBOL = { v:k for k, v in SYMBOL_TO_ATOMIC_NUM.items()}

# Fragment library: maps text labels to (SMILES fragment, SMILES with [*:1] attachment point)
# The attachment point [*:1] marks where the fragment connects to the main molecule
FRAGMENT_LIBRARY = {
    # Simple groups - these expand to actual atoms with proper connectivity
    'Me': ('C', '[*:1]C'),                    # Methyl
    'Et': ('CC', '[*:1]CC'),                  # Ethyl
    'Pr': ('CCC', '[*:1]CCC'),                # Propyl
    'iPr': ('C(C)C', 'C([*:1])(C)C'),         # Isopropyl (attach at middle carbon)
    'i-Pr': ('C(C)C', 'C([*:1])(C)C'),        # Isopropyl variant
    'nPr': ('CCC', '[*:1]CCC'),               # n-Propyl
    'n-Pr': ('CCC', '[*:1]CCC'),              # n-Propyl variant
    'Bu': ('CCCC', '[*:1]CCCC'),              # Butyl
    'nBu': ('CCCC', '[*:1]CCCC'),             # n-Butyl
    'n-Bu': ('CCCC', '[*:1]CCCC'),            # n-Butyl variant
    'iBu': ('CC(C)C', '[*:1]CC(C)C'),         # Isobutyl
    'i-Bu': ('CC(C)C', '[*:1]CC(C)C'),        # Isobutyl variant
    'tBu': ('C(C)(C)C', '[*:1]C(C)(C)C'),     # tert-Butyl
    't-Bu': ('C(C)(C)C', '[*:1]C(C)(C)C'),    # tert-Butyl variant
    'sBu': ('C(C)CC', 'C([*:1])(C)CC'),       # sec-Butyl (attach at middle carbon)
    's-Bu': ('C(C)CC', 'C([*:1])(C)CC'),      # sec-Butyl variant

    # Aromatic groups
    'Ph': ('c1ccccc1', '[*:1]c1ccccc1'),              # Phenyl
    'Bn': ('Cc1ccccc1', '[*:1]Cc1ccccc1'),            # Benzyl
    'Bz': ('C(=O)c1ccccc1', '[*:1]C(=O)c1ccccc1'),    # Benzoyl

    # Oxygen-containing groups
    'MeO': ('CO', 'CO[*:1]'),                  # Methoxy: C-*
    'OMe': ('OC', '[*:1]OC'),                 # Methoxy: *-O-C
    'PhO': ('Oc1ccccc1', '[*:1]Oc1ccccc1'),   # Phenoxy
    'EtO': ('CCO', 'CCO[*:1]'),                # Ethoxy: CC-*
    'OEt': ('OCC', '[*:1]OCC'),               # Ethoxy variant
    'Ac': ('C(=O)C', '[*:1]C(=O)C'),          # Acetyl
    'AcO': ('OC(=O)C', '[*:1]OC(=O)C'),       # Acetoxy
    'OAc': ('OC(=O)C', '[*:1]OC(=O)C'),       # Acetoxy variant
    'OCF3': ('OC(F)(F)F', '[*:1]OC(F)(F)F'),  # Trifluoromethoxy

    # Nitrogen-containing groups
    'MeN': ('NC', 'CN[*:1]'),                 # Methylamino (attach at N)
    'NMe': ('NC', '[*:1]NC'),                 # N-methyl
    'NMe2': ('N(C)C', '[*:1]N(C)C'),          # Dimethylamino
    'NEt2': ('N(CC)CC', '[*:1]N(CC)CC'),      # Diethylamino
    'NO2': ('[N+](=O)[O-]', '[*:1][N+](=O)[O-]'),  # Nitro
    'CN': ('C#N', '[*:1]C#N'),                # Cyano
    'NC': ('C#N', '[*:1]C#N'),                # Cyano variant
    'HCl H2N': ('N', '[*:1]N'),               # Amine hydrochloride (just the amine part)

    # Halogenated groups
    'CF3': ('C(F)(F)F', '[*:1]C(F)(F)F'),           # Trifluoromethyl
    'F3C': ('C(F)(F)F', '[*:1]C(F)(F)F'),           # Trifluoromethyl variant
    'F3CO': ('OC(F)(F)F', '[*:1]OC(F)(F)F'),        # Trifluoromethoxy variant
    'CCl3': ('C(Cl)(Cl)Cl', '[*:1]C(Cl)(Cl)Cl'),    # Trichloromethyl
    'CHF2': ('C(F)F', '[*:1]C(F)F'),                # Difluoromethyl
    'CH2F': ('CF', '[*:1]CF'),                      # Fluoromethyl

    # Sulfur-containing groups
    'SMe': ('SC', '[*:1]SC'),                       # Methylthio
    'MeS': ('CS', 'C[*:1]'),                        # Methylthio variant (attach at S)
    'SO2Me': ('S(=O)(=O)C', '[*:1]S(=O)(=O)C'),     # Methanesulfonyl
    'MeSO2': ('CS(=O)(=O)', '[*:1]CS(=O)(=O)'),     # Methanesulfonyl variant
    'SO2NH2': ('S(=O)(=O)N', '[*:1]S(=O)(=O)N'),    # Sulfonamide
    'SO3H': ('S(=O)(=O)O', '[*:1]S(=O)(=O)O'),      # Sulfonic acid
    # 'SH': ('S', '[*:1]S'),                          # Thiol

    # Boron-containing groups
    'B(OH)2': ('B(O)O', '[*:1]B(O)O'),                # Boronic acid

    # Carbonyl groups
    'CHO': ('C=O', '[*:1]C=O'),                     # Formyl/Aldehyde
    'CO2H': ('C(=O)O', '[*:1]C(=O)O'),              # Carboxylic acid
    'COOH': ('C(=O)O', '[*:1]C(=O)O'),              # Carboxylic acid variant
    'CO2Me': ('C(=O)OC', '[*:1]C(=O)OC'),           # Methyl ester
    'CO2Et': ('C(=O)OCC', '[*:1]C(=O)OCC'),         # Ethyl ester
    'CONH2': ('C(=O)N', '[*:1]C(=O)N'),             # Primary amide
    'CONMe2': ('C(=O)N(C)C', '[*:1]C(=O)N(C)C'),    # Dimethyl amide

    # Common protecting groups
    'Boc': ('C(=O)OC(C)(C)C', '[*:1]C(=O)OC(C)(C)C'),        # tert-Butoxycarbonyl
    'Cbz': ('C(=O)OCc1ccccc1', '[*:1]C(=O)OCc1ccccc1'),      # Carboxybenzyl
    'Fmoc': ('C(=O)OCC1c2ccccc2-c2ccccc12', '[*:1]C(=O)OCC1c2ccccc2-c2ccccc12'),  # Fluorenylmethyloxycarbonyl

    # Complex protecting group abbreviations (ChemDraw style)
    # BnO2CHN = Bn-O-CO-NH- = benzyloxycarbonylamino (Cbz-NH-)
    'BnO2C': ('C(=O)OCc1ccccc1', '[*:1]C(=O)OCc1ccccc1'),  # Benzyl methanoate
    'BnO2CHN': ('NC(=O)OCc1ccccc1', '[*:1]NC(=O)OCc1ccccc1'),  # Benzyloxycarbonylamino
    'CbzHN': ('NC(=O)OCc1ccccc1', '[*:1]NC(=O)OCc1ccccc1'),    # Cbz-amino variant
    'BocHN': ('NC(=O)OC(C)(C)C', '[*:1]NC(=O)OC(C)(C)C'),      # Boc-amino
    'FmocHN': ('NC(=O)OCC1c2ccccc2-c2ccccc12', '[*:1]NC(=O)OCC1c2ccccc2-c2ccccc12'),  # Fmoc-amino
    'AcHN': ('NC(=O)C', '[*:1]NC(=O)C'),                       # Acetamido
    'TfaHN': ('NC(=O)C(F)(F)F', '[*:1]NC(=O)C(F)(F)F'),        # Trifluoroacetamido

    # More carbamate/amide abbreviations
    'NHBoc': ('NC(=O)OC(C)(C)C', '[*:1]NC(=O)OC(C)(C)C'),      # Boc-amino (reversed notation)
    'NHCbz': ('NC(=O)OCc1ccccc1', '[*:1]NC(=O)OCc1ccccc1'),    # Cbz-amino (reversed notation)
    'NHAc': ('NC(=O)C', '[*:1]NC(=O)C'),                       # Acetamido (reversed notation)
    'NHBn': ('NCc1ccccc1', '[*:1]NCc1ccccc1'),                 # Benzylamino
    'BnHN': ('NCc1ccccc1', '[*:1]NCc1ccccc1'),                 # Benzylamino variant

    # Sulfonamide abbreviations
    'TsHN': ('NS(=O)(=O)c1ccc(C)cc1', '[*:1]NS(=O)(=O)c1ccc(C)cc1'),  # Tosylamino
    'MsHN': ('NS(=O)(=O)C', '[*:1]NS(=O)(=O)C'),                      # Mesylamino
    'NHTs': ('NS(=O)(=O)c1ccc(C)cc1', '[*:1]NS(=O)(=O)c1ccc(C)cc1'),  # Tosylamino (reversed)
    'NHMs': ('NS(=O)(=O)C', '[*:1]NS(=O)(=O)C'),                      # Mesylamino (reversed)
}

# Labels that indicate scaffold/variable positions (should result in is_scaffold=True)
SCAFFOLD_LABELS = {
    'R', 'R1', 'R2', 'R3', 'R4', 'R5', 'R6', 'R7', 'R8', 'R9', 'R10',
    'Ra', 'Rb', 'Rc', 'Rd', 'Re', 'Rf', 'Rg',
    'X', 'Y', 'Z', 'A', 'Q', 'W',
    '*', '(R2)0-3', '(R2)0-1',
    'C1-4 alkyl', 'C1-4alkyl', 'C1-4 alkoxy', 'C1-4alkoxy',
    'C1-6 alkyl', 'C1-6alkyl',
}

ATTACHMENT_POINT_TOKEN = "*"
ATTACHMENT_POINT_REGEXP = r"(?:{0}|\[{0}[^\]]*\])".format(re.escape(ATTACHMENT_POINT_TOKEN))


def get_fragment_smiles(label: str) -> tuple:
    """
    Get SMILES fragment and attachment info for a label.

    Returns:
        (smiles, smiles_with_attachment, is_scaffold) or (None, None, False) if not found
    """
    # Direct lookup
    if label in FRAGMENT_LIBRARY:
        smiles, smiles_with_attachment = FRAGMENT_LIBRARY[label]
        return smiles, smiles_with_attachment, False

    # Check scaffold labels
    if label in SCAFFOLD_LABELS:
        return None, None, True

    # Check for scaffold patterns (Rn, etc.)
    # Keep the legacy isdigit coverage (including literal superscript digits)
    # as well as styled/indexed letter variables.
    if R_LABEL.fullmatch(label) or (label.startswith('R') and label[1:].isdigit()):
        return None, None, True

    # Labels with commas, semicolons, or "and" are usually lists of options (scaffold)
    if ',' in label or ';' in label or ' and ' in label.lower() or 'or' in label.lower():
        return None, None, True

    # Check for alkyl/alkoxy patterns with ranges (scaffold indicators)
    if 'alkyl' in label.lower() or 'alkoxy' in label.lower():
        return None, None, True

    return None, None, False


def parse_cdxml(tree, fix=True):
    """
    Parse and fix CDXML molecule.

    Handles both simple element labels (N, O, S) and complex fragment labels
    (Me, MeO, CF3, Ph, etc.) by mapping them to appropriate atoms.

    For fragment labels, sets the atom to the attachment point element
    (e.g., MeO -> O, Me -> C, NH2 -> N).

    Returns:
        Tuple of:
        - List of (node_order, label, node_id) tuples for multi-atom fragments that need expansion.
          node_order is the 0-based index of the node in the order nodes appear in CDXML.
        - Dict mapping node_id -> node_order for all nodes in the CDXML
    """
    return _parse_cdxml(tree, fix=fix)


def _parse_cdxml(tree, fix=True, repair_nested=True):
    root = tree.getroot()
    preserve_selection = has_ambiguous_multimolecule_expansion(root, FRAGMENT_LIBRARY)
    fragments_to_expand = []  # List of (node_order, label, node_id) for multi-atom fragments
    node_id_to_order = {}  # Maps node_id -> order index

    # Standard single-element labels (case sensitive)
    SINGLE_ELEMENTS = {'H', 'C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I', 'Si', 'B'}

    # ChemDraw XML uses 'n' tag for nodes (atoms)
    # First pass: build node_id -> order mapping
    for node_order, node in enumerate(root.iter('n')):
        node_id = node.get('id')
        if node_id:
            node_id_to_order[node_id] = node_order

    # Second pass: process nodes and identify fragments to expand
    for node_order, node in enumerate(root.iter('n')):
        # Get the node's id attribute - this is the reliable identifier
        node_id = node.get('id')

        if fix:
            repair_invalid_nitrogen_hydrogens(node)
            # Prefer explicit internal atoms/bonds over interpreting the outer
            # text. Some ChemDraw files incorrectly mark these as Unspecified.
            if repair_nested and not preserve_selection and normalize_nested_fragment(node):
                continue

        # 'Element' is often implicit (Carbon) or in 'Element' attribute
        # 't' tag inside 'n' often contains labels like "Me", "Ph"
        element = node.get('Element', '100')  # Default to Carbon if missing

        # Check for text label
        label = atom_label(node).replace('\n', '').replace('\r', '')
        canonical_r_label = r_label(node)
        if canonical_r_label:
            label = canonical_r_label

        if not label:
            continue

        declared_variable = variable_label(node, FRAGMENT_LIBRARY)
        if fix and declared_variable:
            node.set('Element', '0')
            continue

        # Skip if element is already correctly set
        is_bad_element = element == '100'
        if not is_bad_element and not fix:
            continue

        # Explicit element identity is stronger evidence than a prefix such as
        # "N" in "Na+". R placeholders and known abbreviations still need the
        # existing interpretation/expansion path. An Unspecified node with no
        # Element (e.g. HN) continues through the label repair below.
        if (node.get('Element') is not None and element not in ('0', '99', '100')
                and explicit_element_is_usable(node)
                and label not in FRAGMENT_LIBRARY
                and label not in SCAFFOLD_LABELS and not get_fragment_smiles(label)[2]):
            continue

        # Determine what element this label should map to
        target_element = None

        # Check if it's a simple element label
        if label in SINGLE_ELEMENTS:
            target_element = label
        # Check fragment library
        elif label in FRAGMENT_LIBRARY:
            frag_smiles, frag_smiles_with_attachment = FRAGMENT_LIBRARY[label]
            # Parse the attached SMILES to get the attachment atom element
            # The atom connected to [*:1] is the attachment point
            fragments_to_expand.append((node_order, label, node_id))

        # Check if it's a scaffold label (R groups, etc.)
        elif label in SCAFFOLD_LABELS or get_fragment_smiles(label)[2]:
            # This is a scaffold/variable position - set to dummy atom
            target_element = '*'
        # Handle common patterns not in library
        else:
            # Check for nitrogen-containing labels
            if label.startswith('NH') or label.startswith('HN') or label.startswith('H2N') or label.startswith('N'):
                target_element = 'N'
            # Check for oxygen-containing labels
            elif label.startswith('HO') or label.startswith('OH') or label.endswith('O') and label not in {'NO2'}:
                target_element = 'O'
            # Check for sulfur-containing labels
            elif label.startswith('S') or label.endswith('S'):
                target_element = 'S'
            # Default to carbon for unknown organic fragments
            elif len(label) >=2:
                logger.error('Invalid fragment label "%s"' % label)
            elif any(c.isalpha() for c in label):
                target_element = 'C'

        # Update the element if we determined a target
        if target_element and fix:
            node.set('Element', ATOMIC_NUM_TO_SYMBOL.get(target_element, '0'))
        # print(label, ATOMIC_NUM_TO_SYMBOL.get(target_element, '0'), target_element)
    return fragments_to_expand, node_id_to_order


def split_mol_at_atom(mol, atom_idx: int):
    """
    Split molecule at given atom, creating fragments with attachment points.

    For each bond from atom_idx to its neighbors, creates a fragment containing
    that neighbor's connected component with an attachment point (*) where the
    bond was cut.

    Args:
        mol: RDKit molecule
        atom_idx: Index of atom to split at

    Returns:
        List of (fragment_mol, attachment_atom_idx) tuples
    """

    # Try to find a bond between index atom and its' nieghbor
    neighbor_idx = atom_idx - 1
    bond = mol.GetBondBetweenAtoms(neighbor_idx, atom_idx) if neighbor_idx >= 0 else None

    if bond is None:
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbor_idx = min([x.GetIdx() for x in atom.GetNeighbors()])
        bond = mol.GetBondBetweenAtoms(neighbor_idx, atom_idx)
    mol_frag = Chem.FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)
    mol_frag = Chem.GetMolFrags(mol_frag, asMols=True)
    return mol_frag


def fragment_matches_library(frag_mol, expected_smiles: str) -> bool:
    """
    Check if a fragment matches the expected SMILES from FRAGMENT_LIBRARY.

    Compares by canonical SMILES after removing attachment points.

    Args:
        frag_mol: Fragment molecule (may contain attachment point)
        expected_smiles: Expected SMILES from FRAGMENT_LIBRARY

    Returns:
        True if fragment matches expected structure
    """
    if frag_mol is None:
        return False

    try:
        # Remove attachment points from fragment for comparison
        frag_emol = RWMol(frag_mol)
        atoms_to_remove = []
        for atom in frag_emol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                atoms_to_remove.append(atom.GetIdx())

        # Remove in reverse order to preserve indices
        for idx in sorted(atoms_to_remove, reverse=True):
            frag_emol.RemoveAtom(idx)

        frag_clean = frag_emol.GetMol()

        # Get expected molecule
        expected_mol = Chem.MolFromSmiles(expected_smiles)
        if expected_mol is None:
            return False

        # Compare by InChI or canonical SMILES (more robust than direct comparison)
        try:
            frag_inchi = Chem.MolToInchi(frag_clean)
            expected_inchi = Chem.MolToInchi(expected_mol)
            if frag_inchi and expected_inchi and frag_inchi == expected_inchi:
                return True
        except:
            pass

        # Fallback: compare canonical SMILES
        try:
            Chem.SanitizeMol(frag_clean)
            frag_can = Chem.MolToSmiles(frag_clean)
            expected_can = Chem.MolToSmiles(expected_mol)
            return frag_can == expected_can
        except:
            pass

        # Last resort: compare heavy atom count and formula
        if frag_clean.GetNumHeavyAtoms() == expected_mol.GetNumHeavyAtoms():
            frag_atoms = sorted([a.GetSymbol() for a in frag_clean.GetAtoms()])
            exp_atoms = sorted([a.GetSymbol() for a in expected_mol.GetAtoms()])
            return frag_atoms == exp_atoms

        return False

    except Exception as e:
        logger.debug(f"Error comparing fragments: {e}")
        return False


def find_atom_by_cdx_id(mol, cdx_atom_id: str) -> int:
    """
    Find atom index by its CDX atom ID property.

    RDKit preserves the CDXML node 'id' attribute as '_CDX_ATOM_ID' property.

    Args:
        mol: RDKit molecule
        cdx_atom_id: The CDXML node id to find

    Returns:
        Atom index if found, -1 otherwise
    """
    target_id = int(cdx_atom_id) if cdx_atom_id else -1
    for atom in mol.GetAtoms():
        if atom.HasProp("_CDX_ATOM_ID"):
            try:
                if int(atom.GetProp("_CDX_ATOM_ID")) == target_id:
                    return atom.GetIdx()
            except ValueError:
                pass
    return -1


def expand_fragments_in_mol(mol, fragments_to_expand):
    """
    Expand multi-atom fragment abbreviations in a molecule.

    Uses split_mol_at_atom to split the molecule and check if one of the fragments
    matches the expected fragment from FRAGMENT_LIBRARY:
    1. For each fragment label, find the atom by its CDX node ID (_CDX_ATOM_ID property)
    2. Split the molecule at that atom
    3. Check if any resulting fragment matches the expected structure (smallest first)
    4. If match found -> fragment already expanded, keep original molecule
    5. If no match -> join the scaffold and expected fragment using join_molecule_fragments

    Args:
        mol: RDKit Mol object
        fragments_to_expand: Tuple of (fragments_list, node_id_to_order) from parse_cdxml
                            where fragments_list contains (node_order, label, node_id) tuples

    Returns:
        New molecule with fragments expanded
    """
    if not fragments_to_expand or mol is None:
        return mol

    # Handle the tuple return from parse_cdxml
    if isinstance(fragments_to_expand, tuple):
        fragments_list, node_id_to_order = fragments_to_expand
    else:
        fragments_list = fragments_to_expand
        node_id_to_order = {}

    if not fragments_list:
        return mol

    # Process fragments in reverse order to avoid index shifting issues
    for fragment_info in sorted(fragments_list, key=lambda x: -x[0]):
        node_order, label, node_id = fragment_info

        # Find the atom by its CDX ID (node_id from CDXML becomes _CDX_ATOM_ID in RDKit)
        atom_idx = find_atom_by_cdx_id(mol, node_id)

        if atom_idx < 0:
            logger.debug(f"Could not find atom with CDX ID {node_id} for {label}")
            continue

        if label not in FRAGMENT_LIBRARY:
            logger.debug(f"Fragment {label} not in FRAGMENT_LIBRARY")
            continue

        frag_smiles, frag_smiles_with_attachment = FRAGMENT_LIBRARY[label]
        frag_mol = Chem.MolFromSmiles(frag_smiles)

        if frag_mol is None or frag_mol.GetNumAtoms() <= 1:
            continue

        try:
            # Split molecule at the labeled atom using FragmentOnBonds
            split_frags = split_mol_at_atom(mol, atom_idx)

            if not split_frags or len(split_frags) < 2:
                logger.debug(f"Fragment {label} at atom {atom_idx}: could not split molecule")
                continue

            # Sort fragments by size (smallest first) to check the likely fragment first
            split_frags_sorted = sorted(split_frags, key=lambda m: m.GetNumAtoms())

            # Check if any fragment matches the expected structure from FRAGMENT_LIBRARY
            fragment_already_present = False
            scaffold_mol = None

            for i, frag in enumerate(split_frags_sorted):
                if fragment_matches_library(frag, frag_smiles):
                    logger.debug(f"Fragment {label} at atom {atom_idx} (cdx_id={node_id}) found in split fragments, skipping expansion")
                    fragment_already_present = True
                    break
                else:
                    # Keep track of the scaffold (the larger fragment that's not the expected one)
                    if scaffold_mol is None or frag.GetNumAtoms() > scaffold_mol.GetNumAtoms():
                        scaffold_mol = frag

            if fragment_already_present:
                continue

            # Fragment not present - need to expand using join_molecule_fragments
            logger.debug(f"Expanding fragment {label} at atom {atom_idx} (cdx_id={node_id})")

            # Use the larger fragment as scaffold (it has the attachment point from split)
            # The smaller one is likely the placeholder atom that needs replacement
            if scaffold_mol is None:
                scaffold_mol = split_frags_sorted[-1]  # Use largest fragment as scaffold

            # Only FragmentOnBonds' new cut atoms belong to the joining
            # namespace. Source R atoms keep their labels/isotopes, and a cut
            # next to atom zero is still a valid attachment (isotope zero).
            scaffold_rw = RWMol(scaffold_mol)
            for atom in scaffold_rw.GetAtoms():
                if atom.GetAtomicNum() == 0 and not atom.HasProp('_CDX_ATOM_ID'):
                    atom.SetAtomMapNum(1)
                    atom.SetIsotope(0)
            scaffold_mol = scaffold_rw.GetMol()

            # Use the pre-defined attached SMILES from FRAGMENT_LIBRARY
            # It already has [*:1] attachment point, just need to parse it
            decoration_mol = Chem.MolFromSmiles(frag_smiles_with_attachment)

            if decoration_mol is None:
                # Fallback: try constructing it manually
                decoration_smi = add_first_attachment_point_number(f"[*]{frag_smiles}", 1)
                decoration_mol = Chem.MolFromSmiles(decoration_smi)

            if decoration_mol is None:
                logger.debug(f"Could not create decoration molecule for {label}")
                continue

            # Join scaffold and decoration (regenerate 2D coords for clean layout)
            joined_mol = join_molecule_fragments(
                scaffold_mol,
                decoration_mol,
                regenerate_2d_coords=True
            )

            if joined_mol is not None:
                mol = joined_mol
                logger.debug(f"Successfully expanded fragment {label} at atom {atom_idx} (cdx_id={node_id})")
            else:
                logger.debug(f"join_molecule_fragments failed for {label} at atom {atom_idx}")

        except Exception as e:
            logger.debug(f"Failed to expand fragment {label} at atom {atom_idx}: {e}")
            continue

    return mol

def read_cdx_molecule(path_or_bytes, as_single=True):
    """Read a CDX file and return one or more RDKit molecules.

    Fixes bad atom types from ChemDraw fragment labels and expands multi-atom
    abbreviations (Me, Ph, Boc, CF3, ...) into their full structures.

    Args:
        path_or_bytes: Path (str / Path) to a CDX file, or raw CDX bytes.
        as_single: If True (default), combine all fragments into a single molecule.
            If False, return a list of molecules.
    """
    return _read_cdx_molecule(path_or_bytes, as_single, repair_nested=True)


def _read_cdx_molecule(path_or_bytes, as_single, repair_nested):
    # Redirect stderr to suppress font warnings from pycdxml
    old_stderr = sys.stderr
    sys.stderr = StringIO()

    try:
        if isinstance(path_or_bytes, bytes):
            source = BytesIO(path_or_bytes)
        elif isinstance(path_or_bytes, (str, Path)):
            source = str(path_or_bytes)
        else:
            raise ValueError(f"Invalid type {type(path_or_bytes)}")

        doc: cdxml_converter.ChemDrawDocument = cdxml_converter.read_cdx(source)
        labels = capture_r_labels(doc.cdxml)
        all_labels = capture_atom_labels(doc.cdxml)
        original_types = {n.get('id'): n.get('NodeType') for n in doc.cdxml.getroot().iter('n')}
        fragments_to_expand = _parse_cdxml(doc.cdxml, fix=True, repair_nested=repair_nested)
        repaired_fragment = any(n.get('NodeType') == 'Fragment'
                                and original_types.get(n.get('id')) != 'Fragment'
                                for n in doc.cdxml.getroot().iter('n'))
        mols = Chem.rdmolfiles.MolsFromCDXML(doc.to_cdxml(), sanitize=False)
        for mol in mols:
            detach_variable_atom_maps(mol, all_labels)

        def finalize(mol):
            original = Chem.Mol(mol)
            mol = fix_valence(mol)
            if mol is None:
                mol = sanitize_aromatic_variable_graph(original, all_labels)
            if mol is not None and repaired_fragment:
                # Expanding CDXML nodes changes atom indices. Recompute cached
                # stereo from retained bond directions before serializing.
                # Cleaning treats R3/R4 as identical wildcard substituents and
                # can erase their drawn tetrahedral orientation. Preserve the
                # imported tags for the downstream agent choosing those groups.
                drawn_tags = [(atom, atom.GetChiralTag()) for atom in mol.GetAtoms()
                              if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]
                Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
                for atom, tag in drawn_tags:
                    atom.SetChiralTag(tag)
            return restore_r_labels(restore_atom_labels(mol, all_labels), labels)

        if as_single:
            combined = functools.reduce(Chem.CombineMols, mols)
            if fragments_to_expand:
                combined = expand_fragments_in_mol(combined, fragments_to_expand)
            result = finalize(combined)
            if result is None and repaired_fragment:
                # A valid nested group can expose an unrelated error in the
                # surrounding drawing which the legacy valence fixer cannot
                # repair (e.g. excess H next to a newly restored nitro N+).
                # Preserve resolution without broadening hydrogen repairs.
                return _read_cdx_molecule(path_or_bytes, as_single, repair_nested=False)
            return result

        result = []
        for mol in mols:
            if fragments_to_expand:
                mol = expand_fragments_in_mol(mol, fragments_to_expand)
            result.append(finalize(mol))
        if repaired_fragment and any(mol is None for mol in result):
            return _read_cdx_molecule(path_or_bytes, as_single, repair_nested=False)
        return result
    finally:
        sys.stderr = old_stderr


def add_first_attachment_point_number(smi, num):
    """
    Changes/adds a number to the first attachment point.
    :param smi: SMILES string with the molecule.
    :param num: Number to add.
    :return: A SMILES string with the number added.
    """
    return re.sub(
        ATTACHMENT_POINT_REGEXP,
        "[{}:{}]".format(ATTACHMENT_POINT_TOKEN, num),
        smi,
        count=1,
    )



def join_molecule_fragments(scaffold: Mol, decoration: Mol, regenerate_2d_coords=False):
    """
    Joins a RDKit MOL scaffold with a decoration. They must be labelled.
    :param scaffold: RDKit MOL of the scaffold.
    :param decoration: RDKit MOL of the decoration.
    This is useful when debugging, but it can give problems.
    :param regenerate_2d_coords: If True, regenerate 2D coordinates after joining.
    :return: A Mol object of the joined scaffold.
    """

    if scaffold and decoration:
        # obtain id in the decoration
        try:
            attachment_points = [
                atom.GetProp("molAtomMapNumber")
                for atom in decoration.GetAtoms()
                if atom.GetSymbol() == ATTACHMENT_POINT_TOKEN
            ]
            if len(attachment_points) != 1:
                return None  # more than one attachment point...
            attachment_point = attachment_points[0]
        except KeyError:
            return None

        combined_scaffold = RWMol(CombineMols(decoration, scaffold))
        attachments = [
            atom
            for atom in combined_scaffold.GetAtoms()
            if atom.GetSymbol() == ATTACHMENT_POINT_TOKEN
            and atom.HasProp("molAtomMapNumber")
            and atom.GetProp("molAtomMapNumber") == attachment_point
        ]
        if len(attachments) != 2:
            return None  # something weird

        neighbors = []
        for atom in attachments:
            if atom.GetDegree() != 1:
                return None  # the attachment is wrongly generated
            neighbors.append(atom.GetNeighbors()[0])

        bonds = [atom.GetBonds()[0] for atom in attachments]
        bond_type = BondType.SINGLE
        if any(bond for bond in bonds if bond.GetBondType() == BondType.DOUBLE):
            bond_type = BondType.DOUBLE

        combined_scaffold.AddBond(neighbors[0].GetIdx(), neighbors[1].GetIdx(), bond_type)
        combined_scaffold.RemoveAtom(attachments[0].GetIdx())
        combined_scaffold.RemoveAtom(attachments[1].GetIdx())

        # Label the atoms in the bond
        bondNumbers = [
            int(atom.GetProp("bondNum"))
            for atom in combined_scaffold.GetAtoms()
            if atom.HasProp("bondNum")
        ]

        if bondNumbers:
            bondNum = max(bondNumbers) + 1
        else:
            bondNum = 0

        for neighbor in neighbors:
            idx = neighbor.GetIdx()
            atom = combined_scaffold.GetAtomWithIdx(idx)
            atom.SetIntProp("bondNum", bondNum)

        scaffold = combined_scaffold.GetMol()

        # Regenerate 2D coordinates if requested
        if regenerate_2d_coords and scaffold is not None:
            AllChem.Compute2DCoords(scaffold)

    else:
        return None

    return scaffold
