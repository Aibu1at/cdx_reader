"""
cdx-reader: A Python library for reading ChemDraw CDX/CDXML files
into RDKit molecule objects.
"""

from .cdx_reader import (
    read_cdx_molecule,
    parse_cdxml,
    expand_fragments_in_mol,
    join_molecule_fragments,
    get_fragment_smiles,
    add_first_attachment_point_number,
    FRAGMENT_LIBRARY,
    SCAFFOLD_LABELS,
)

__version__ = "0.1.0"
__all__ = [
    "read_cdx_molecule",
    "parse_cdxml",
    "expand_fragments_in_mol",
    "join_molecule_fragments",
    "get_fragment_smiles",
    "add_first_attachment_point_number",
    "FRAGMENT_LIBRARY",
    "SCAFFOLD_LABELS",
]
