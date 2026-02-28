"""PNGS (N-linked glycosylation site) detection for protein alignments.

Detects NXS/T motifs (X != P) and maps them to alignment columns.
Supports both amino-acid and nucleotide alignments (NT sequences are
translated internally before motif scanning).
"""

from __future__ import annotations

from tpixel.hxb2 import HxB2Position
from tpixel.models import Marker

# Standard genetic code -------------------------------------------------------

CODON_TABLE: dict[str, str] = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


def translate(nt_seq: str) -> str:
    """Translate an ungapped nucleotide sequence to protein.

    Codons containing ambiguous bases (e.g. ``N``) are rendered as ``X``.
    A trailing incomplete codon is silently ignored.

    Args:
        nt_seq: Ungapped nucleotide sequence string.

    Returns:
        Translated protein sequence.

    Examples:
        >>> translate("ATGAATGCTTCT")
        'MNAS'
        >>> translate("ATGNNN")
        'MX'
        >>> translate("ATG")
        'M'
        >>> translate("ATGA")
        'M'
    """
    protein: list[str] = []
    seq = nt_seq.upper().replace("U", "T")
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i : i + 3]
        protein.append(CODON_TABLE.get(codon, "X"))
    return "".join(protein)


def find_pngs_in_sequence(ungapped_seq: str) -> list[int]:
    """Find PNGS motifs (N-X-S/T, X!=P) in an ungapped protein sequence.

    Args:
        ungapped_seq: Ungapped protein sequence string.

    Returns:
        0-based positions of the N residue in each PNGS motif.

    Examples:
        >>> find_pngs_in_sequence("NAS")
        [0]
        >>> find_pngs_in_sequence("NAT")
        [0]
        >>> find_pngs_in_sequence("NPS")
        []
        >>> find_pngs_in_sequence("MNASNPS")
        [1]
    """
    sites: list[int] = []
    for i in range(len(ungapped_seq) - 2):
        if ungapped_seq[i] != "N":
            continue
        if ungapped_seq[i + 1] != "P" and ungapped_seq[i + 2] in ("S", "T"):
            sites.append(i)
    return sites


def _ungapped_to_alignment_map(aligned_seq: str) -> list[int]:
    """Map ungapped positions to alignment column indices.

    Args:
        aligned_seq: Aligned sequence string (may contain ``-`` gaps).

    Returns:
        List of column indices for non-gap positions.

    Examples:
        >>> _ungapped_to_alignment_map("A-C-G")
        [0, 2, 4]
        >>> _ungapped_to_alignment_map("ACG")
        [0, 1, 2]
        >>> _ungapped_to_alignment_map("---")
        []
    """
    return [col for col, res in enumerate(aligned_seq) if res != "-"]


def find_pngs_markers(
    aligned_ref_seq: str,
    hxb2_map: list[HxB2Position] | None = None,
) -> list[Marker]:
    """Detect PNGS in aligned reference and return Marker annotations.

    Labels use HxB2 positions if hxb2_map is provided, otherwise column indices.
    """
    ungapped = aligned_ref_seq.replace("-", "")
    pngs_positions = find_pngs_in_sequence(ungapped)
    col_map = _ungapped_to_alignment_map(aligned_ref_seq)

    markers: list[Marker] = []
    for pos in pngs_positions:
        col = col_map[pos]
        if hxb2_map and col < len(hxb2_map) and hxb2_map[col].hxb2_aa_pos is not None:
            label = f"N{hxb2_map[col].hxb2_aa_pos}"
        else:
            label = f"N{pos + 1}"
        markers.append(Marker(col=col, label=label))

    return markers


def find_pngs_markers_nt(
    aligned_ref_seq: str,
    hxb2_map: list[HxB2Position] | None = None,
) -> list[Marker]:
    """Detect PNGS from an aligned *nucleotide* reference sequence.

    Workflow:
    1. Strip gaps to get ungapped NT sequence.
    2. Translate to protein via :func:`translate`.
    3. Find NXS/T motifs with :func:`find_pngs_in_sequence`.
    4. Map each AA position back to the alignment column of the first
       nucleotide of the corresponding codon.
    5. Label with HxB2 AA position when available.
    """
    ungapped_nt = aligned_ref_seq.replace("-", "").replace(".", "")
    protein = translate(ungapped_nt)
    pngs_aa_positions = find_pngs_in_sequence(protein)

    # Map: ungapped NT index → alignment column
    nt_col_map = _ungapped_to_alignment_map(aligned_ref_seq)

    markers: list[Marker] = []
    for aa_pos in pngs_aa_positions:
        # First NT of the codon in ungapped coordinates
        nt_ungapped_idx = aa_pos * 3
        if nt_ungapped_idx >= len(nt_col_map):
            continue
        col = nt_col_map[nt_ungapped_idx]

        if hxb2_map and col < len(hxb2_map) and hxb2_map[col].hxb2_aa_pos is not None:
            label = f"N{hxb2_map[col].hxb2_aa_pos}"
        else:
            label = f"N{aa_pos + 1}"
        markers.append(Marker(col=col, label=label))

    return markers
