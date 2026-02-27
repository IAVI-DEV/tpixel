"""PNGS (N-linked glycosylation site) detection for protein alignments.

Detects NXS/T motifs (X != P) and maps them to alignment columns.
"""

from __future__ import annotations

from tpixel.hxb2 import HxB2Position
from tpixel.models import Marker


def find_pngs_in_sequence(ungapped_seq: str) -> list[int]:
    """Find PNGS motifs (N-X-S/T, X!=P) in an ungapped protein sequence.

    Returns 0-based positions of the N residue.
    """
    sites: list[int] = []
    for i in range(len(ungapped_seq) - 2):
        if ungapped_seq[i] != "N":
            continue
        if ungapped_seq[i + 1] != "P" and ungapped_seq[i + 2] in ("S", "T"):
            sites.append(i)
    return sites


def _ungapped_to_alignment_map(aligned_seq: str) -> list[int]:
    """Map ungapped positions to alignment column indices."""
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
