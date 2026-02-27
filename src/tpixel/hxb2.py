"""HxB2 coordinate mapping for HIV Env gp160 protein alignments.

Maps alignment columns to HxB2 amino acid positions and Env structural
regions using the LANL convention boundaries.
"""

from __future__ import annotations

from dataclasses import dataclass

from tpixel.models import Region

ENV_REGIONS: list[tuple[str, int, int]] = [
    ("SP",    1,  30),
    ("C1",   31, 130),
    ("V1",  131, 157),
    ("V2",  158, 196),
    ("C2",  197, 295),
    ("V3",  296, 331),
    ("C3",  332, 384),
    ("V4",  385, 418),
    ("C4",  419, 459),
    ("V5",  460, 469),
    ("C5",  470, 511),
    ("gp41", 512, 856),
]

REGION_COLORS: dict[str, str] = {
    "SP":   "#FFF9C4",
    "C1":   "#EEEEEE",
    "V1":   "#BBDEFB",
    "V2":   "#BBDEFB",
    "C2":   "#EEEEEE",
    "V3":   "#BBDEFB",
    "C3":   "#EEEEEE",
    "V4":   "#BBDEFB",
    "C4":   "#EEEEEE",
    "V5":   "#BBDEFB",
    "C5":   "#EEEEEE",
    "gp41": "#F8BBD0",
}

_REGION_LOOKUP: dict[int, str] = {}
for _name, _start, _end in ENV_REGIONS:
    for _pos in range(_start, _end + 1):
        _REGION_LOOKUP[_pos] = _name


def get_env_region(hxb2_aa_pos: int) -> str | None:
    return _REGION_LOOKUP.get(hxb2_aa_pos)


@dataclass
class HxB2Position:
    alignment_col: int
    hxb2_aa_pos: int | None
    region: str | None
    hxb2_residue: str


def build_hxb2_map(aligned_seqs: list[tuple[str, str]], hxb2_id: str = "HxB2") -> list[HxB2Position]:
    """Walk the HxB2 row and map every alignment column to HxB2 coordinates.

    Args:
        aligned_seqs: List of (name, sequence) from read_fasta.
        hxb2_id: Sequence ID of HxB2 in the alignment.

    Returns:
        One HxB2Position per alignment column.
    """
    hxb2_seq = None
    for name, seq in aligned_seqs:
        if name == hxb2_id or name.split()[0] == hxb2_id:
            hxb2_seq = seq
            break

    if hxb2_seq is None:
        raise ValueError(f"HxB2 sequence '{hxb2_id}' not found in alignment")

    positions: list[HxB2Position] = []
    aa_counter = 0

    for col_idx, residue in enumerate(hxb2_seq):
        if residue == "-":
            positions.append(HxB2Position(col_idx, None, None, "-"))
        else:
            aa_counter += 1
            positions.append(HxB2Position(col_idx, aa_counter, get_env_region(aa_counter), residue))

    return positions


def hxb2_col_labels(hxb2_map: list[HxB2Position], step: int = 50) -> list[tuple[int, str]]:
    """Build x-axis tick labels at regular HxB2 AA intervals."""
    max_pos = max((p.hxb2_aa_pos for p in hxb2_map if p.hxb2_aa_pos is not None), default=0)
    labels: list[tuple[int, str]] = []
    for target in range(step, max_pos + 1, step):
        for p in hxb2_map:
            if p.hxb2_aa_pos == target:
                labels.append((p.alignment_col, str(target)))
                break
    return labels


def hxb2_regions(hxb2_map: list[HxB2Position]) -> list[Region]:
    """Build Region annotations from HxB2 position map."""
    region_spans: list[tuple[str, int, int]] = []
    current: str | None = None
    span_start = 0

    for p in hxb2_map:
        r = p.region
        if r != current:
            if current is not None:
                region_spans.append((current, span_start, p.alignment_col))
            current = r
            span_start = p.alignment_col

    if current is not None:
        region_spans.append((current, span_start, len(hxb2_map)))

    return [
        Region(name, start, end, REGION_COLORS.get(name, "#EEEEEE"))
        for name, start, end in region_spans
        if name is not None
    ]
