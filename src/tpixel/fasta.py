"""FASTA parsing and panel construction."""

from __future__ import annotations

from pathlib import Path

from tpixel.models import Panel


def read_fasta(path: str | Path) -> list[tuple[str, str]]:
    """Parse a FASTA file into a list of (name, sequence) tuples.

    Args:
        path: Path to the FASTA file.

    Returns:
        List of (header_name, concatenated_sequence) tuples.
    """
    seqs: list[tuple[str, str]] = []
    name: str | None = None
    buf: list[str] = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name = line[1:].strip()
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def fasta_panel(
    path: str | Path,
    col_start: int | None = None,
    col_end: int | None = None,
    ref_positions: list[int] | None = None,
) -> Panel:
    """Build a Panel from an aligned FASTA.

    Args:
        path: Path to the aligned FASTA file.
        col_start: 1-based inclusive start column for slicing the alignment.
        col_end: 1-based inclusive end column for slicing the alignment.
        ref_positions: 1-based positions of reference sequences. Last is
            the primary reference; earlier ones become extra reference rows.
            Defaults to [1].

    Returns:
        A Panel with reference row, sequence rows, and column labels.

    Raises:
        ValueError: If the FASTA file contains no sequences.
    """
    if ref_positions is None:
        ref_positions = [1]

    seqs = read_fasta(path)
    if not seqs:
        raise ValueError(f"No sequences in {path}")

    # Primary reference is the last position in ref_positions
    primary_idx = ref_positions[-1] - 1
    _ref_name, ref_seq = seqs[primary_idx]

    # Slice columns if requested (1-based inclusive)
    if col_start is not None or col_end is not None:
        cs = (col_start or 1) - 1
        ce = col_end or len(ref_seq)
        ref_seq = ref_seq[cs:ce]
        seqs = [(n, s[cs:ce]) for n, s in seqs]

    aln_len = len(ref_seq)
    ref_row = list(ref_seq.upper())

    # Extra reference rows (all ref positions except the last)
    extra_ref_rows: list[tuple[str, list[str]]] = []
    for pos in ref_positions[:-1]:
        idx = pos - 1
        name, seq = seqs[idx]
        row = list(seq.upper()[:aln_len])
        row += ["-"] * (aln_len - len(row))
        extra_ref_rows.append((name, row))

    # Sample sequences: everything not in ref_positions
    ref_indices = {pos - 1 for pos in ref_positions}
    seq_rows: list[tuple[str, list[str]]] = []
    for i, (name, seq) in enumerate(seqs):
        if i in ref_indices:
            continue
        row = list(seq.upper()[:aln_len])
        row += ["-"] * (aln_len - len(row))
        seq_rows.append((name, row))

    # Column labels: 1-based position in the reference (skip gap columns)
    col_labels: list[tuple[int, str]] = []
    ref_pos = 0
    for i, base in enumerate(ref_row):
        if base != "-":
            ref_pos += 1
            if ref_pos == 1 or ref_pos % 10 == 0:
                col_labels.append((i, str(ref_pos)))

    label = Path(path).stem
    return Panel(label, ref_row, seq_rows, aln_len, col_labels,
                 extra_ref_rows=extra_ref_rows or None)
