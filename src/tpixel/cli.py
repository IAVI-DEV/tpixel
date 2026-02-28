"""Click CLI for tpixel."""

import sys

import click

from tpixel.fasta import fasta_panel, read_fasta
from tpixel.renderer import render_panels


def _expand_stdin(paths: list[str]) -> list[str]:
    """If paths is ``['-']``, read file paths from stdin (one per line).

    Args:
        paths: List of file path strings. A single ``'-'`` triggers stdin reading.

    Returns:
        Expanded list of file paths.

    Examples:
        >>> _expand_stdin(["file1.fasta", "file2.fasta"])
        ['file1.fasta', 'file2.fasta']
        >>> _expand_stdin([])
        []
    """
    if paths and len(paths) == 1 and paths[0] == "-":
        return [line.strip() for line in sys.stdin if line.strip()]
    return list(paths)


def _auto_detect_hiv(seqs: list[tuple[str, str]]) -> bool:
    """Check if alignment contains HxB2 and a ``*_ref`` sequence.

    Args:
        seqs: Parsed (name, sequence) tuples from :func:`read_fasta`.

    Returns:
        ``True`` if both HxB2 and a ``*_ref`` sequence are present.

    Examples:
        >>> _auto_detect_hiv([("HxB2", "ACG"), ("animal1_ref", "ACG")])
        True
        >>> _auto_detect_hiv([("seq1", "ACG"), ("seq2", "ACG")])
        False
    """
    names = {n.split()[0] for n, _ in seqs}
    has_hxb2 = "HxB2" in names
    has_ref = any(n.endswith("_ref") for n in names)
    return has_hxb2 and has_ref


@click.command(
    context_settings={"help_option_names": ["-h", "--help"]},
    epilog="Use '-' to read file paths from stdin, e.g.:\n\n"
    "  find . -name '*.fasta' | tpixel --fasta - -o out.png",
)
@click.option(
    "--fasta",
    multiple=True,
    help="Aligned FASTA file(s) — each becomes a panel. Use '-' for stdin.",
)
@click.option("--columns", help="Column range for FASTA, 1-based inclusive (e.g. 1-120).")
@click.option(
    "-o",
    "--output",
    default="pixel.png",
    show_default=True,
    help="Output image path.",
)
@click.option("--dpi", type=int, default=300, show_default=True, help="Image resolution.")
@click.option(
    "--cell",
    type=float,
    default=None,
    help="Cell size in inches (reserved for future layout modes).",
)
@click.option(
    "--hiv/--no-hiv",
    default=None,
    help="Force HIV mode (HxB2 regions, PNGS, animal grouping). Auto-detected if omitted.",
)
@click.option(
    "--nt/--aa",
    default=None,
    help="Force nucleotide or amino-acid mode. Auto-detected if omitted.",
)
@click.option(
    "--ref-pos",
    default="1,2",
    show_default=True,
    help="Comma-separated 1-based positions of reference sequences. "
    "Last position is the primary reference; earlier ones are extra reference rows.",
)
@click.option(
    "--title",
    default=None,
    help="Title displayed above the plot.",
)
def main(fasta, columns, output, dpi, cell, hiv, nt, ref_pos, title):
    """Pixel-block alignment viewer for hundreds of sequences.

    Renders Roark-style PIXEL plots: grey=match, red=substitution, black=gap.
    Each sequence is a thin row of colored blocks — no text in cells.

    HIV mode is auto-detected when the alignment contains HxB2 and a *_ref
    sequence. Force with --hiv or --no-hiv.
    """
    fasta_paths = _expand_stdin(list(fasta))

    if not fasta_paths:
        raise click.UsageError("Provide --fasta")

    ref_positions = [int(x) for x in ref_pos.split(",")]

    panels = []
    col_start, col_end = None, None
    if columns:
        parts = columns.split("-")
        if len(parts) != 2 or not parts[0].strip().isdigit() or not parts[1].strip().isdigit():
            raise click.UsageError("--columns must be START-END (e.g. 1-120)")
        col_start = int(parts[0].strip())
        col_end = int(parts[1].strip())

    for fasta_path in fasta_paths:
        seqs = read_fasta(fasta_path)
        use_hiv = hiv if hiv is not None else _auto_detect_hiv(seqs)

        if use_hiv:
            from tpixel.hiv import hiv_panel

            seq_type = None
            if nt is True:
                seq_type = "NT"
            elif nt is False:
                seq_type = "AA"
            panel = hiv_panel(fasta_path, ref_positions=ref_positions, seq_type=seq_type)
        else:
            panel = fasta_panel(fasta_path, col_start, col_end, ref_positions=ref_positions)

        if title:
            panel.title = title
        panels.append(panel)

    render_panels(panels, output, dpi=dpi, cell=cell)
