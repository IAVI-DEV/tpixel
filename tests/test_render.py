"""Tests for tpixel pixel-block rendering."""

from __future__ import annotations

import random

from tpixel.fasta import fasta_panel, read_fasta
from tpixel.models import Marker, Panel, Region, SeqGroup
from tpixel.renderer import render_panels


class TestReadFasta:
    def test_roundtrip(self, write_fasta):
        seqs = [("ref", "ACGTACGT"), ("s1", "ACCTACGT")]
        path = write_fasta(seqs)
        result = read_fasta(path)
        assert len(result) == 2
        assert result[0] == ("ref", "ACGTACGT")

    def test_empty_file(self, tmp_path):
        path = tmp_path / "empty.fasta"
        path.write_text("")
        assert read_fasta(str(path)) == []


class TestFastaPanel:
    def test_panel_shape(self, write_fasta):
        seqs = [("ref", "ACGTACGT"), ("s1", "ACCTACGT"), ("s2", "ACGTACGA")]
        path = write_fasta(seqs)
        panel = fasta_panel(path)
        assert panel.total_cols == 8
        assert len(panel.ref_row) == 8
        assert len(panel.seq_rows) == 2

    def test_column_slicing(self, write_fasta):
        seqs = [("ref", "ACGTACGTACGT"), ("s1", "ACCTACGTACGT")]
        path = write_fasta(seqs)
        panel = fasta_panel(path, col_start=2, col_end=5)
        assert panel.total_cols == 4
        assert panel.ref_row == ["C", "G", "T", "A"]


class TestRenderPanels:
    def test_basic_pixel_render(self, write_fasta, output_dir):
        seqs = [("ref", "ACGTACGT"), ("s1", "ACCTACGT"), ("s2", "ACGTACGA")]
        path = write_fasta(seqs)
        panel = fasta_panel(path)
        out = output_dir / "basic_pixel.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_large_alignment(self, write_fasta, output_dir):
        """Render 200 sequences to verify it scales."""
        random.seed(42)
        ref = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        lines = [("ref", ref)]
        for i in range(200):
            mutant = list(ref)
            for j in random.sample(range(len(ref)), k=8):
                mutant[j] = random.choice([b for b in "ACGT" if b != ref[j]])
            for j in random.sample(range(len(ref)), k=3):
                mutant[j] = "-"
            lines.append((f"seq_{i:03d}", "".join(mutant)))
        path = write_fasta(lines)
        panel = fasta_panel(path)
        out = output_dir / "large_200_seqs.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_stacked_panels(self, write_fasta, output_dir):
        f1 = write_fasta(
            [("ref", "ACGTACGT"), ("s1", "ACCTACGT")], "group1.fasta"
        )
        f2 = write_fasta(
            [("ref", "ACGTACGT"), ("s2", "ACGTACGA")], "group2.fasta"
        )
        p1 = fasta_panel(f1)
        p2 = fasta_panel(f2)
        out = output_dir / "stacked_pixel.png"
        render_panels([p1, p2], str(out), dpi=150)
        assert out.exists()

    def test_custom_cell_size(self, write_fasta, output_dir):
        seqs = [("ref", "ACGTACGT"), ("s1", "ACCTACGT")]
        path = write_fasta(seqs)
        panel = fasta_panel(path)
        out = output_dir / "custom_cell.png"
        render_panels([panel], str(out), dpi=150, cell=0.05)
        assert out.exists()

    def test_all_gaps(self, write_fasta, output_dir):
        seqs = [("ref", "ACGT"), ("s1", "----")]
        path = write_fasta(seqs)
        panel = fasta_panel(path)
        out = output_dir / "all_gaps_pixel.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()

    def test_all_matches(self, write_fasta, output_dir):
        seqs = [("ref", "ACGTACGT"), ("s1", "ACGTACGT"), ("s2", "ACGTACGT")]
        path = write_fasta(seqs)
        panel = fasta_panel(path)
        out = output_dir / "all_matches_pixel.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()

    def test_output_formats(self, write_fasta, output_dir):
        seqs = [("ref", "ACGT"), ("s1", "ACGA")]
        path = write_fasta(seqs)
        panel = fasta_panel(path)
        for ext in ["svg", "pdf"]:
            out = output_dir / f"pixel_format.{ext}"
            render_panels([panel], str(out), dpi=150)
            assert out.exists()
            assert out.stat().st_size > 0


class TestRoarkLayout:
    """Test the full 7-layer Roark layout with regions, markers, groups."""

    def _make_roark_panel(self) -> Panel:
        """Build a panel that exercises all 7 layers."""
        random.seed(99)
        ref = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        ref_row = list(ref)
        aln_len = len(ref)

        groups = []
        for g in range(3):
            seqs = []
            for i in range(10):
                mutant = list(ref)
                for j in random.sample(range(aln_len), k=5):
                    mutant[j] = random.choice([b for b in "ACGT-" if b != ref[j]])
                seqs.append((f"g{g}_s{i}", mutant))
            groups.append(SeqGroup(name=f"animal_{g}", seqs=seqs))

        regions = [
            Region("R1", 0, 10, "#BBDEFB"),
            Region("R2", 10, 25, "#EEEEEE"),
            Region("R3", 25, 40, "#F8BBD0"),
        ]

        markers = [
            Marker(col=5, label="M5"),
            Marker(col=15, label="M15"),
            Marker(col=30, label="M30"),
        ]

        col_labels = [(i, str(i + 1)) for i in range(0, aln_len, 10)]

        return Panel(
            label="TestRef",
            ref_row=ref_row,
            seq_rows=[],
            total_cols=aln_len,
            col_labels=col_labels,
            regions=regions,
            markers=markers,
            marker_color="#4CAF50",
            groups=groups,
        )

    def test_full_roark_png(self, output_dir):
        panel = self._make_roark_panel()
        out = output_dir / "roark_full.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_full_roark_pdf(self, output_dir):
        panel = self._make_roark_panel()
        out = output_dir / "roark_full.pdf"
        render_panels([panel], str(out), dpi=300)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_regions_only(self, output_dir):
        """Panel with regions but no markers or groups."""
        ref_row = list("ACGTACGT")
        regions = [
            Region("A", 0, 4, "#BBDEFB"),
            Region("B", 4, 8, "#F8BBD0"),
        ]
        panel = Panel(
            label="ref",
            ref_row=ref_row,
            seq_rows=[("s1", list("ACCTACGT"))],
            total_cols=8,
            col_labels=[(0, "1")],
            regions=regions,
        )
        out = output_dir / "regions_only.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()

    def test_markers_only(self, output_dir):
        """Panel with markers but no regions or groups."""
        ref_row = list("ACGTACGT")
        markers = [Marker(col=2, label="X2"), Marker(col=5, label="X5")]
        panel = Panel(
            label="ref",
            ref_row=ref_row,
            seq_rows=[("s1", list("ACCTACGT"))],
            total_cols=8,
            col_labels=[(0, "1")],
            markers=markers,
        )
        out = output_dir / "markers_only.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()

    def test_grouped_sequences(self, output_dir):
        """Panel with grouped sequences but no regions or markers."""
        ref_row = list("ACGTACGT")
        groups = [
            SeqGroup("group_A", [("s1", list("ACCTACGT")), ("s2", list("ACGTACGA"))]),
            SeqGroup("group_B", [("s3", list("A-GTACGT"))]),
        ]
        panel = Panel(
            label="ref",
            ref_row=ref_row,
            seq_rows=[],
            total_cols=8,
            col_labels=[(0, "1")],
            groups=groups,
        )
        out = output_dir / "grouped_seqs.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()

    def test_large_roark(self, output_dir):
        """Stress test: 500-col alignment, 5 groups, 200 seqs total."""
        random.seed(77)
        ref = "".join(random.choices("ACGT", k=500))
        ref_row = list(ref)

        groups = []
        for g in range(5):
            seqs = []
            for i in range(40):
                mutant = list(ref)
                for j in random.sample(range(500), k=30):
                    mutant[j] = random.choice([b for b in "ACGT-" if b != ref[j]])
                seqs.append((f"g{g}_s{i}", mutant))
            groups.append(SeqGroup(name=f"animal_{g}", seqs=seqs))

        regions = [
            Region("V1", 0, 100, "#BBDEFB"),
            Region("C1", 100, 200, "#EEEEEE"),
            Region("V2", 200, 300, "#BBDEFB"),
            Region("C2", 300, 400, "#EEEEEE"),
            Region("V3", 400, 500, "#F8BBD0"),
        ]

        markers = [Marker(col=i, label=f"N{i}") for i in range(10, 500, 40)]

        col_labels = [(i, str(i + 1)) for i in range(0, 500, 50)]

        panel = Panel(
            label="StressRef",
            ref_row=ref_row,
            seq_rows=[],
            total_cols=500,
            col_labels=col_labels,
            regions=regions,
            markers=markers,
            groups=groups,
        )
        out = output_dir / "roark_stress.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()
        assert out.stat().st_size > 0


class TestTitleAndStats:
    """Test that title only appears when explicitly set, and stats appear in legend."""

    def test_no_title_by_default(self, output_dir):
        """Panel renders without title when title is not set."""
        ref_row = list("ACGTACGT")
        groups = [
            SeqGroup("sample_A", [("s1", list("ACCTACGT")), ("s2", list("ACGTACGA"))]),
            SeqGroup("sample_B", [("s3", list("A-GTACGT"))]),
        ]
        panel = Panel(
            label="ref",
            ref_row=ref_row,
            seq_rows=[],
            total_cols=8,
            col_labels=[(0, "1")],
            groups=groups,
        )
        assert panel.title is None
        out = output_dir / "no_title.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_title_when_set(self, output_dir):
        """Panel renders with title when explicitly provided."""
        ref_row = list("ACGTACGT")
        groups = [
            SeqGroup("sample_A", [("s1", list("ACCTACGT"))]),
            SeqGroup("sample_B", [("s2", list("A-GTACGT"))]),
        ]
        panel = Panel(
            label="ref",
            ref_row=ref_row,
            seq_rows=[],
            total_cols=8,
            col_labels=[(0, "1")],
            groups=groups,
            title="My Custom Title",
        )
        assert panel.title == "My Custom Title"
        out = output_dir / "with_title.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_hiv_panel_no_auto_title(self, tmp_path, output_dir):
        """HIV panel should not generate an auto title."""
        from tpixel.hiv import hiv_panel

        # Build a minimal HIV-like FASTA
        hxb2 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        ref = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        s1 = "ACCTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
        s2 = "ACGTACGAACGTACGTACGTACGTACGTACGTACGTACGT"

        fasta_path = tmp_path / "hiv_test.fasta"
        lines = []
        for name, seq in [("HxB2", hxb2), ("animal1_ref", ref), ("animal1_s1", s1), ("animal2_s1", s2)]:
            lines.append(f">{name}\n{seq}\n")
        fasta_path.write_text("".join(lines))

        panel = hiv_panel(str(fasta_path))
        assert panel.title is None

        # Explicitly set title to simulate --title
        panel.title = "Custom HIV Title"
        out = output_dir / "hiv_custom_title.png"
        render_panels([panel], str(out), dpi=150)
        assert out.exists()
        assert out.stat().st_size > 0


class TestPanel:
    def test_defaults(self):
        p = Panel("test", ["A", "C"], [], 2, [], ins_columns=None)
        assert p.ins_columns == set()
        assert p.label == "test"

    def test_total_seqs_flat(self):
        p = Panel("t", ["A"], [("s1", ["A"]), ("s2", ["A"])], 1, [])
        assert p.total_seqs == 2

    def test_total_seqs_grouped(self):
        g = [SeqGroup("g1", [("s1", ["A"])]), SeqGroup("g2", [("s2", ["A"]), ("s3", ["A"])])]
        p = Panel("t", ["A"], [], 1, [], groups=g)
        assert p.total_seqs == 3

    def test_effective_groups_from_flat(self):
        p = Panel("t", ["A"], [("s1", ["A"])], 1, [])
        eg = p.effective_groups
        assert len(eg) == 1
        assert eg[0].name == ""
        assert len(eg[0].seqs) == 1

    def test_seq_type_nt(self):
        p = Panel("t", list("ACGTACGT"), [], 8, [])
        assert p.seq_type == "NT"

    def test_seq_type_aa(self):
        p = Panel("t", list("MWLKFHRD"), [], 8, [])
        assert p.seq_type == "AA"

    def test_seq_type_nt_with_gaps(self):
        p = Panel("t", list("ACG-T.NU"), [], 8, [])
        assert p.seq_type == "NT"
