"""Tests for tpixel CLI."""

from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from tpixel.cli import _auto_detect_hiv, _expand_stdin, main


class TestExpandStdin:
    def test_passthrough(self):
        assert _expand_stdin(["a.fasta", "b.fasta"]) == ["a.fasta", "b.fasta"]

    def test_empty(self):
        assert _expand_stdin([]) == []


class TestAutoDetectHiv:
    def test_detects_hiv(self):
        seqs = [("HxB2", "ACGT"), ("animal1_ref", "ACGT"), ("animal1_s1", "ACGT")]
        assert _auto_detect_hiv(seqs) is True

    def test_no_hxb2(self):
        seqs = [("ref", "ACGT"), ("animal1_ref", "ACGT")]
        assert _auto_detect_hiv(seqs) is False

    def test_no_ref_suffix(self):
        seqs = [("HxB2", "ACGT"), ("s1", "ACGT")]
        assert _auto_detect_hiv(seqs) is False

    def test_empty_seqs(self):
        assert _auto_detect_hiv([]) is False


class TestCLI:
    def _write_fasta(
        self, tmp_path: Path, seqs: list[tuple[str, str]], name: str = "test.fasta"
    ) -> Path:
        path = tmp_path / name
        with open(path, "w", encoding="utf-8") as fh:
            for seq_name, seq in seqs:
                fh.write(f">{seq_name}\n{seq}\n")
        return path

    def test_no_fasta_error(self):
        runner = CliRunner()
        result = runner.invoke(main, [])
        assert result.exit_code != 0
        assert "Provide --fasta" in result.output

    def test_basic_render(self, tmp_path):
        fasta = self._write_fasta(tmp_path, [("ref", "ACGTACGT"), ("s1", "ACCTACGT")])
        out = tmp_path / "out.png"
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta", str(fasta), "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()
        assert out.stat().st_size > 0

    def test_columns_flag(self, tmp_path):
        fasta = self._write_fasta(tmp_path, [("ref", "ACGTACGTACGT"), ("s1", "ACCTACGTACGT")])
        out = tmp_path / "cols.png"
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta", str(fasta), "--columns", "2-5", "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_invalid_columns_error(self, tmp_path):
        fasta = self._write_fasta(tmp_path, [("ref", "ACGT"), ("s1", "ACGT")])
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta", str(fasta), "--columns", "bad"])
        assert result.exit_code != 0
        assert "--columns must be START-END" in result.output

    def test_title_flag(self, tmp_path):
        fasta = self._write_fasta(tmp_path, [("ref", "ACGTACGT"), ("s1", "ACCTACGT")])
        out = tmp_path / "title.png"
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta", str(fasta), "--title", "My Title", "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_hiv_mode_forced(self, tmp_path):
        seqs = [
            ("HxB2", "ACGTACGTAC" * 9),
            ("animal1_ref", "ACGTACGTAC" * 9),
            ("animal1_s1", "ACCTACGTAC" * 9),
        ]
        fasta = self._write_fasta(tmp_path, seqs)
        out = tmp_path / "hiv.png"
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta", str(fasta), "--hiv", "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_no_hiv_mode_forced(self, tmp_path):
        seqs = [
            ("HxB2", "ACGTACGT"),
            ("animal1_ref", "ACGTACGT"),
            ("animal1_s1", "ACCTACGT"),
        ]
        fasta = self._write_fasta(tmp_path, seqs)
        out = tmp_path / "nohiv.png"
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta", str(fasta), "--no-hiv", "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_hiv_auto_detect(self, tmp_path):
        seqs = [
            ("HxB2", "ACGTACGTAC" * 9),
            ("animal1_ref", "ACGTACGTAC" * 9),
            ("animal1_s1", "ACCTACGTAC" * 9),
        ]
        fasta = self._write_fasta(tmp_path, seqs)
        out = tmp_path / "autodetect.png"
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta", str(fasta), "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_multiple_fasta_panels(self, tmp_path):
        f1 = self._write_fasta(tmp_path, [("ref", "ACGTACGT"), ("s1", "ACCTACGT")], "a.fasta")
        f2 = self._write_fasta(tmp_path, [("ref", "ACGTACGT"), ("s2", "ACGTACGA")], "b.fasta")
        out = tmp_path / "multi.png"
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta", str(f1), "--fasta", str(f2), "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()
        assert out.stat().st_size > 0

    def test_nt_flag(self, tmp_path):
        seqs = [
            ("HxB2", "ACGTACGTAC" * 9),
            ("animal1_ref", "ACGTACGTAC" * 9),
            ("animal1_s1", "ACCTACGTAC" * 9),
        ]
        fasta = self._write_fasta(tmp_path, seqs)
        out = tmp_path / "nt.png"
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta", str(fasta), "--hiv", "--nt", "-o", str(out)])
        assert result.exit_code == 0, result.output

    def test_aa_flag(self, tmp_path):
        seqs = [
            ("HxB2", "MWLKFHRD" * 5),
            ("animal1_ref", "MWLKFHRD" * 5),
            ("animal1_s1", "MWLKFHRE" * 5),
        ]
        fasta = self._write_fasta(tmp_path, seqs)
        out = tmp_path / "aa.png"
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta", str(fasta), "--hiv", "--aa", "-o", str(out)])
        assert result.exit_code == 0, result.output

    def test_svg_output(self, tmp_path):
        fasta = self._write_fasta(tmp_path, [("ref", "ACGTACGT"), ("s1", "ACCTACGT")])
        out = tmp_path / "out.svg"
        runner = CliRunner()
        result = runner.invoke(main, ["--fasta", str(fasta), "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()
