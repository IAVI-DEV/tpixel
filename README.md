# tpixel

Pixel-block alignment viewer for hundreds of sequences. Renders Roark-style PIXEL plots: grey=match, red=substitution, black=gap.

## Install

```bash
pip install tpixel
```

## Usage

```bash
# Protein alignment (AA mode)
tpixel --fasta examples/env_protein_aligned.fasta -o examples/test_hiv_pixel.png

# Nucleotide alignment (NT mode — auto-detected)
tpixel --fasta examples/env_codon_aligned.fasta -o examples/test_hiv_nt_pixel.png
```

HIV mode is auto-detected when the alignment contains HxB2 and a `*_ref` sequence. Force with `--hiv` or `--no-hiv`. Nucleotide vs amino-acid mode is auto-detected from sequence content. Force with `--nt` or `--aa`.

## Examples

Both alignments contain the same 33 sequences (HxB2, animal1_ref, and 31 sample sequences from 7 animals):

| Sequence | Type |
|----------|------|
| HxB2 | Coordinate reference |
| animal1_ref | Parental lineage reference |
| animal1_s1 .. s7 | animal1 samples (7) |
| animal2_s1 .. s5 | animal2 samples (5) |
| animal3_s1 .. s8 | animal3 samples (8) |
| animal4_s1 | animal4 sample (1) |
| animal5_s1 .. s5 | animal5 samples (5) |
| animal6_s1 .. s2 | animal6 samples (2) |
| animal7_s1 .. s3 | animal7 samples (3) |

### Protein alignment (AA)

`env_protein_aligned.fasta` — 887 columns of HIV-1 Env protein. Each AA = 1 pixel.

![PIXEL plot — protein](https://raw.githubusercontent.com/tmsincomb/tpixel/main/examples/test_hiv_pixel.png)

### Nucleotide alignment (NT)

`env_codon_aligned.fasta` — 2661 columns of codon-aligned HIV-1 Env DNA. Each NT = 1 pixel. PNGS markers are detected by translating internally and mapped back to NT columns.

![PIXEL plot — nucleotide](https://raw.githubusercontent.com/tmsincomb/tpixel/main/examples/test_hiv_nt_pixel.png)
