# orthoseq

**Design DNA sequences that won't be mistaken for anything else in your sample.**

Tell it how many sequences you want and how long they should be. It returns a set
where

- **no two sequences resemble each other** — every pair differs at at least half
  its positions, and you can require more,
- **none of them appear in a genome you care about** — nothing survives that
  `blastn` matches to human, mouse, or any database you name,
- **none of them collide with a list you supply** — not equal to, contained in, or
  containing any sequence you blacklist, in either direction and on either strand.

Use it for barcodes, spike-ins, and synthetic targets. It also extends an
existing set, scores sequences from elsewhere, and designs 10x Flex v2 probes.

h/t https://github.com/pachterlab/qcbc

One file, Python 3.9+, standard library only. Requires NCBI BLAST+ (`blastn`,
`makeblastdb`) on `PATH`.

```bash
# 100 sequences of 100 nt that don't occur in human or mouse
python3 orthoseq.py -n 100 -l 100 --ref-db human_genome --ref-db mouse_genome -o designed.fa

# grow that set to 150 later
python3 orthoseq.py -n 150 --resume designed.fa --ref-db human_genome

# score sequences from elsewhere
python3 orthoseq.py report barcodes.txt --ref-db human_genome -o barcodes.tsv

# design 10x Flex v2 probes against a set of targets
python3 orthoseq.py probes designed.fa -o probes
```

Each command writes its output and **exits `1` if it fell short** — fewer than
`n` sequences found, or a probe target with no usable pair.

📖 **[DETAILED_README.md](DETAILED_README.md)** — full documentation: reference
genomes, off-target screening, probe design rules, troubleshooting.

---

## The metrics

The columns of the `report` TSV, in order. Generation filters on the same
measurements.

| Column | What it is | How it's computed |
|---|---|---|
| `name` | sequence name | from the FASTA header, or the non-sequence field of the line; otherwise `seq_00001`, `seq_00002`, … |
| `length` | number of bases | character count |
| `n_count` | number of `N`s | character count. `N` means a mixed synthesis position — all four bases, equal amounts. [What each metric assumes about it](DETAILED_README.md#how-n-is-handled) |
| `gc` | fraction G and C | (G + C) ÷ length. Each `N` counts 0.5 |
| `max_homopolymer` | longest run of one base | longest stretch of identical bases. `AAAAT` → 4 |
| `max_tandem_copies` | most back-to-back copies of a 2- or 3-base unit | at each position, take the next 2 and 3 bases and count immediate repeats. `ATATATAT` → 4. No repeat → 1, not 0 |
| `entropy3` | complexity, in bits | Shannon entropy of the counts of each distinct overlapping 3-mer. 3-mers containing `N` are excluded |
| `min_hamming` | differing positions vs. the closest other sequence | compare position by position against every other sequence, count differences, take the smallest. `NA` if fewer than two sequences, or lengths differ |
| `nearest` | name of that closest sequence | — |
| `min_hamming_rc` | same, vs. reverse complements | each sequence is also compared to the others' reverse complements |
| `nearest_rc` | name of that closest sequence | — |
| `best_ref_db` | which reference it hit | highest-bitscore `blastn` hit across every `--ref-db` and `--ref-fasta`. `NA` if no hit or no reference given |
| `best_ref_bitscore` | strength of that hit | from `blastn` |
| `best_ref_evalue` | how often a hit that strong occurs by chance | from `blastn`; accounts for database size |
| `best_ref_pident` | percent identity of that hit | from `blastn` |
| `best_ref_alnlen` | bases the hit covered | from `blastn` |
| `sequence` | the bases | as read |

### Three thresholds that trip people up

**`entropy3` limits depend on length.** A sequence of length *L* has *L*−2
3-mers, capping `entropy3` at log2(*L*−2): 5.58 for a 50-mer, 4.52 for a 25 nt
probe half. Random sequence averages 4.95 at 50 nt, 5.46 at 100 nt. Reusing a
threshold across lengths does not work.
[Details →](DETAILED_README.md#--min-entropy-lower-it-for-short-sequences)

**`min_hamming` is already high by chance.** Two random sequences differ at ~75%
of positions — 75 of 100, sd 4.3. Requiring 50 is 5.8 sd below that and rejects
almost nothing. Useful thresholds start around 0.65.
[Details →](DETAILED_README.md#--min-hamming-compare-it-against-what-chance-already-gives-you)

**`gc` is a whole-sequence average, so halves can still fail.** In a 50-mer with
exactly 25 G/C, each half has 12.5 ± 1.79. 26% of these sequences have a half
outside the 0.44–0.72 probe range despite a perfect 0.50 overall. Use
`--gc-window` to bound windows instead of the average.
[Details →](DETAILED_README.md#making-targets-that-are-probe-designable-from-the-start)

---

## Defaults

### Generating sequences

| Flag | Default | What it does |
|---|---|---|
| `-n` / `--num` | *required* | Sequences in the finished set |
| `-l` / `--length` | *required*, unless `--resume` supplies it | Bases per sequence |
| `-o` / `--out` | `orthoseq.fa`, or the `--resume` file | Output FASTA; report goes to `<out>.tsv` |
| `--gc` / `--gc-tol` | `0.5` / `0.05` | Target GC and tolerance, so 0.45–0.55 |
| `--gc-window` | off | Also bound GC across every W-base window (`25:0.44:0.72`) |
| `--max-homopolymer` | `4` | Rejects runs of 5+ identical bases |
| `--max-tandem-copies` | `4` | Rejects 5+ back-to-back copies of a 2- or 3-base unit |
| `--min-entropy` | `5.0` | Minimum `entropy3`. **Exceeds the 4.95 random average at 50 nt**, so it rejects most candidates below ~60 nt. Lower it there |
| `--min-hamming` | `0.5` | Positions two sequences must differ at. <1 = fraction of `-l`; ≥1 = base count; `0` = off |
| `--no-hamming-rc` | off | Compare forward strands only |
| `--ref-db` / `--ref-fasta` | none | Reference to avoid matching; both repeatable |
| `--ref-evalue` | `1e-3` | Rejects on a reference hit this strong or stronger — roughly a 25 bp perfect match |
| `--exclude` | none | File of sequences to stay clear of. Reads FASTA, `.txt`, `.csv`, `.tsv` |
| `--exclude-orientations` | all four | Which of `fwd rc rev comp` to check that file in |
| `--exclude-k` | `14` | Rejects on a shared 14-mer; `0` leaves only the exact equal/contains check |
| `--self-k` | `0` (off) | Rejects a candidate sharing a k-mer with an accepted one. Catches matches at an offset, which `min_hamming` cannot see |
| `--self-blast` / `--self-bitscore` | off / `25` | Blasts candidates against each other. Catches matches with insertions or deletions |
| `--avoid-motif` | none | Forbids a motif and its reverse complement; repeatable |
| `--fix-base` | none | Pins one position to one base, counting from 0; repeatable |
| `--oversample` | `4.0` | Candidates prepared per sequence still needed, per round |
| `--max-rounds` | `8` | Stops after this many rounds |
| `--threads` | every core | Passed to `blastn` |
| `--seed` | none | Fixes the random seed |

### Designing probes (`probes`)

Each check applies to one 25 nt half, as 10x specifies. **Source** distinguishes
defaults 10x states numerically from those where 10x gives only a qualitative
rule and the number is orthoseq's.

| Flag | Default | What it does | Source |
|---|---|---|---|
| `--gc-min` / `--gc-max` | `0.44` / `0.72` | GC range per half | **10x**: "GC content should be between 44 − 72% for each 25 bp probe half" |
| `--cross-mismatch-min` | `5` | Mismatches required against every other target supplied; `0` disables | **10x**: "at least five mismatches in at least one of the LHS or RHS probes" |
| `--pairs-per-target` | `3` | Non-overlapping pairs per target | **10x**: three recommended; pairs "should not overlap with each other" |
| `--no-require-tn` | off, T required | Drops the junction rule. **Rejects more targets than any other check**: T at probe position 25 requires an A at target position 25, and a 50 nt target has one possible site | **10x**: recommended, not mandatory — other motifs "can also function effectively" |
| `--workflow` | `multiplex` | RHS tail: partial Constant Sequence for multiplex (CG000834/835), partial Capture Sequence 1 for `singleplex` (CG000841) | **10x** Table 1. **Not interchangeable** |
| `--max-homopolymer` | `4` | Longest single-base run per half | orthoseq's number. 10x states only "Avoid homopolymer repeats" |
| `--max-tandem-copies` | `4` | Most 2- or 3-base tandem copies per half | orthoseq's number. 10x states only "Avoid overlap with annotated repeat or low complexity sequences" |
| `--min-entropy` | `4.0` | Minimum `entropy3` per half | orthoseq's number, same low-complexity rule. The 25 nt cap is 4.52, so this binds |
| `--min-gap` | `0` | Extra spacing between pairs, beyond non-overlap | orthoseq's; no 10x basis |
| `--ref-db` / `--ref-fasta` | none | Blasts each half against a reference | orthoseq's flags; 10x advises BLAST-ing probes against the reference transcriptome |

**Excluding the WTA probe set.** 10x requires custom probes not to overlap the
whole-transcriptome probes. Download the Chromium Human or Mouse Transcriptome
Probe Set from the 10x support site and pass the CSV to `--exclude` when
generating targets — orthoseq reads the `probe_seq` column directly, no
preprocessing:

```bash
python3 orthoseq.py -n 100 -l 50 --fix-base 25:A --gc-window 25:0.44:0.72 \
  --exclude Chromium_Human_Transcriptome_Probe_Set_v2.0.0_GRCh38-2024-A.csv \
  --ref-db human_genome -o targets.fa
```

**Two 10x rules orthoseq does not check**, as both concern real transcripts:
design against coding regions rather than UTRs, and avoid common SNPs (or keep
them ≥4 bp from the ligation junction).

### Scoring sequences (`report`)

Pass `--ref-db` or `--ref-fasta`
to fill the `best_ref_*` columns.

`--no-hamming-rc` skips the reverse-complement comparison (slows down larger sets).
