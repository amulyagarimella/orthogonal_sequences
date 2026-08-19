# orthoseq

**Design DNA sequences that won't be mistaken for anything else in your sample.**

Tell it how many sequences you want and how long they should be, and it builds a
set where

- **no two sequences resemble each other** — any pair differs at at least half
  their positions, and you can ask for more,
- **none of them appear in a genome you care about** — nothing survives that
  `blastn` can match to human, mouse, or any other database you name,
- **none of them collide with a list you already have** — not equal to, contained
  in, or containing any sequence you blacklist, read in either direction and on
  either strand.

Use it for barcodes, spike-ins, and synthetic targets. It can also add to a set
you built earlier, score sequences that came from somewhere else, and design 10x
Flex v2 probes against them.

```bash
python3 orthoseq.py -n 100 -l 100 --ref-db human_genome -o designed.fa
```

> This is the full reference. For the summary, the metrics, and the defaults on
> one page, see **[README.md](README.md)**.

---

## Contents

- [Install](#install)
- [Quick start](#quick-start)
- [Usage](#usage)
  - [Generating sequences](#generating-sequences)
  - [Adding to a set you already have](#adding-to-a-set-you-already-have)
  - [Scoring sequences you already have](#scoring-sequences-you-already-have)
  - [Designing Flex v2 probes](#designing-flex-v2-probes)
- [Flag reference](#flag-reference)
- [Picking thresholds](#picking-thresholds)
- [Installing reference genomes](#installing-reference-genomes)
- [Troubleshooting](#troubleshooting)
- [How it works](#how-it-works)

---

## Install

orthoseq itself needs no install: it is one file, Python 3.9 or newer, and
imports nothing outside the standard library. Clone the repo and run
`orthoseq.py`.

You do need **NCBI BLAST+** — the `blastn` and `makeblastdb` programs — on your
`PATH`. orthoseq checks for both at startup and stops with an error if either is
missing. Reference genomes are a separate download; see
[Installing reference genomes](#installing-reference-genomes).

```bash
python3 orthoseq.py --help            # generating sequences
python3 orthoseq.py report --help     # scoring existing sequences
python3 orthoseq.py probes --help     # designing probes
```

> **On macOS, use NCBI's BLAST+ build rather than Homebrew's.** Homebrew's
> `blast` cannot open the databases NCBI distributes. If you see an
> `MDB_INVALID` error, go to [Troubleshooting](#the-lmdb-error-on-macos).

---

## Quick start

**Make 100 sequences of 100 nt that don't occur in human or mouse:**

```bash
python3 orthoseq.py -n 100 -l 100 \
  --ref-db human_genome --ref-db mouse_genome \
  --exclude my_existing_seqs.fa \
  --gc 0.5 --threads 10 \
  -o designed.fa
```

You get `designed.fa` and a companion `designed.fa.tsv`. The TSV has one row per
sequence: its GC, its longest homopolymer, its 3-mer entropy, how far it is from
its closest neighbour in the set, and the best reference hit it still has — so
you can see exactly how much margin you ended up with.

**Grow that set to 150 later:**

```bash
python3 orthoseq.py -n 150 --resume designed.fa --ref-db human_genome
```

**Score sequences that came from somewhere else:**

```bash
python3 orthoseq.py report barcodes.txt --ref-db human_genome -o barcodes.tsv
```

**Design 10x Flex v2 probes against a set of targets:**

```bash
python3 orthoseq.py probes designed.fa -o probes
```

Every command **exits `1` if it fell short** — fewer than `n` sequences found,
or at least one probe target that produced no usable pair — and `0` otherwise.
Output files are written either way, so a shell script can check the exit code
and still inspect what came out.

---

## Usage

Three commands: the default one generates sequences, `report` scores sequences,
and `probes` designs probes.

### Generating sequences

The default command — no subcommand word. You must pass `-n` (how many
sequences) and `-l` (how many bases each). Everything else has a default that
works.

```bash
python3 orthoseq.py -n 100 -l 100 -o designed.fa
```

Add constraints as you need them:

| To do this | Use |
|---|---|
| Reject anything matching a genome | `--ref-db human_genome` (repeatable) |
| Reject anything matching a small FASTA | `--ref-fasta plasmid.fa` (temporary database, deleted afterwards) |
| Stay clear of sequences you already have | `--exclude existing.fa` |
| Hit a particular GC | `--gc 0.5 --gc-tol 0.05` |
| Push the sequences further apart | `--min-hamming 0.65` |
| Keep out a restriction site | `--avoid-motif GAATTC` (repeatable; its reverse complement is blocked too) |
| Pin one position to one base | `--fix-base 25:A` (position counts from 0; repeatable) |
| Get the same output every run | `--seed 1` |

If it cannot find `n` sequences it still writes the ones it found, prints a
warning naming the settings to loosen, and exits `1`. See
[Not enough sequences found](#not-enough-sequences-found).

---

### Adding to a set you already have

`--resume` takes a set you made earlier and tops it up. The sequences already in
the file are kept exactly as they are and count toward `-n`. New sequences are
checked against them the same way they are checked against each other, so the
finished file — old and new together — still meets `--min-hamming` on both
strands.

```bash
# 100 sequences yesterday, 150 wanted today
python3 orthoseq.py -n 150 --resume designed.fa \
  --ref-db human_genome --min-hamming 0.6 --oversample 8 --max-rounds 12
```

What to expect:

- **`-n` is how many you want in the end**, not how many to add. The command
  above adds 50.
- **Do not pass `-l`.** The length is read from the file, and passing one that
  disagrees is an error.
- **Without `-o`, the `--resume` file is rewritten in place**, now holding its
  original sequences plus the new ones. The originals keep their names and their
  order; new sequences get numbered names that avoid the ones already taken.
- **Change any other setting you like between runs** — GC, rounds, oversampling,
  exclusion list, reference databases all behave exactly as in a fresh run.
- **Sequences already in the file are never re-checked or dropped.** If they do
  not meet the `--min-hamming` you are asking for now, orthoseq names the closest
  pair in a warning and keeps going.
- **The file does not have to be FASTA.** A `.txt`, `.csv`, or `.tsv` list works
  too, so you can extend a set that came from another tool.

---

### Scoring sequences you already have

You do not need to generate anything to get a report. `orthoseq report` writes
the same TSV for any set of sequences:

```bash
python3 orthoseq.py report barcodes.txt --ref-db human_genome -o barcodes.tsv
```

**Input format.** A file starting with `>` is read as FASTA. Anything else is
read one record per line, split on commas and tabs: the first field made only of
`ACGTN` is the sequence, and the first field that isn't is its name. That covers
`.txt`, `.csv`, and `.tsv`; header rows are skipped, and sequences with no name
get `seq_00001`, `seq_00002`, and so on.

**Columns:**

| Column | What it holds |
|---|---|
| `name`, `length`, `sequence` | as read from the file |
| `n_count` | how many `N`s the sequence contains — see [How N is handled](#how-n-is-handled) |
| `gc` | fraction of G and C |
| `max_homopolymer` | longest run of one base |
| `max_tandem_copies` | most back-to-back copies of a 2- or 3-base unit |
| `entropy3` | Shannon entropy of its 3-mers, in bits |
| `min_hamming` / `nearest` | positions differing from the closest other sequence, and which one that is |
| `min_hamming_rc` / `nearest_rc` | the same, comparing against reverse complements |
| `best_ref_db` … `best_ref_alnlen` | strongest reference hit: which database, bitscore, E-value, percent identity, alignment length |

The four neighbour columns are `NA` when the file holds fewer than two
sequences, or when the sequences are not all the same length — you cannot count
differing positions between sequences of different lengths.

**Speed.** With no `--ref-db` or `--ref-fasta` nothing is blasted, but the
neighbour columns still compare every sequence against every other one, so the
time grows with the square of the set size. Measured on 50 nt sequences:

| sequences | time |
|---|---|
| 100 | 0.4 s |
| 1,000 | 0.9 s |
| 5,000 | 14 s |
| 10,000 | 57 s |

`--no-hamming-rc` skips the reverse-complement half of that work.

#### How N is handled

An `N` here means a **mixed synthesis position** — a spot where the oligo pool
carries all four bases in equal amounts — not a base nobody measured. So every
column has to decide what to do with it, and the rule is: **columns you read
report the average; filters that reject use the worst case.** A filter that says
"probably fine" is not doing its job.

| Column | What an N does | Reasoning |
|---|---|---|
| `gc` | adds 0.5, the average | a mixed position really is half G/C across the pool |
| `max_homopolymer` | extends a run — worst case | some molecule in the pool does have that run, so `AAANAAA` is a 7-base risk. One N cannot stand in for two different bases at once, so `AANCC` stays at 3 |
| `max_tandem_copies` | worst case, but each N is locked once a copy uses it | keeps `ANACAG` at 2 copies instead of 3 |
| `entropy3` | 3-mers containing an N are not counted | see below |
| `min_hamming` | counts as differing from everything, N included | worst case |

<details>
<summary>Why 3-mers with an N are dropped rather than averaged out</summary>

Averaging is the tempting choice here and it backfires. If you spread an
N-containing 3-mer across all the 3-mers it could be, an all-N 50-mer scores
**6.00 bits — higher than the 5.58-bit ceiling a real 50-mer can reach at all**,
and far above the 4.95 bits a random 50-mer averages. The one sequence you know
nothing about would come out as the most complex thing in the file. Dropping
ambiguous 3-mers instead gives an all-N sequence a score of 0.00, which trips the
low-complexity filter — the behaviour you want.

The Hamming rule closes a similar hole. Comparing packed sequences with XOR
alone treats N against N as a match, which made two entirely ambiguous sequences
look identical and got one of them thrown out as a duplicate.
</details>

Since all of this is inference, the report carries `n_count` so no row can hide
how much of it was inferred. Sequences without any `N` are scored exactly as
before. Generation never emits an `N`, and probe design skips any window
containing one, since you cannot target a base that isn't fixed.

---

### Designing Flex v2 probes

`orthoseq probes` turns target sequences into 10x GEM-X Flex v2 custom probe
pairs. Each pair is a 25 nt LHS probe and a 25 nt RHS probe that together make
up the reverse complement of one 50 nt stretch of the target; they bind end to
end and are then ligated into one molecule. The rules come from the 10x GEM-X
Flex v2 Custom Probe Design Technical Note, CG000839 Rev B.

```bash
python3 orthoseq.py probes designed.fa -o probes
```

Three files come out:

| File | What's in it |
|---|---|
| `probes.tsv` | one row per pair: where the site is, both halves, GC, homopolymer, entropy, junction, off-target margin, and the finished oligos |
| `probes_order.csv` | `Pool name,Sequence Name,Sequence` — ready to upload to IDT as an oPool, with `/5Phos/` written into every RHS sequence |
| `probes_probes.fa` | the same oligos as plain DNA, for blasting or for running back through `orthoseq report` |

#### Which half binds where

The ligated 50 nt probe is the reverse complement of the target site. Reverse
complementing flips the order, so the **first** 25 nt of the probe — the LHS —
binds the **3' half** of the target site, and the RHS binds the 5' half. That is
the easy thing to get backwards, so it is checked against 10x's own worked
example: run on the EGFP region from Appendix B, `orthoseq probes` returns their
published `EGFP-LHS-1` and `EGFP-RHS-1` sequences base for base.

Fixed handles are attached when the oligos are written out, following Table 1:

```
LHS   5'-CCTTGGCACCCGAGAATTCCA-[25 nt]-3'          partial Read 2S
RHS   /5Phos/-[25 nt]-CCCATATAAGAAA-3'             partial Constant Sequence
```

**Check you are using the right workflow.** That RHS tail belongs to the
**multiplex** workflow (CG000834 / CG000835) and is the default. For the 4
Samples Kit (CG000841) pass `--workflow singleplex`, which substitutes the
partial Capture Sequence 1, `CGGTCCTAGCAA`. The two tails are not
interchangeable. Order RHS probes 5'-phosphorylated.

#### The checks

Each check is applied to one 25 nt half at a time, because that is how 10x
specifies them.

| Check | Default | Where it comes from |
|---|---|---|
| GC of each half | `--gc-min 0.44` / `--gc-max 0.72` | the range 10x states |
| Longest single-base run | `--max-homopolymer 4` | "avoid homopolymer repeats" |
| Copies of a 2- or 3-base unit | `--max-tandem-copies 4` | "avoid … low complexity" |
| 3-mer entropy | `--min-entropy 4.0` | same. A 25 nt half holds 23 3-mers, so its **ceiling is log2(23) = 4.52** — do not carry over the 6.0 figure that applies to long sequences |
| T at the ligation junction | on; `--no-require-tn` turns it off | the LHS must end in T |
| Distance from off-targets | `--cross-mismatch-min 5` | at least 5 mismatches in one half or the other |
| Pairs per target | `--pairs-per-target 3`, non-overlapping | 10x recommends 3 |

**The junction T rejects more targets than every other check, and it is a
property of your targets, not something probe design can fix.** A T at position
25 of the probe requires an **A at position 25 of the 50 nt target site**. A 50
nt target has exactly one possible site, so about one target in four is even
eligible. In a test set of 100 random 50-mers, 22 had an A at position 25, and
16 of those also cleared the composition checks. Adding `--no-require-tn`
raises that to 67, and 10x does note that other junction motifs "can also
function effectively" — but if you are designing the targets yourself, it is
cheaper to put an A at position 25 up front and keep the recommended junction.

Targets of about 150 nt and up have room for the full three non-overlapping
pairs. `--min-gap` adds spacing on top of the required non-overlap.

#### Off-target screening

Two separate screens that answer different questions.

**`--cross-mismatch-min` compares probes against the other targets you supplied.**
It is on by default. Every 50 nt window of every target is checked on both
strands, and each pair is scored as the *larger* of its two half-mismatch counts
— because ligation needs both halves bound next to each other, so one
well-mismatched half is enough to stop it. A margin at or above the threshold
prints as a plain number; anything below prints as `mm@target±:offset`, naming
the target, strand, and position it comes too close to. Cost measured on 100
targets of 400 nt (300 pairs): 5.9 s with the screen, 0.7 s without. Short
targets are effectively free. `--cross-mismatch-min 0` turns it off.

**`--ref-db` / `--ref-fasta` compares probes against a real reference.** Off by
default. Each half is blasted, and the best hit is reported as a mismatch count
— substitutions plus every probe base the alignment never reached. It does not
subtract the site you meant to hit, so blasting against a database that contains
your own targets will correctly report 0 mismatches.

#### Making targets that are probe-designable from the start

`--gc` controls the GC of the whole sequence. The probe rules judge each 25 nt
half on its own. Nothing connects the two, so a 50-mer can sit at a fine 0.52
overall while its halves split 0.40 and 0.64 — passing generation, then failing
probe design.

This is common, not bad luck. In a 50-mer with exactly 25 G/C bases, the first
half holds 12.5 ± 1.79 of them, and the 0.44 floor sits at 11 — only 0.84
standard deviations away. Working it out exactly, **26% of such sequences have at
least one half outside the 0.44–0.72 range.**

**`--gc-window` fixes this** by holding every sub-window in range instead of just
the overall average:

```bash
python3 orthoseq.py -n 100 -l 50 --fix-base 25:A --gc-window 25:0.44:0.72 -o designed.fa
```

Write it as `W:MIN:MAX` to set the window width and its bounds, or as a bare `W`
to reuse the `--gc` / `--gc-tol` range. It is a rolling count, so it adds one
pass over each candidate.

Measured on sets of 100 × 50 nt, all generated with `--fix-base 25:A` so the
junction T is satisfied, then counting how many targets yielded at least one
probe pair. Ranges are over 7 seeds:

| generation settings | targets with a probe pair |
|---|---|
| `--gc 0.50 --gc-tol 0.05` (default) | 59–76 of 100 |
| `--gc 0.54 --gc-tol 0.04` (shifting the average) | 81–91 of 100 |
| `--gc-window 25:0.44:0.72` | 94–99 of 100 |
| `--gc-window 25:0.44:0.72 --min-entropy 5.2` | **99–100 of 100** |

Shifting the average stalls around 90 because it does nothing about the *gap*
between the two halves. Constraining the window removes GC as a reason for
failure altogether.

The last few need `--min-entropy` for the same reason, one level down:
generation measures 3-mer entropy across the whole 50-mer and requires 5.0,
while probe design measures each 25 nt half and requires 4.0. A sequence can
score 5.05 overall and still have a half at 3.88. Raising the whole-sequence
requirement to 5.2 drags the weak halves up with it. There is no
`--entropy-window` yet; it would work the same way `--gc-window` does.

---

## Flag reference

Flags for generating sequences — the default command. Run
`python3 orthoseq.py report --help` or `probes --help` for those.

| Flag | Default | What it does |
|---|---|---|
| `-n` / `--num` | *required* | How many sequences the finished set should hold |
| `-l` / `--length` | *required*, unless `--resume` supplies it | Length of each sequence, in bases |
| `-o` / `--out` | `orthoseq.fa`, or the `--resume` file | Output FASTA. The report goes to `<out>.tsv` unless `--report` says otherwise |
| `--ref-db` | — | BLAST database to avoid matching; repeatable |
| `--ref-fasta` | — | FASTA to avoid matching; a database is built for it and deleted afterwards |
| `--ref-evalue` | `1e-3` | Reject a candidate on any reference hit with an E-value this low or lower |
| `--min-hamming` | `0.5` | How many positions any two sequences must differ at. Below 1 it is a fraction of `-l`; 1 or above it is a base count; `0` turns the check off |
| `--no-hamming-rc` | off | Compare forward strands only, skipping reverse complements |
| `--self-k` | `0` (off) | Also reject a candidate sharing any k-mer with an accepted one |
| `--self-blast` | off | Also blast candidates against each other |
| `--self-bitscore` | `25` | With `--self-blast`, the alignment bitscore that counts as too similar |
| `--exclude` | — | File of sequences the output must stay clear of |
| `--exclude-orientations` | all four | Which of `fwd rc rev comp` to check the exclusion list in |
| `--exclude-k` | `14` | Reject a candidate sharing a k-mer this long with the exclusion list; `0` leaves only the exact equal/contains/contained-in check |
| `--gc` / `--gc-tol` | `0.5` / `0.05` | Target GC and how far from it is acceptable |
| `--gc-window` | — | Also hold GC in range across every W-base window (`--gc-window 25:0.44:0.72`) |
| `--max-homopolymer` | `4` | Longest run of a single base allowed |
| `--max-tandem-copies` | `4` | Most back-to-back copies of a 2- or 3-base unit allowed |
| `--min-entropy` | `5.0` | Lowest 3-mer entropy allowed, measured over the whole sequence |
| `--avoid-motif` | — | Forbid a motif and its reverse complement; repeatable |
| `--fix-base` | — | Hold one position at one base, counting from 0; repeatable |
| `--oversample` | `4.0` | Candidates to prepare per sequence still needed, each round. Truncated to a whole number |
| `--max-rounds` | `8` | Stop after this many rounds even if short |
| `--resume` | — | Top up an existing set to `-n`; supplies `-l` and the default `-o` |
| `--threads` | every core | Passed through to `blastn` |
| `--seed` | — | Fix the random seed so the run repeats exactly |

---

## Picking thresholds

Two settings decide whether a run succeeds: `--min-hamming` and `--ref-evalue`.
Both are easy to set to a number that sounds reasonable and is impossible.

### `--min-hamming`: compare it against what chance already gives you

Two random sequences of the same length already differ at about **three quarters
of their positions**, because at each position they agree only when they happen
to draw the same base — a 1-in-4 chance. For `l=100` that means an expected 75
differing positions with a standard deviation of 4.3.

So the default `--min-hamming 0.5`, which is 50 positions at `l=100`, sits 5.8
standard deviations below that average. Essentially every random pair clears it.
**It is a safety floor, not a design target.** Raise it toward 0.65–0.7 to
genuinely push the set apart, and expect more rounds. Past about 0.75 you are
asking for sequences that random DNA does not supply and the run will come up
short.

Below 1, `--min-hamming` is a fraction of `--length`; at 1 or above it is a
count of bases:

```bash
--min-hamming 0.5   # default: differ at half the positions or more
--min-hamming 30    # differ at 30 positions or more, whatever the length
--min-hamming 0     # turn the check off
```

**The same fraction is much harsher at shorter lengths**, because the rounding
can push it past the average. At `l=100`, `0.75` gives exactly 75 — the average
itself. At `l=50` it gives 37.5, which rounds to 38, landing just above the
average of 37.5. (Rounding is half-to-even, so `0.65 × 50 = 32.5` rounds *down*
to 32.)

Measured at `l=50`, `n=100`, three seeds each:

| `--min-hamming` | positions required | chance a random pair clears it | sequences found |
|---|---|---|---|
| 0.60 | 30 | 99.4% | 100, 100, 100 |
| 0.65 | 32 | 97.1% | 100, 100, 100 |
| 0.70 | 35 | 83.7% | 22, 25, 26 |
| 0.75 | 38 | 51.1% | 7, 6, 7 |

Notice how far apart the last two columns are at 0.70: five random pairs in six
clear the bar, yet the run finds only about 25 sequences. That is because every
accepted sequence has to clear *every* other accepted one at once, and the odds
of that collapse quickly as the set grows. **At `l=50`, 0.65 is the practical
ceiling** — not the 0.7 you would guess from the `l=100` numbers.

### `--ref-evalue`: why an E-value rather than a bitscore

Random sequence is not free of matches to a 3 Gbp genome. Against human plus
mouse together, the longest exact match a random 100-mer hits by pure chance is
around log₄(100 × 6×10⁹) ≈ **20 bp**, worth roughly 37 bits. Set a bitscore
cutoff below that and you reject every candidate you ever generate.

`--ref-evalue` avoids the problem, because an E-value already accounts for how
big the database is: the same setting means the same thing whether you screen
against a plasmid or against hg38. The default `1e-3` is about a 25 bp perfect
match, or roughly 30 bp at 90% identity. Use `1e-6` to be stricter and expect
more rounds; use `0.1` if candidates are running out.

Comparing sequences to each other is cheap next to this. **The genome screen is
what costs you.**

### `--min-entropy`: lower it for short sequences

`--min-entropy` defaults to 5.0, measured across the whole sequence. Shorter
sequences hold fewer 3-mers, so they cannot score as high — which means the same
default gets steadily harsher as `-l` shrinks, and below about 60 nt it is doing
most of the rejecting. Measured on random sequence:

| length | average `entropy3` | rejected by the default 5.0 |
|---|---|---|
| 40 | 4.73 | 98% |
| 50 | 4.95 | 61% |
| 60 | 5.11 | 22% |
| 70 | 5.23 | 5% |
| 100 | 5.46 | 0% |

You can see it in the run summary: a default 100 × 50 nt run reported 474
`low_complexity` rejections, while the same run at `l=100` reported none.

So **at `l=50` use about `--min-entropy 4.6`, and at `l=40` about `4.3`** — a bit
under the average for that length — unless you actually want to select for
high-complexity sequences, in which case the default is doing exactly what you
asked, just slowly. Probe design uses a separate, lower default of 4.0 because it
scores 25 nt halves.

### When counting differing positions isn't enough

Hamming distance compares position 1 to position 1, position 2 to position 2,
and so on. It cannot see a match that is shifted along, and it cannot handle
insertions or deletions at all. Two sequences can differ at 75 of 100 positions
and still share a 20-mer at an offset.

| Check | Flag | What it catches | qcbc equivalent |
|---|---|---|---|
| Differing positions | `--min-hamming` (on) | substitutions, at matching positions only | `pdist` |
| Shared k-mer | `--self-k 12` | identical stretches at any offset | `ambiguous` |
| blastn, all against all | `--self-blast` | similarity with insertions or deletions | — |

The last two rarely fire on random sequence: two random 100-mers typically share
an exact stretch of only about 6 bp, and over 300 trials the longest was 10. That
is why differing positions alone is the default. **Turn the other two on for
structured sequences or sets that came from elsewhere**, where shifted repeats
are a real possibility.

---

## Installing reference genomes

### Recommended: download pre-built databases from NCBI

This skips `makeblastdb`, which is the slow part.

```bash
mkdir -p ~/blastdb && cd ~/blastdb
update_blastdb.pl --source ncbi --decompress human_genome    # GRCh38, ~10 GB
update_blastdb.pl --source ncbi --decompress mouse_genome    # GRCm39,  ~8 GB
echo 'export BLASTDB=$HOME/blastdb' >> ~/.zshrc
```

Once `BLASTDB` is set you can name them directly: `--ref-db human_genome`. Allow
about 18 GB of disk for the two together.

<details>
<summary>If the download arrives under an accession name instead</summary>

Older versions of `update_blastdb.pl` unpack under the assembly accession
(`GCF_000001405.39_top_level`) rather than `human_genome`. Add an alias file in
`$BLASTDB` to get the short name back:

```
# ~/blastdb/human_genome.nal
TITLE Homo sapiens GRCh38.p13 [GCF_000001405.39] top level
DBLIST GCF_000001405.39_top_level
```
</details>

### Alternative: download a FASTA and build the database yourself

Use this for any organism NCBI does not pre-build.

```bash
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz                       # makeblastdb cannot read gzip
makeblastdb -in hg38.fa -dbtype nucl -out hg38 -title hg38   # 20-40 min
```

Use the **full** assembly — every chromosome plus the unplaced contigs and the
alternate haplotypes. Each extra base is one more thing your sequences get to
avoid. Don't worry about the soft-masking in UCSC files: orthoseq turns masking
off when it blasts.

For a small one-off reference, skip the database entirely and pass
`--ref-fasta plasmid.fa`. orthoseq builds one in a temporary directory and
deletes it when the run finishes.

---

## Troubleshooting

### Not enough sequences found

orthoseq writes the sequences it did find, warns, and exits `1`. Try these in
order:

1. **Lower `--min-hamming`.** This is usually the cause, especially below
   `l=100`. See [Picking thresholds](#picking-thresholds); at `l=50` anything
   above 0.65 will struggle.
2. **Lower `--min-entropy`** if your sequences are shorter than about 60 nt. The
   default of 5.0 is above what a random 50-mer averages, so it rejects most
   candidates at that length — see
   [`--min-entropy`](#--min-entropy-lower-it-for-short-sequences).
3. **Raise `--ref-evalue`** toward `0.1` when screening against a large genome.
4. **Widen `--gc-tol`.**
5. **Raise `--oversample` or `--max-rounds`.** These cost time, not quality: more
   candidates per round, or more rounds.

### The LMDB error on macOS

Homebrew's `blast` is built against a version of LMDB whose on-disk format
differs from the one NCBI uses to build the databases it distributes. It can
read databases it built itself, but it **fails on every pre-built NCBI
database**:

```
BLAST Database error: LMDB runtime error: mdb_env_open: MDB_INVALID: File is not an LMDB file
```

Downloading again does not help — a fresh copy fails the same way. Install
NCBI's own build and put it ahead of Homebrew on your `PATH`:

```bash
curl -O https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-aarch64-macosx.tar.gz
tar xzf ncbi-blast-2.17.0+-aarch64-macosx.tar.gz -C ~
echo 'export PATH="$HOME/ncbi-blast-2.17.0+/bin:$PATH"' >> ~/.zshrc
```

**Then check that it took.** The error message never says which `blastn`
produced it, so confirm the right one is winning:

```bash
which blastn        # must NOT be /opt/homebrew/bin/blastn
blastn -version     # 2.17.0+ from NCBI
```

Two ways this quietly fails:

- **The `.pkg` installer alone is not enough.** It puts `blastn` in
  `/usr/local/ncbi/blast/bin`, which on a stock Apple Silicon machine comes
  *after* `/opt/homebrew/bin` on the `PATH` — so Homebrew still wins and the LMDB
  error continues. Prepend the directory explicitly, as above.
- **Shells that aren't login shells may never read `~/.zshrc`.** Cron jobs, an
  editor's built-in terminal, a CI step, or a `subprocess` launched by another
  program can all see the old `PATH` while your own terminal looks fine. If
  orthoseq reports a BLAST failure you cannot reproduce by hand, run
  `which blastn` in both places before you start suspecting the database.

### A probe target produced no pairs

Nearly always the junction T — see [The checks](#the-checks). A 50 nt target has
only one possible probe site, so it needs an A at position 25. Either generate
targets with `--fix-base 25:A`, or drop the requirement with `--no-require-tn`.

---

## How it works

Each round generates candidates, filters them, and keeps the ones that survive.
Rounds repeat until `n` sequences are accepted or `--max-rounds` runs out.

1. **Generate** random candidates at the target GC, filling in any positions
   `--fix-base` has pinned.
2. **Cheap local checks** — GC (overall, and per window if `--gc-window` is set),
   homopolymer runs, 2- and 3-base tandem repeats, 3-mer entropy, and any
   forbidden motifs.
3. **Exclusion list** — reject a candidate equal to, contained in, or containing
   any listed sequence, in **all four orientations**, plus a shared-14-mer check
   that catches near misses the exact test would let through.
4. **Reference screen** — `blastn` against each reference database. Any hit with
   an E-value at or below `--ref-evalue` disqualifies the candidate. Run with
   `-dust no -soft_masking false`, so repetitive candidates get screened rather
   than skipped.
5. **Compare against each other** — count differing positions against every
   already-accepted sequence, and against their reverse complements unless
   `--no-hamming-rc` is set.

Steps 1–3 repeat until the round has `max(need × --oversample, 32)` candidates in
hand, so only sequences that already passed the free checks are ever blasted.

### The four orientations of the exclusion list

`--exclude` rejects a candidate that equals, is contained in, or contains any
listed sequence, in any of these four forms:

| Orientation | `ACGTT` becomes | What it means |
|---|---|---|
| `fwd` | `ACGTT` | exactly as listed |
| `rc` | `AACGT` | reverse complement — the opposite strand |
| `rev` | `TTGCA` | read backwards, bases unchanged |
| `comp` | `TGCAA` | bases swapped for their partners, order unchanged |

All four are checked by default and cost nothing extra per candidate: applying
any two of these transforms gives you another one from the same set (reverse
complement is just reverse followed by complement, and each undoes itself), so
the exclusion list is expanded once at startup instead of transforming every
candidate.

Narrow it with `--exclude-orientations fwd rc` if you only care about the two
strands. `rev` and `comp` do not correspond to anything DNA does on its own;
check them when the list is a set of *identifiers* whose bases must not be
reachable by rearranging them in any simple way, which is the usual reason to
care.

### Why counting differing positions, and why it is fast

The measure is the minimum number of differing positions between any two
sequences — the same one [qcbc](https://github.com/pachterlab/qcbc) calls
`pdist`, and `--no-hamming-rc` matches its `-rc` flag.

Each sequence is packed into a single Python integer and compared with XOR, so
the work happens in C rather than one base at a time. Comparing every sequence
against every other, on both strands, takes 7 ms for 100 sequences of 100 nt and
177 ms for 500 — though it does grow with the square of the count, which is why
the `report` timings above climb the way they do.

<details>
<summary>Why this beats generating 10× and pruning</summary>

It is generate-and-prune, but ordered so the expensive step sees the fewest
candidates. The free checks — composition, differing positions, shared k-mers —
run first and reject a lot: in a default `l=100` run, 346 candidates were thrown
out to produce a batch of 400 for blasting, so roughly 46% never reached BLAST.
What does reach it goes in one batch per round per database, not one query per
sequence.

It also **loops**, which a single oversampled pass cannot. Generate 10× up front
and you either waste most of the work or still come up short; rounds request
exactly the shortfall and check each new candidate against everything already
accepted.

`bedtools random` is the wrong tool for this job entirely: it samples intervals
*out of* the genome, so every sequence it hands you is already a perfect match to
the reference. You want DNA that was never in the genome, not DNA taken from it.
</details>

<details>
<summary>Is blastn the right tool here?</summary>

At this scale, yes. 100 sequences of 100 nt is 10 kb of query; the search takes
seconds and the database is built once. The alternatives only pay off if you
scale up a lot.

| Tool | Verdict |
|---|---|
| `blastn -task blastn` | **Use this.** Its word size of 11 will seed on any 11 bp exact match, and its E-values account for database size. |
| `bbduk.sh` (BBTools) | A faster k-mer screen, but it only finds exact k-mers, so it misses the diverged matches blastn catches. Needs about 30 GB of RAM to hold human plus mouse at k=31. |
| `mmseqs2 easy-search` | Fast on large inputs, but less sensitive than blastn on short nucleotide queries. Worth considering above ~10⁵ candidates. |
| `bowtie2` / `bwa` / `minimap2` | Wrong tool. Read aligners assume the query nearly matches, and will silently skip the weak matches that matter here. |
| `megablast` | Far too insensitive — a word size of 28 misses everything short. |

For structure rather than sequence — secondary structure, matched melting
temperatures — run the output through NUPACK or ViennaRNA. BLAST cannot tell you
about any of that.
</details>
