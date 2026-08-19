#!/usr/bin/env python3
"""
orthoseq — generate a set of mutually dissimilar DNA sequences that also avoid
one or more reference genomes and an optional user-supplied exclusion list.

Pipeline (per round):
    1. sample random candidates at a target GC
    2. cheap local filters (homopolymer / dinucleotide repeat / entropy / motifs)
    3. exclusion-list filter (exact identity / substring / superstring against
       every orientation of each listed sequence -- forward, reverse complement,
       reverse and complement -- plus a shared-k-mer near-match test)
    4. reference screen: blastn against each reference DB; any hit at or below
       --ref-evalue disqualifies the candidate
    5. mutual-dissimilarity screen: minimum pairwise Hamming distance between
       accepted sequences (optionally also against their reverse complements),
       with an optional shared-k-mer test and an optional all-vs-all blastn

Rounds repeat until n sequences are accepted or --max-rounds is exhausted.

`--resume FILE` picks up from a set you already have: its sequences are kept,
count toward -n, and constrain the new ones the same way freshly accepted ones
do. The sequence length is read off the file.

`orthoseq report INPUT` skips generation and writes the same per-sequence report
for sequences you already have (FASTA, or one per line in .txt/.csv/.tsv).

`orthoseq probes INPUT` designs 10x GEM-X Flex v2 custom probe pairs against a
set of sequences, following the guidelines in CG000839.

Requires: blastn, makeblastdb (NCBI BLAST+) on PATH. Python 3.9+, stdlib only.
"""

from __future__ import annotations

import argparse
import os
import random
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter
from dataclasses import dataclass, field
from math import log2

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

BLAST_FIELDS = "qseqid sseqid pident length mismatch gapopen qstart qend bitscore evalue"


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


#: The four strand/direction orientations of a sequence. They form a group under
#: composition (rc == rev o comp, each is its own inverse), which is why the
#: exclusion list can be expanded once instead of transforming every candidate.
ORIENTATIONS = ("fwd", "rc", "rev", "comp")


def orient(seq: str, which: str) -> str:
    if which == "fwd":
        return seq
    if which == "rc":
        return revcomp(seq)
    if which == "rev":
        return seq[::-1]
    if which == "comp":
        return seq.translate(COMPLEMENT)
    raise ValueError(f"unknown orientation {which!r}")


# --------------------------------------------------------------------------
# Hamming distance (equal-length sequences only)
# --------------------------------------------------------------------------


def encode(seq: str) -> int:
    """Pack a sequence into one big integer so XOR compares it in C."""
    return int.from_bytes(seq.encode(), "big")


def hamming_enc(x: int, y: int, nbytes: int) -> int:
    """Hamming distance between two encode()d sequences of nbytes each."""
    return nbytes - (x ^ y).to_bytes(nbytes, "big").count(0)


def n_mask(seq: str) -> int:
    """Bitmask of the ambiguous positions, one bit per base."""
    if "N" not in seq:
        return 0
    m = 0
    for i, c in enumerate(seq):
        if c == "N":
            m |= 1 << i
    return m


def both_n(mask_a: int, mask_b: int) -> int:
    """Positions ambiguous in both sequences -- the ones XOR alone misses."""
    return bin(mask_a & mask_b).count("1")


def hamming(a: str, b: str) -> int:
    """Worst-case distance: an N never counts as agreeing with anything.

    This is the readable definition of the rule. The hot paths inline it with
    encode()/n_mask() hoisted out of their loops -- see pairwise_hamming --
    so any change to the convention has to land here and there together.

    XOR already treats N vs a real base as a mismatch, but scores N vs N as a
    match, which would make two all-ambiguous sequences look identical and get
    one of them rejected as a near-duplicate. Those positions are added back.
    """
    if len(a) != len(b):
        raise ValueError("Hamming distance requires equal-length sequences")
    return (hamming_enc(encode(a), encode(b), len(a))
            + both_n(n_mask(a), n_mask(b)))


def resolve_min_hamming(value: float, length: int) -> int:
    """<1 means a fraction of the sequence length, >=1 an absolute count."""
    if value <= 0:
        return 0
    if value < 1:
        return max(1, round(value * length))
    return int(value)


class HammingGuard:
    """Rejects any sequence within min_dist of one already kept.

    Reverse-complement checking is free: the complement of each kept sequence is
    stored alongside it, and hamming(rc(a), b) == hamming(a, rc(b)) because
    reverse-complementing both operands preserves the distance.
    """

    def __init__(self, length: int, min_dist: int, check_rc: bool):
        self.n = length
        self.min_dist = min_dist
        self.check_rc = check_rc
        self.kept: list[int] = []

    def too_close(self, seq: str) -> int | None:
        """The offending distance, or None if the sequence is far enough away."""
        if self.min_dist <= 0:
            return None
        x = encode(seq)
        for y in self.kept:
            d = hamming_enc(x, y, self.n)
            if d < self.min_dist:
                return d
        return None

    def add(self, seq: str) -> None:
        self.kept.append(encode(seq))
        if self.check_rc:
            self.kept.append(encode(revcomp(seq)))

    def clone(self) -> "HammingGuard":
        g = HammingGuard(self.n, self.min_dist, self.check_rc)
        g.kept = list(self.kept)
        return g


def pairwise_hamming(
    records: list[tuple[str, str]], check_rc: bool
) -> tuple[dict[str, tuple[int, str]], dict[str, tuple[int, str]]]:
    """Nearest neighbour by Hamming distance, forward and reverse-complement.

    Returns (forward, revcomp) maps of name -> (distance, partner name). Both are
    empty when the set has fewer than two sequences or the lengths differ, since
    Hamming distance is undefined then.
    """
    if len(records) < 2:
        return {}, {}
    n = len(records[0][1])
    if any(len(s) != n for _, s in records):
        return {}, {}

    enc = [(name, encode(s)) for name, s in records]
    rcseqs = [revcomp(s) for _, s in records] if check_rc else []
    rcs = [encode(s) for s in rcseqs]
    # The N-vs-N correction is a no-op on unambiguous input, so the masks are
    # only built -- and only consulted -- when some sequence actually needs them.
    has_n = any("N" in s for _, s in records)
    masks = [n_mask(s) for _, s in records] if has_n else []
    rc_masks = [n_mask(s) for s in rcseqs] if has_n else []
    fwd: dict[str, tuple[int, str]] = {}
    rc: dict[str, tuple[int, str]] = {}

    def offer(table, who, dist, partner):
        cur = table.get(who)
        if cur is None or dist < cur[0]:
            table[who] = (dist, partner)

    for i in range(len(enc)):
        ni, xi = enc[i]
        for j in range(i + 1, len(enc)):
            nj, xj = enc[j]
            d = hamming_enc(xi, xj, n)
            if has_n:
                d += both_n(masks[i], masks[j])
            offer(fwd, ni, d, nj)
            offer(fwd, nj, d, ni)
            if check_rc:
                d = hamming_enc(xi, rcs[j], n)
                if has_n:
                    d += both_n(masks[i], rc_masks[j])
                offer(rc, ni, d, nj)
                offer(rc, nj, d, ni)
    return fwd, rc


# --------------------------------------------------------------------------
# sequence generation and local filters
# --------------------------------------------------------------------------


def random_seq(
    length: int,
    gc: float,
    rng: random.Random,
    fixed: dict[int, str] | None = None,
) -> str:
    at = (1.0 - gc) / 2.0
    g = gc / 2.0
    bases = rng.choices("ACGT", weights=[at, g, g, at], k=length)
    for pos, base in (fixed or {}).items():
        bases[pos] = base
    return "".join(bases)


def parse_gc_window(
    spec: str | None, gc: float, gc_tol: float, length: int
) -> tuple[int, float, float] | None:
    """`W` or `W:MIN:MAX` -> (window, min, max); W alone reuses --gc/--gc-tol."""
    if spec is None:
        return None
    parts = spec.split(":")
    if len(parts) not in (1, 3) or not parts[0].isdigit():
        sys.exit(f"[orthoseq] --gc-window expects W or W:MIN:MAX, got {spec!r}")
    w = int(parts[0])
    if not 1 <= w <= length:
        sys.exit(f"[orthoseq] --gc-window {w} must be between 1 and --length {length}")
    if len(parts) == 3:
        try:
            lo, hi = float(parts[1]), float(parts[2])
        except ValueError:
            sys.exit(f"[orthoseq] --gc-window bounds must be numbers, got {spec!r}")
    else:
        lo, hi = gc - gc_tol, gc + gc_tol
    if not 0.0 <= lo <= hi <= 1.0:
        sys.exit(f"[orthoseq] --gc-window needs 0 <= MIN <= MAX <= 1, got {lo}..{hi}")
    return w, lo, hi


def parse_fixed_bases(specs: list[str], length: int) -> dict[int, str]:
    """`POS:BASE` strings -> {position: base}, with 0-based positions."""
    fixed: dict[int, str] = {}
    for spec in specs:
        pos, _, base = spec.partition(":")
        base = base.upper()
        if not pos.lstrip("-").isdigit() or base not in ("A", "C", "G", "T"):
            sys.exit(f"[orthoseq] --fix-base expects POS:BASE, got {spec!r}")
        if not 0 <= int(pos) < length:
            sys.exit(f"[orthoseq] --fix-base position {pos} is outside 0..{length - 1}")
        fixed[int(pos)] = base
    return fixed


def gc_fraction(seq: str) -> float:
    """Expected GC fraction, counting each N as the degenerate position it is.

    An N is an equimolar ACGT mix, so it contributes 0.5 of a GC base on
    average; the denominator stays the full length. This is the one statistic
    where the expectation is the number you want -- the repeat and complexity
    checks below are gates, and take the worst case instead.
    """
    return (seq.count("G") + seq.count("C") + 0.5 * seq.count("N")) / len(seq)


def window_gc_range(seq: str, w: int) -> tuple[float, float]:
    """Lowest and highest GC fraction over every w-nt sub-window.

    Counted the same way gc_fraction does -- an N is a degenerate position
    worth half a GC base -- and accumulated as integer halves so a long rolling
    sum cannot drift.
    """
    weight = {"G": 2, "C": 2, "N": 1}
    run = sum(weight.get(c, 0) for c in seq[:w])
    lo = hi = run
    for i in range(w, len(seq)):
        run += weight.get(seq[i], 0) - weight.get(seq[i - w], 0)
        if run < lo:
            lo = run
        elif run > hi:
            hi = run
    return lo / (2 * w), hi / (2 * w)


def max_homopolymer(seq: str) -> int:
    """Longest single-base run, resolving each N to whatever extends it.

    Some molecule in the pool has that run, so it is the length worth flagging:
    AAANAAA is a 7-base risk, not two 3-base ones. An N cannot serve two bases
    at once, though, so AANCC stays at 3 -- the trailing Ns of a broken run are
    the only ones the next base may reclaim.
    """
    best = run = ntail = 0
    base: str | None = None  # the single non-N base in the current run
    for ch in seq:
        if ch == "N" or ch == base:
            run += 1
        else:
            run = ntail + 1  # restart, reclaiming the Ns just before this base
        if ch != "N":
            base = ch
        ntail = ntail + 1 if ch == "N" else 0
        best = max(best, run)
    return best


def max_tandem_repeat(seq: str, unit_sizes=(2, 3)) -> int:
    """Longest run of a tandem repeat, in copies, over the given unit sizes."""
    best = 1
    n = len(seq)
    # The unambiguous case is the overwhelmingly common one and the hot path for
    # both the generator and the probe scan, so it skips the pinning machinery
    # below -- with no N to resolve, that reduces to plain string equality.
    if "N" not in seq:
        for u in unit_sizes:
            for i in range(n - u):
                unit = seq[i : i + u]
                copies = 1
                j = i + u
                while seq[j : j + u] == unit:
                    copies += 1
                    j += u
                best = max(best, copies)
        return best

    for u in unit_sizes:
        for i in range(n - u):
            unit = list(seq[i : i + u])
            copies = 1
            j = i + u
            while j + u <= n:
                nxt = seq[j : j + u]
                if not all(a == b or "N" in (a, b) for a, b in zip(unit, nxt)):
                    break
                # An N is free to be any base, but only once: pinning it as soon
                # as a copy resolves it keeps ANACAG at 2 copies rather than 3.
                unit = [b if a == "N" else a for a, b in zip(unit, nxt)]
                copies += 1
                j += u
            best = max(best, copies)
    return best


def kmer_entropy(seq: str, k: int = 3) -> float:
    """Shannon entropy of the k-mer distribution, in bits per k-mer."""
    counts: dict[str, int] = {}
    total = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer:
            # Spreading an N-containing k-mer over its 4^j resolutions would
            # score an all-N sequence at the 6.0-bit ceiling -- more complex
            # than any real sequence can reach. Unknown is not evidence of
            # complexity, so ambiguous k-mers simply do not vote.
            continue
        counts[kmer] = counts.get(kmer, 0) + 1
        total += 1
    if total == 0:
        return 0.0
    return -sum((c / total) * log2(c / total) for c in counts.values())


def canonical_kmers(seq: str, k: int) -> set[str]:
    """Strand-independent k-mer set: each k-mer stored as min(kmer, revcomp).

    Ambiguous k-mers are skipped. An N-containing k-mer can never equal one from
    an ACGT-only sequence, so it is inert against a candidate; all it can do is
    match the same ambiguous k-mer in another N-containing sequence and report a
    conflict that the shared uncertainty, not shared sequence, produced.
    """
    out = set()
    for i in range(len(seq) - k + 1):
        km = seq[i : i + k]
        if "N" in km:
            continue
        out.add(min(km, revcomp(km)))
    return out


def composition_reason(
    seq: str, max_hp: int, max_tandem: int, min_entropy: float
) -> str | None:
    """Shared non-GC composition gates, or None if the sequence passes.

    LocalFilters and ProbeRules both apply these; only their GC rules differ,
    and those stay separate because target-plus-tolerance and an explicit
    min/max window disagree at the boundary in floating point.
    """
    if max_homopolymer(seq) > max_hp:
        return "homopolymer"
    if max_tandem_repeat(seq) > max_tandem:
        return "tandem_repeat"
    if kmer_entropy(seq, 3) < min_entropy:
        return "low_complexity"
    return None


@dataclass
class LocalFilters:
    gc_target: float = 0.5
    gc_tol: float = 0.05
    max_homopolymer: int = 4
    max_tandem_copies: int = 4
    min_entropy: float = 5.0  # bits over 3-mers; max is log2(64) = 6
    avoid_motifs: list[str] = field(default_factory=list)
    #: Optional per-sub-window GC rule. --gc/--gc-tol constrain the sequence as
    #: a whole, which leaves the two halves of a 50-mer free to split 0.40/0.64;
    #: this pins every window of gc_window nt, so downstream rules that judge a
    #: fragment (10x Flex reads each 25 nt probe half separately) are satisfied
    #: by construction rather than by luck.
    gc_window: int = 0
    gc_window_min: float = 0.0
    gc_window_max: float = 1.0

    def check(self, seq: str) -> str | None:
        """Return a rejection reason, or None if the sequence passes."""
        gc = gc_fraction(seq)
        if abs(gc - self.gc_target) > self.gc_tol:
            return f"gc={gc:.3f}"
        if self.gc_window and len(seq) >= self.gc_window:
            lo, hi = window_gc_range(seq, self.gc_window)
            if lo < self.gc_window_min or hi > self.gc_window_max:
                return f"gc_window={lo:.3f}-{hi:.3f}"
        reason = composition_reason(seq, self.max_homopolymer,
                                    self.max_tandem_copies, self.min_entropy)
        if reason:
            return reason
        for motif in self.avoid_motifs:
            if motif in seq or revcomp(motif) in seq:
                return f"motif:{motif}"
        return None


# --------------------------------------------------------------------------
# FASTA I/O
# --------------------------------------------------------------------------


def read_fasta(path: str) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks).upper()))
                name, chunks = line[1:].split()[0] or f"seq{len(records)}", []
            else:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    return records


def read_seqs(path: str) -> list[tuple[str, str]]:
    """Read sequences from FASTA, or from any one-record-per-line file.

    A file starting with '>' is FASTA. Otherwise each line is split on commas
    and tabs; the first ACGTN-only field is the sequence and the first other
    field, if any, is its name. Lines with no sequence-looking field (headers,
    blanks) are skipped. This covers .txt, .csv and .tsv without a flag.
    """
    with open(path) as fh:
        for line in fh:
            if line.strip():
                if line.lstrip().startswith(">"):
                    return read_fasta(path)
                break
        else:
            return []

    records: list[tuple[str, str]] = []
    with open(path) as fh:
        for line in fh:
            fields = [f.strip() for f in line.replace(",", "\t").split("\t") if f.strip()]
            seq = next((f.upper() for f in fields if not set(f.upper()) - set("ACGTN")), None)
            if seq is None:
                continue
            name = next((f.split()[0] for f in fields if f.upper() != seq), None)
            records.append((name or f"seq_{len(records) + 1:05d}", seq))
    return records


def infer_length(records: list[tuple[str, str]], source: str) -> int:
    """The single sequence length in a set, or exit if it is not single."""
    lengths = sorted({len(s) for _, s in records})
    if len(lengths) != 1:
        sys.exit(f"[orthoseq] {source} holds sequences of {len(lengths)} different "
                 f"lengths ({', '.join(str(x) for x in lengths[:5])}"
                 f"{' ...' if len(lengths) > 5 else ''}); resuming needs one length")
    return lengths[0]


def dedupe_names(records: list[tuple[str, str]]) -> list[tuple[str, str]]:
    """Make names unique, so a resumed set can be a BLAST subject safely."""
    used: set[str] = set()
    out = []
    for name, seq in records:
        base, i = name, 1
        while name in used:
            i += 1
            name = f"{base}_{i}"
        used.add(name)
        out.append((name, seq))
    return out


def assign_names(records: list[tuple[str, str]], taken=()) -> list[tuple[str, str]]:
    """Rename to seq_0001, seq_0002, ... skipping any name already spoken for."""
    used = set(taken)
    out = []
    i = 0
    for _, seq in records:
        while True:
            i += 1
            name = f"seq_{i:04d}"
            if name not in used:
                break
        used.add(name)
        out.append((name, seq))
    return out


def write_fasta(path: str, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


# --------------------------------------------------------------------------
# BLAST wrappers
# --------------------------------------------------------------------------


@dataclass
class Hit:
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    bitscore: float
    evalue: float


def _run_blastn(args: list[str], label: str, verbose: bool) -> list[Hit]:
    t0 = time.time()
    proc = subprocess.run(args, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"[orthoseq] blastn failed ({label}):\n{proc.stderr.strip()}")
    hits = []
    for line in proc.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 10:
            continue
        hits.append(Hit(f[0], f[1], float(f[2]), int(f[3]), int(f[4]),
                        float(f[8]), float(f[9])))
    if verbose:
        print(
            f"[orthoseq]   blastn {label}: {len(hits)} hits in {time.time() - t0:.1f}s",
            file=sys.stderr,
        )
    return hits


def blast_vs_db(
    query_fa: str,
    db: str,
    evalue: float,
    word_size: int,
    threads: int,
    verbose: bool,
) -> list[Hit]:
    args = [
        "blastn",
        "-task", "blastn",
        "-query", query_fa,
        "-db", db,
        "-word_size", str(word_size),
        "-evalue", str(evalue),
        "-dust", "no",
        "-soft_masking", "false",
        "-max_target_seqs", "5",
        "-max_hsps", "1",
        "-num_threads", str(threads),
        "-outfmt", f"6 {BLAST_FIELDS}",
    ]
    return _run_blastn(args, os.path.basename(db), verbose)


def blast_vs_subject(
    query_fa: str,
    subject_fa: str,
    word_size: int,
    threads: int,
    verbose: bool,
    evalue: float = 10.0,
) -> list[Hit]:
    args = [
        "blastn",
        "-task", "blastn",
        "-query", query_fa,
        "-subject", subject_fa,
        "-word_size", str(word_size),
        "-evalue", str(evalue),
        "-dust", "no",
        "-soft_masking", "false",
        "-max_hsps", "1",
        "-outfmt", f"6 {BLAST_FIELDS}",
    ]
    # -subject mode does not accept -num_threads in some BLAST+ builds.
    return _run_blastn(args, "self", verbose)


def require_tools(*tools: str) -> None:
    for tool in tools:
        if shutil.which(tool) is None:
            sys.exit(f"[orthoseq] {tool} not found on PATH")


def make_ref_dbs(
    ref_dbs: list[str], ref_fastas: list[str], tmpdir: str, verbose: bool
) -> list[str]:
    """Existing BLAST DBs, plus a temporary DB built for each reference FASTA."""
    dbs = list(ref_dbs)
    for fa in ref_fastas:
        db = os.path.join(tmpdir, os.path.basename(fa) + ".db")
        if verbose:
            print(f"[orthoseq] makeblastdb {fa}", file=sys.stderr)
        proc = subprocess.run(
            ["makeblastdb", "-in", fa, "-dbtype", "nucl", "-out", db],
            capture_output=True, text=True,
        )
        if proc.returncode != 0:
            sys.exit(f"[orthoseq] makeblastdb failed:\n{proc.stderr.strip()}")
        dbs.append(db)
    return dbs


# --------------------------------------------------------------------------
# report
# --------------------------------------------------------------------------


def write_report(
    path: str,
    records: list[tuple[str, str]],
    dbs: list[str],
    tmpdir: str,
    word_size: int,
    threads: int,
    verbose: bool,
    hamming_rc: bool = True,
) -> None:
    """Per-sequence composition stats, nearest neighbour, and best reference hit."""
    best_hit: dict[str, tuple[str, float, float, float, int]] = {}
    if dbs and records:
        query_fa = os.path.join(tmpdir, "report_query.fa")
        write_fasta(query_fa, records)
        for db in dbs:
            for h in blast_vs_db(query_fa, db, 10.0, word_size, threads, verbose):
                prev = best_hit.get(h.qseqid)
                if prev is None or h.bitscore > prev[1]:
                    best_hit[h.qseqid] = (
                        os.path.basename(db), h.bitscore, h.evalue, h.pident, h.length,
                    )
    fwd, rc = pairwise_hamming(records, hamming_rc)
    with open(path, "w") as fh:
        fh.write("name\tlength\tn_count\tgc\tmax_homopolymer\tmax_tandem_copies\tentropy3\t"
                 "min_hamming\tnearest\tmin_hamming_rc\tnearest_rc\t"
                 "best_ref_db\tbest_ref_bitscore\tbest_ref_evalue\t"
                 "best_ref_pident\tbest_ref_alnlen\tsequence\n")
        for name, s in records:
            b = best_hit.get(name)
            ref_cols = (
                [b[0], f"{b[1]:.1f}", f"{b[2]:.2g}", f"{b[3]:.1f}", str(b[4])]
                if b else ["NA"] * 5
            )
            ham_cols = []
            for table in (fwd, rc):
                hit = table.get(name)
                ham_cols += [str(hit[0]), hit[1]] if hit else ["NA", "NA"]
            fh.write("\t".join([name, str(len(s)), str(s.count("N")),
                                f"{gc_fraction(s):.3f}",
                                str(max_homopolymer(s)), str(max_tandem_repeat(s)),
                                f"{kmer_entropy(s, 3):.2f}"]
                               + ham_cols + ref_cols + [s]) + "\n")


# --------------------------------------------------------------------------
# exclusion list
# --------------------------------------------------------------------------


class ExclusionSet:
    """Exact substring/superstring membership plus a shared-k-mer near-match test.

    Every listed sequence is expanded into the requested orientations up front,
    so a candidate is rejected when it equals, is contained in, or contains any
    orientation of any listed sequence. Expanding the list is equivalent to
    transforming each candidate -- the orientations form a group -- and is
    cheaper, because the expansion happens once.
    """

    def __init__(
        self,
        records: list[tuple[str, str]],
        k: int,
        orientations: tuple[str, ...] = ORIENTATIONS,
    ):
        self.k = k
        self.orientations = tuple(orientations)
        self.variants: list[tuple[str, str]] = []
        seen: set[str] = set()
        for _, s in records:
            for o in self.orientations:
                v = orient(s, o)
                if v not in seen:
                    seen.add(v)
                    self.variants.append((v, o))
        self.kmers: set[str] = set()
        for v, _ in self.variants:
            self.kmers |= canonical_kmers(v, k)

    def reject(self, cand: str) -> str | None:
        for e, o in self.variants:
            if cand == e:
                return f"exclusion:identical_{o}"
            if cand in e:
                return f"exclusion:substring_{o}"
            if e in cand:
                return f"exclusion:superstring_{o}"
        if self.kmers and (canonical_kmers(cand, self.k) & self.kmers):
            return f"exclusion:shared_{self.k}mer"
        return None


# --------------------------------------------------------------------------
# selection
# --------------------------------------------------------------------------


def greedy_independent_set(
    ids: list[str],
    conflicts: dict[str, set[str]],
    limit: int,
    priority: dict[str, float],
) -> list[str]:
    """Greedy MIS: repeatedly take the lowest-priority-score, lowest-degree node."""
    remaining = set(ids)
    chosen: list[str] = []
    while remaining and len(chosen) < limit:
        node = min(
            remaining,
            key=lambda i: (
                len(conflicts.get(i, set()) & remaining),
                priority.get(i, 0.0),
                i,
            ),
        )
        chosen.append(node)
        remaining.discard(node)
        remaining -= conflicts.get(node, set())
    return chosen


# --------------------------------------------------------------------------
# 10x GEM-X Flex v2 probe design
# --------------------------------------------------------------------------

#: Handle sequences from CG000839 Rev B, Table 1. The LHS probe carries a
#: partial Read 2S on its 5' end; the RHS probe is 5'-phosphorylated and
#: carries a workflow-specific tail on its 3' end.
LHS_HANDLE = "CCTTGGCACCCGAGAATTCCA"
RHS_TAILS = {
    # multiplex Flex v2 (CG000834 / CG000835): partial Constant Sequence
    "multiplex": "CCCATATAAGAAA",
    # singleplex Flex v2, 4 Samples Kit (CG000841): partial Capture Sequence 1
    "singleplex": "CGGTCCTAGCAA",
}

PROBE_HALF = 25
PROBE_SITE = 2 * PROBE_HALF


@dataclass
class ProbeRules:
    """Per-half acceptance rules, defaulting to the CG000839 guidelines."""

    gc_min: float = 0.44
    gc_max: float = 0.72
    max_homopolymer: int = 4
    max_tandem_copies: int = 4
    #: bits over 3-mers. A 25 nt half has 23 3-mers, so the ceiling is
    #: log2(23) = 4.52, not the 6.0 that applies to a long sequence.
    min_entropy: float = 4.0
    require_tn: bool = True

    def check_half(self, half: str) -> str | None:
        gc = gc_fraction(half)
        if not (self.gc_min <= gc <= self.gc_max):
            return f"gc={gc:.3f}"
        return composition_reason(half, self.max_homopolymer,
                                  self.max_tandem_copies, self.min_entropy)


@dataclass
class ProbePair:
    """One LHS/RHS pair against a 50 nt site of one target."""

    target: str
    start: int          # 0-based offset of the 50 nt site in the target
    lhs: str            # 25 nt, reverse complement of site[25:50]
    rhs: str            # 25 nt, reverse complement of site[0:25]
    penalty: float      # tie-break score between valid windows; see scan_sites
    index: int = 0      # 1-based, assigned by select_pairs
    #: fewest mismatches at any unintended site, and where that site was.
    #: None when the cross-target screen was disabled.
    cross_mm: int | None = None
    cross_hit: str = ""

    @property
    def site(self) -> str:
        """The 50 nt target site, 5'->3' -- the halves are its reverse complement."""
        return revcomp(self.lhs + self.rhs)

    @property
    def tn(self) -> bool:
        """Whether the LHS ends in T, the recommended ligation junction."""
        return self.lhs.endswith("T")

    @property
    def name(self) -> str:
        return f"{self.target}_p{self.index}"

    def oligos(self, workflow: str) -> tuple[str, str]:
        """(LHS, RHS) as ordered, with handles; RHS needs a 5' phosphate."""
        return LHS_HANDLE + self.lhs, self.rhs + RHS_TAILS[workflow]


def probe_halves(site: str) -> tuple[str, str]:
    """Split the reverse complement of a 50 nt target site into LHS and RHS.

    The ligated probe is revcomp(site), so its first 25 nt -- the LHS probe --
    pair with the *3'* half of the site and the RHS probe with the 5' half.
    """
    probe = revcomp(site)
    return probe[:PROBE_HALF], probe[PROBE_HALF:]


def scan_sites(name: str, seq: str, rules: ProbeRules) -> list[ProbePair]:
    """Every 50 nt window of one target that yields an acceptable probe pair."""
    out: list[ProbePair] = []
    mid_gc = (rules.gc_min + rules.gc_max) / 2
    for start in range(len(seq) - PROBE_SITE + 1):
        # TN junction: T is the 3'-most base of the LHS probe, so the opposing
        # target base -- site[25] -- must be an A. It is one character and it
        # rejects roughly three windows in four, so it goes first.
        if rules.require_tn and seq[start + PROBE_HALF] != "A":
            continue
        site = seq[start : start + PROBE_SITE]
        if set(site) - set("ACGT"):
            continue
        lhs, rhs = probe_halves(site)
        if rules.check_half(lhs) or rules.check_half(rhs):
            continue
        tn = lhs.endswith("T")
        # Prefer mid-range GC, balanced halves, and high complexity; the site
        # itself is arbitrary, so this only breaks ties between valid windows.
        penalty = (
            abs(gc_fraction(lhs) - mid_gc)
            + abs(gc_fraction(rhs) - mid_gc)
            + abs(gc_fraction(lhs) - gc_fraction(rhs))
            + (0.0 if tn else 0.1)
            - 0.05 * (kmer_entropy(lhs, 3) + kmer_entropy(rhs, 3))
        )
        out.append(ProbePair(target=name, start=start, lhs=lhs, rhs=rhs,
                             penalty=penalty))
    return out


def select_pairs(cands: list[ProbePair], limit: int, min_gap: int) -> list[ProbePair]:
    """Best non-overlapping windows, greedily by penalty.

    Probe pairs against the same target must not overlap (CG000839), so a
    chosen window blocks everything within PROBE_SITE + min_gap of it.
    """
    chosen: list[ProbePair] = []
    for cand in sorted(cands, key=lambda c: (c.penalty, c.start)):
        if len(chosen) >= limit:
            break
        if all(abs(cand.start - c.start) >= PROBE_SITE + min_gap for c in chosen):
            chosen.append(cand)
    chosen.sort(key=lambda c: c.start)
    for i, c in enumerate(chosen, 1):
        c.index = i
    return chosen


def cross_screen(
    pairs: list[ProbePair],
    targets: list[tuple[str, str]],
    min_mm: int,
) -> None:
    """Annotate each pair with its closest unintended site in the input set.

    Ligation needs both halves bound at adjacent positions, so the unit of
    off-target risk is a competing 50 nt window rather than a lone 25-mer.
    CG000839 asks for at least five mismatches in at least one of the LHS or
    RHS probes, so each window scores as max(mm_lhs, mm_rhs) and the worst
    window over the whole set is what gets recorded.
    """
    if min_mm <= 0:
        return
    # Every 50 nt window of the set, packed once: each becomes the two 25 nt
    # halves a probe pair would have to match to bind there. Both strands are
    # included, in case the construct is present as duplex DNA rather than as
    # the transcript alone.
    windows: list[tuple[str, str, int, int, int]] = []
    for name, seq in targets:
        for sign in ("+", "-"):
            s = seq if sign == "+" else revcomp(seq)
            for j in range(len(s) - PROBE_SITE + 1):
                windows.append((
                    name, sign, j,
                    encode(s[j : j + PROBE_HALF]),
                    encode(s[j + PROBE_HALF : j + PROBE_SITE]),
                ))

    half, ham = PROBE_HALF, hamming_enc  # hoisted out of a many-million-iteration loop
    for pair in pairs:
        # what a target must look like for each half of this pair to bind
        want_rhs, want_lhs = encode(revcomp(pair.rhs)), encode(revcomp(pair.lhs))
        start, target = pair.start, pair.target
        worst, where = half, ""
        for name, sign, j, first, second in windows:
            # max(a, b) >= a, so a window whose first half is already no better
            # than the incumbent cannot improve it -- skip the second hamming.
            a = ham(want_rhs, first, half)
            if a >= worst:
                continue
            b = ham(want_lhs, second, half)
            mm = a if a > b else b
            if mm < worst:
                if j == start and sign == "+" and name == target:
                    continue  # the intended site, which trivially matches
                worst, where = mm, f"{name}{sign}:{j}"
                if worst == 0:
                    break
        pair.cross_mm, pair.cross_hit = worst, where


def cross_cell(pair: ProbePair, min_mm: int) -> str:
    """The off-target column: bare margin, or margin@where when it is too close."""
    if pair.cross_mm is None:
        return "NA"
    if pair.cross_mm >= min_mm:
        return str(pair.cross_mm)
    return f"{pair.cross_mm}@{pair.cross_hit}"


def write_probe_outputs(
    stem: str,
    pairs: list[ProbePair],
    workflow: str,
    min_mm: int,
    ref_hits: dict[str, tuple[str, int]],
) -> tuple[str, str, str]:
    """Write the report TSV, the ordering CSV, and a plain FASTA of the oligos."""
    tsv, csv, fa = stem + ".tsv", stem + "_order.csv", stem + "_probes.fa"
    rows = [(p, *p.oligos(workflow)) for p in pairs]

    with open(tsv, "w") as fh:
        fh.write("\t".join([
            "probe_pair", "target", "pair", "site_start", "site_end",
            "target_site", "lhs_probe", "rhs_probe", "lhs_gc", "rhs_gc",
            "lhs_homopolymer", "rhs_homopolymer", "lhs_entropy3", "rhs_entropy3",
            "tn_junction", "cross_target_mm", "lhs_best_ref", "lhs_ref_mm",
            "rhs_best_ref", "rhs_ref_mm", "lhs_oligo", "rhs_oligo_5phos",
        ]) + "\n")
        for p, lhs_o, rhs_o in rows:
            cols = [p.name, p.target, str(p.index), str(p.start),
                    str(p.start + PROBE_SITE), p.site, p.lhs, p.rhs,
                    f"{gc_fraction(p.lhs):.3f}", f"{gc_fraction(p.rhs):.3f}",
                    str(max_homopolymer(p.lhs)), str(max_homopolymer(p.rhs)),
                    f"{kmer_entropy(p.lhs, 3):.2f}", f"{kmer_entropy(p.rhs, 3):.2f}",
                    "yes" if p.tn else "no", cross_cell(p, min_mm)]
            for side in ("LHS", "RHS"):
                h = ref_hits.get(f"{p.name}_{side}")
                cols += [h[0], str(h[1])] if h else ["NA", "NA"]
            cols += [lhs_o, "/5Phos/" + rhs_o]
            fh.write("\t".join(cols) + "\n")

    pool = os.path.basename(stem)
    with open(csv, "w") as fh:
        fh.write("Pool name,Sequence Name,Sequence\n")
        for p, lhs_o, rhs_o in rows:
            fh.write(f"{pool},{p.name}_LHS,{lhs_o}\n")
            fh.write(f"{pool},{p.name}_RHS,/5Phos/{rhs_o}\n")

    # the FASTA is for BLAST, so it carries no /5Phos/ modification
    write_fasta(fa, [rec for p, lhs_o, rhs_o in rows
                     for rec in ((f"{p.name}_LHS", lhs_o), (f"{p.name}_RHS", rhs_o))])
    return tsv, csv, fa


def probes_main(argv) -> int:
    """`orthoseq probes INPUT` -- 10x Flex v2 probe pairs against given targets."""
    p = argparse.ArgumentParser(
        prog="orthoseq probes",
        description="Design 10x GEM-X Flex v2 custom probe pairs (25 nt LHS + "
        "25 nt RHS, reverse complementary to a 50 nt target site) against a set "
        "of sequences: FASTA, or one sequence per line (.txt/.csv/.tsv).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", help="target sequences to design probes against")
    p.add_argument("-o", "--out", default=None, metavar="STEM",
                   help="output stem; writes <stem>.tsv, <stem>_order.csv and "
                        "<stem>_probes.fa (default: <input> without extension)")
    p.add_argument("--workflow", choices=sorted(RHS_TAILS), default="multiplex",
                   help="RHS 3' tail: 'multiplex' = partial Constant Sequence "
                        "(CG000834/CG000835), 'singleplex' = partial Capture "
                        "Sequence 1 (CG000841, 4 Samples Kit)")
    p.add_argument("--pairs-per-target", type=int, default=3, metavar="N",
                   help="probe pairs to design per target; 10x recommends 3")
    p.add_argument("--min-gap", type=int, default=0, metavar="NT",
                   help="extra spacing required between the 50 nt sites of two "
                        "pairs on the same target (they are never allowed to overlap)")

    g = p.add_argument_group("per-half rules (CG000839)")
    g.add_argument("--gc-min", type=float, default=0.44)
    g.add_argument("--gc-max", type=float, default=0.72)
    g.add_argument("--max-homopolymer", type=int, default=4)
    g.add_argument("--max-tandem-copies", type=int, default=4)
    g.add_argument("--min-entropy", type=float, default=4.0,
                   help="min Shannon entropy over the 3-mers of a 25 nt half, in "
                        "bits; the ceiling for 25 nt is log2(23) = 4.52")
    g.add_argument("--no-require-tn", action="store_true",
                   help="accept sites without the recommended T at the ligation "
                        "junction (an A at offset 25 of the target site)")

    g = p.add_argument_group("off-target screening")
    g.add_argument("--cross-mismatch-min", type=int, default=5, metavar="MM",
                   help="flag a pair unless both halves are this many mismatches "
                        "from every unintended site in the input set; 0 disables")
    g.add_argument("--ref-db", action="append", default=[], metavar="DB",
                   help="BLAST nucleotide DB (e.g. a transcriptome) to report "
                        "each probe half's best hit against (repeatable). The "
                        "intended site is not subtracted, so a DB that contains "
                        "the targets themselves reports 0 mismatches")
    g.add_argument("--ref-fasta", action="append", default=[], metavar="FASTA",
                   help="reference FASTA to screen against (repeatable)")
    g.add_argument("--ref-word-size", type=int, default=7)
    g.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    g.add_argument("--tmpdir", default=None)
    p.add_argument("-q", "--quiet", action="store_true")
    a = p.parse_args(argv)
    verbose = not a.quiet

    targets = read_seqs(a.input)
    if not targets:
        sys.exit(f"[orthoseq] no sequences read from {a.input}")
    if a.pairs_per_target < 1:
        sys.exit("[orthoseq] --pairs-per-target must be at least 1")

    rules = ProbeRules(
        gc_min=a.gc_min,
        gc_max=a.gc_max,
        max_homopolymer=a.max_homopolymer,
        max_tandem_copies=a.max_tandem_copies,
        min_entropy=a.min_entropy,
        require_tn=not a.no_require_tn,
    )

    short = [n for n, s in targets if len(s) < PROBE_SITE]
    if short and verbose:
        print(f"[orthoseq] {len(short)} target(s) shorter than {PROBE_SITE} nt, "
              f"no probe possible: {', '.join(short[:5])}"
              f"{' ...' if len(short) > 5 else ''}", file=sys.stderr)

    pairs: list[ProbePair] = []
    empty: list[str] = []
    for name, seq in targets:
        picked = select_pairs(scan_sites(name, seq, rules),
                              a.pairs_per_target, a.min_gap)
        pairs.extend(picked)
        if not picked and len(seq) >= PROBE_SITE:
            empty.append(name)

    if not pairs:
        sys.exit("[orthoseq] no target yielded an acceptable probe pair; try "
                 "--no-require-tn or a wider --gc-min/--gc-max")

    cross_screen(pairs, targets, a.cross_mismatch_min)

    ref_hits: dict[str, tuple[str, int]] = {}
    stem = a.out or os.path.splitext(a.input)[0]
    tmpdir = tempfile.mkdtemp(prefix="orthoseq.", dir=a.tmpdir)
    try:
        if a.ref_db or a.ref_fasta:
            require_tools("blastn", *(["makeblastdb"] if a.ref_fasta else []))
            dbs = make_ref_dbs(a.ref_db, a.ref_fasta, tmpdir, verbose)
            query_fa = os.path.join(tmpdir, "probe_halves.fa")
            write_fasta(query_fa, [(f"{p.name}_{side}", half)
                                   for p in pairs
                                   for side, half in (("LHS", p.lhs), ("RHS", p.rhs))])
            for db in dbs:
                for h in blast_vs_db(query_fa, db, 10.0, a.ref_word_size,
                                     a.threads, verbose):
                    # mismatches the hybrid would carry: substitutions in the
                    # alignment, plus every probe base the alignment never covered
                    mm = h.mismatch + (PROBE_HALF - h.length)
                    prev = ref_hits.get(h.qseqid)
                    if prev is None or mm < prev[1]:
                        ref_hits[h.qseqid] = (os.path.basename(db), mm)
        out = write_probe_outputs(stem, pairs, a.workflow, a.cross_mismatch_min, ref_hits)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    if verbose:
        counts = Counter(p.target for p in pairs)
        designed = len(counts)
        full = sum(1 for c in counts.values() if c == a.pairs_per_target)
        print(f"[orthoseq] {len(pairs)} probe pairs for {designed}/{len(targets)} "
              f"targets ({full} with the full {a.pairs_per_target})", file=sys.stderr)
        if empty:
            print(f"[orthoseq] no acceptable site in {len(empty)} target(s): "
                  f"{', '.join(empty[:5])}{' ...' if len(empty) > 5 else ''}",
                  file=sys.stderr)
        flagged = [p.name for p in pairs
                   if p.cross_mm is not None and p.cross_mm < a.cross_mismatch_min]
        if flagged:
            print(f"[orthoseq] {len(flagged)} pair(s) within "
                  f"{a.cross_mismatch_min} mismatches of another target site: "
                  f"{', '.join(flagged[:5])}{' ...' if len(flagged) > 5 else ''}",
                  file=sys.stderr)
        print(f"[orthoseq] -> {', '.join(out)}", file=sys.stderr)
    return 0 if not empty and not short else 1


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="orthoseq",
        description="Generate n mutually dissimilar sequences of length l that "
        "avoid reference genomes and a user exclusion list.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("-n", "--num", type=int, required=True,
                   help="number of sequences in the final set; with --resume the "
                        "existing sequences count toward it")
    p.add_argument("-l", "--length", type=int, default=None,
                   help="sequence length (nt); read off the file when --resume is given")
    p.add_argument("-o", "--out", default=None,
                   help="output FASTA (default: orthoseq.fa, or the --resume file "
                        "itself, which is rewritten with its own sequences plus the new ones)")
    p.add_argument("--report", default=None, help="output TSV report (default: <out>.tsv)")
    p.add_argument("--resume", default=None, metavar="FILE",
                   help="pick up from an existing set: FASTA, or one sequence per "
                        "line (.txt/.csv/.tsv). Its sequences are kept as-is, count "
                        "toward -n, and constrain the new ones exactly as freshly "
                        "accepted sequences do. --length comes from the file")

    g = p.add_argument_group("references")
    g.add_argument(
        "--ref-db",
        action="append",
        default=[],
        metavar="DB",
        help="BLAST nucleotide DB to avoid (repeatable), e.g. human_genome",
    )
    g.add_argument(
        "--ref-fasta",
        action="append",
        default=[],
        metavar="FASTA",
        help="reference FASTA to avoid; a temporary BLAST DB is built (repeatable)",
    )
    g.add_argument(
        "--ref-evalue",
        type=float,
        default=1e-3,
        help="reject a candidate if it has any reference hit at or below this E-value",
    )
    g.add_argument("--ref-word-size", type=int, default=11, help="blastn word size vs references")

    g = p.add_argument_group("mutual dissimilarity")
    g.add_argument("--min-hamming", type=float, default=0.5, metavar="D",
                   help="minimum pairwise Hamming distance between accepted "
                        "sequences; below 1 it is a fraction of --length, 1 or "
                        "more an absolute base count, 0 disables the check")
    g.add_argument("--no-hamming-rc", action="store_true",
                   help="only compare sequences in the given orientation, rather "
                        "than also against each other's reverse complements")
    g.add_argument("--self-k", type=int, default=0, metavar="K",
                   help="also require that no two accepted sequences share a "
                        "canonical k-mer, catching similarity at an offset that "
                        "position-locked Hamming distance cannot see; 0 disables")
    g.add_argument("--self-blast", action="store_true",
                   help="also run an all-vs-all blastn and reject conflicting "
                        "pairs; catches gapped similarity, and is the only check "
                        "here that tolerates indels")
    g.add_argument("--self-bitscore", type=float, default=25.0,
                   help="with --self-blast, reject a pair whose HSP bitscore is >= this")
    g.add_argument("--self-word-size", type=int, default=7,
                   help="with --self-blast, blastn word size for the all-vs-all")

    g = p.add_argument_group("exclusion list")
    g.add_argument("--exclude", default=None, metavar="FILE",
                   help="sequences to avoid: FASTA, or one per line (.txt/.csv/.tsv). "
                        "Candidates identical to, contained in, or containing any "
                        "listed sequence are rejected")
    g.add_argument("--exclude-orientations", nargs="+", default=list(ORIENTATIONS),
                   choices=ORIENTATIONS, metavar="ORI",
                   help="orientations of each listed sequence to enforce, from: "
                        + " ".join(ORIENTATIONS))
    g.add_argument("--exclude-k", type=int, default=14,
                   help="also reject candidates sharing a canonical k-mer with the list; 0 disables")

    g = p.add_argument_group("composition")
    g.add_argument("--gc", type=float, default=0.5, help="target GC fraction")
    g.add_argument("--gc-tol", type=float, default=0.05, help="allowed GC deviation")
    g.add_argument("--max-homopolymer", type=int, default=4)
    g.add_argument("--max-tandem-copies", type=int, default=4,
                   help="max copies of a 2-3 nt tandem repeat unit")
    g.add_argument("--min-entropy", type=float, default=5.0,
                   help="min Shannon entropy over 3-mers, in bits (max 6)")
    g.add_argument("--avoid-motif", action="append", default=[], metavar="MOTIF",
                   help="forbid this motif and its reverse complement (repeatable)")
    g.add_argument("--gc-window", default=None, metavar="W[:MIN:MAX]",
                   help="also require every W-nt sub-window to have GC within "
                        "[MIN, MAX], defaulting to the --gc/--gc-tol window. "
                        "--gc-window 25:0.44:0.72 makes every sequence satisfy "
                        "the 10x Flex per-probe-half GC rule by construction")
    g.add_argument("--fix-base", action="append", default=[], metavar="POS:BASE",
                   help="hold one 0-based position at one base (repeatable), e.g. "
                        "--fix-base 25:A to give every sequence the A that the "
                        "10x Flex TN ligation junction needs at that offset")

    g = p.add_argument_group("run control")
    g.add_argument("--oversample", type=float, default=4.0,
                   help="candidates generated per still-needed sequence, per round")
    g.add_argument("--max-rounds", type=int, default=8)
    g.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    g.add_argument("--seed", type=int, default=None)
    g.add_argument("--tmpdir", default=None)
    g.add_argument("-q", "--quiet", action="store_true")
    return p.parse_args(argv)


def report_main(argv) -> int:
    """`orthoseq report INPUT` — stats for sequences you already have."""
    p = argparse.ArgumentParser(
        prog="orthoseq report",
        description="Report per-sequence stats (and the best reference hit) for an "
        "existing set of sequences: FASTA, or one sequence per line (.txt/.csv/.tsv).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", help="FASTA, or a file with one sequence per line")
    p.add_argument("-o", "--report", default=None, help="output TSV (default: <input>.tsv)")
    p.add_argument("--ref-db", action="append", default=[], metavar="DB",
                   help="BLAST nucleotide DB to screen against (repeatable)")
    p.add_argument("--ref-fasta", action="append", default=[], metavar="FASTA",
                   help="reference FASTA to screen against (repeatable)")
    p.add_argument("--ref-word-size", type=int, default=11)
    p.add_argument("--no-hamming-rc", action="store_true",
                   help="skip the reverse-complement nearest-neighbour columns")
    p.add_argument("--threads", type=int, default=os.cpu_count() or 4)
    p.add_argument("--tmpdir", default=None)
    p.add_argument("-q", "--quiet", action="store_true")
    a = p.parse_args(argv)
    verbose = not a.quiet

    records = read_seqs(a.input)
    if not records:
        sys.exit(f"[orthoseq] no sequences read from {a.input}")
    if a.ref_db or a.ref_fasta:
        require_tools("blastn", *(["makeblastdb"] if a.ref_fasta else []))

    report_path = a.report or (a.input + ".tsv")
    tmpdir = tempfile.mkdtemp(prefix="orthoseq.", dir=a.tmpdir)
    try:
        dbs = make_ref_dbs(a.ref_db, a.ref_fasta, tmpdir, verbose)
        write_report(report_path, records, dbs, tmpdir, a.ref_word_size,
                     a.threads, verbose, hamming_rc=not a.no_hamming_rc)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    if verbose:
        print(f"[orthoseq] {len(records)} sequences -> {report_path}", file=sys.stderr)
    return 0


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == "report":
        return report_main(argv[1:])
    if argv and argv[0] == "probes":
        return probes_main(argv[1:])

    a = parse_args(argv)
    verbose = not a.quiet
    rng = random.Random(a.seed)

    require_tools("blastn", "makeblastdb")

    seed_records: list[tuple[str, str]] = []
    if a.resume:
        seed_records = dedupe_names(read_seqs(a.resume))
        if not seed_records:
            sys.exit(f"[orthoseq] no sequences read from {a.resume}")
        seed_len = infer_length(seed_records, a.resume)
        if a.length is None:
            a.length = seed_len
        elif a.length != seed_len:
            sys.exit(f"[orthoseq] --length {a.length} disagrees with the {seed_len} nt "
                     f"sequences in {a.resume}")
        if verbose:
            print(f"[orthoseq] resuming from {a.resume}: {len(seed_records)} sequences "
                  f"of {a.length} nt, need {max(0, a.num - len(seed_records))} more",
                  file=sys.stderr)
    elif a.length is None:
        sys.exit("[orthoseq] --length is required unless --resume supplies it")
    if a.out is None:
        a.out = a.resume or "orthoseq.fa"

    if a.self_k > a.length:
        sys.exit("[orthoseq] --self-k must be <= --length")

    min_hamming = resolve_min_hamming(a.min_hamming, a.length)
    if min_hamming > a.length:
        sys.exit(f"[orthoseq] --min-hamming resolves to {min_hamming}, which exceeds "
                 f"--length {a.length}")
    if verbose and min_hamming:
        rc_note = "" if a.no_hamming_rc else ", also vs reverse complements"
        print(f"[orthoseq] min pairwise Hamming distance: "
              f"{min_hamming}/{a.length}{rc_note}", file=sys.stderr)

    # A resumed set is kept whole even if it was built under looser settings, so
    # say when it does not meet the threshold the new sequences are held to.
    if min_hamming and len(seed_records) > 1:
        tables = pairwise_hamming(seed_records, not a.no_hamming_rc)
        worst = min((d for t in tables for d, _ in t.values()), default=min_hamming)
        if worst < min_hamming:
            print(f"[orthoseq] warning: {a.resume} already contains a pair {worst} apart, "
                  f"below --min-hamming {min_hamming}; it is kept as given, and the new "
                  f"sequences are still held to the threshold", file=sys.stderr)

    win = parse_gc_window(a.gc_window, a.gc, a.gc_tol, a.length)
    if verbose and win:
        print(f"[orthoseq] per-window GC: every {win[0]} nt within "
              f"{win[1]:.3f}-{win[2]:.3f}", file=sys.stderr)

    fixed = parse_fixed_bases(a.fix_base, a.length)
    if verbose and fixed:
        print("[orthoseq] fixed bases: "
              + ", ".join(f"{p}={b}" for p, b in sorted(fixed.items())), file=sys.stderr)

    filters = LocalFilters(
        gc_target=a.gc,
        gc_tol=a.gc_tol,
        max_homopolymer=a.max_homopolymer,
        max_tandem_copies=a.max_tandem_copies,
        min_entropy=a.min_entropy,
        avoid_motifs=[m.upper() for m in a.avoid_motif],
        gc_window=win[0] if win else 0,
        gc_window_min=win[1] if win else 0.0,
        gc_window_max=win[2] if win else 1.0,
    )

    tmpdir = tempfile.mkdtemp(prefix="orthoseq.", dir=a.tmpdir)
    try:
        # reference FASTAs -> temporary BLAST DBs
        dbs = make_ref_dbs(a.ref_db, a.ref_fasta, tmpdir, verbose)
        if not dbs and verbose:
            print("[orthoseq] warning: no reference DB given; skipping genome screen",
                  file=sys.stderr)

        excl = None
        if a.exclude:
            recs = read_seqs(a.exclude)
            if not recs:
                sys.exit(f"[orthoseq] no sequences read from {a.exclude}")
            excl = ExclusionSet(recs, a.exclude_k if a.exclude_k > 0 else 10**9,
                                orientations=tuple(a.exclude_orientations))
            if a.exclude_k <= 0:
                excl.kmers = set()
            if verbose:
                print(f"[orthoseq] exclusion list: {len(recs)} sequences -> "
                      f"{len(excl.variants)} oriented variants "
                      f"({', '.join(excl.orientations)})", file=sys.stderr)

        accepted: list[tuple[str, str]] = list(seed_records)
        accepted_kmers: set[str] = set()
        accepted_guard = HammingGuard(a.length, min_hamming, not a.no_hamming_rc)
        seen: set[str] = set()
        stats: dict[str, int] = {}
        for _, s in seed_records:
            seen.add(s)
            accepted_guard.add(s)
            if a.self_k > 0:
                accepted_kmers |= canonical_kmers(s, a.self_k)

        for rnd in range(1, a.max_rounds + 1):
            need = a.num - len(accepted)
            if need <= 0:
                break
            batch = max(need * int(a.oversample), 32)
            if verbose:
                print(f"[orthoseq] round {rnd}: need {need}, generating {batch}",
                      file=sys.stderr)

            # --- generate + local filters + exclusion + intra-batch k-mer dedup
            cands: dict[str, str] = {}
            cand_kmers: dict[str, set[str]] = {}
            batch_kmers = set(accepted_kmers)
            # Seeded from the accepted set only, so candidates discarded later in
            # this round do not permanently reserve sequence space.
            batch_guard = accepted_guard.clone()
            tries = 0
            while len(cands) < batch and tries < batch * 200:
                tries += 1
                s = random_seq(a.length, a.gc, rng, fixed)
                if s in seen:
                    continue
                seen.add(s)
                reason = filters.check(s)
                if reason:
                    stats[reason.split(":")[0].split("=")[0]] = (
                        stats.get(reason.split(":")[0].split("=")[0], 0) + 1
                    )
                    continue
                if excl and (reason := excl.reject(s)):
                    stats[reason] = stats.get(reason, 0) + 1
                    continue
                if batch_guard.too_close(s) is not None:
                    stats["min_hamming"] = stats.get("min_hamming", 0) + 1
                    continue
                km = canonical_kmers(s, a.self_k) if a.self_k > 0 else set()
                if km and (km & batch_kmers):
                    stats["self_kmer"] = stats.get("self_kmer", 0) + 1
                    continue
                batch_kmers |= km
                batch_guard.add(s)
                cid = f"c{rnd}_{len(cands):05d}"
                cands[cid] = s
                cand_kmers[cid] = km

            if not cands:
                if verbose:
                    print("[orthoseq] round produced no candidates; loosening may be needed",
                          file=sys.stderr)
                continue

            cand_fa = os.path.join(tmpdir, f"cand_{rnd}.fa")
            write_fasta(cand_fa, list(cands.items()))

            # --- reference screen
            best_ref: dict[str, tuple[str, float, float]] = {}
            for db in dbs:
                for h in blast_vs_db(cand_fa, db, a.ref_evalue, a.ref_word_size,
                                     a.threads, verbose):
                    prev = best_ref.get(h.qseqid)
                    if prev is None or h.bitscore > prev[1]:
                        best_ref[h.qseqid] = (os.path.basename(db), h.bitscore, h.evalue)
            for cid in list(cands):
                if cid in best_ref:
                    del cands[cid]
                    stats["reference_hit"] = stats.get("reference_hit", 0) + 1
            if verbose:
                print(f"[orthoseq]   {len(cands)} candidates survive reference screen",
                      file=sys.stderr)
            if not cands:
                continue

            # --- mutual dissimilarity: Hamming and k-mer tests are already
            # applied during generation, so only the optional blastn runs here
            conflicts: dict[str, set[str]] = {c: set() for c in cands}
            drop_vs_accepted: set[str] = set()
            if a.self_blast:
                surv_fa = os.path.join(tmpdir, f"surv_{rnd}.fa")
                write_fasta(surv_fa, list(cands.items()))
                subj = list(cands.items()) + accepted
                subj_fa = os.path.join(tmpdir, f"subj_{rnd}.fa")
                write_fasta(subj_fa, subj)
                accepted_ids = {name for name, _ in accepted}
                for h in blast_vs_subject(surv_fa, subj_fa, a.self_word_size,
                                          a.threads, verbose):
                    if h.qseqid == h.sseqid or h.bitscore < a.self_bitscore:
                        continue
                    if h.sseqid in accepted_ids:
                        drop_vs_accepted.add(h.qseqid)
                    elif h.sseqid in conflicts:
                        conflicts[h.qseqid].add(h.sseqid)
                        conflicts[h.sseqid].add(h.qseqid)
                for cid in drop_vs_accepted:
                    cands.pop(cid, None)
                    conflicts.pop(cid, None)
                    stats["self_blast"] = stats.get("self_blast", 0) + 1
                for cid in conflicts:
                    conflicts[cid] &= set(cands)

            # --- greedy selection
            priority = {c: 0.0 for c in cands}
            picked = greedy_independent_set(
                list(cands), conflicts, a.num - len(accepted), priority
            )
            for cid in picked:
                accepted.append((cid, cands[cid]))
                accepted_kmers |= cand_kmers[cid]
                accepted_guard.add(cands[cid])
            if verbose:
                print(f"[orthoseq]   accepted {len(picked)} "
                      f"(total {len(accepted)}/{a.num})", file=sys.stderr)

        if len(accepted) < a.num:
            print(f"[orthoseq] WARNING: only {len(accepted)}/{a.num} sequences found "
                  f"after {a.max_rounds} rounds. Loosen --min-hamming, --ref-evalue, "
                  f"--gc-tol, or raise --max-rounds/--oversample.", file=sys.stderr)

        # Resumed sequences keep the names they came with; only the new ones are
        # numbered, and around whatever names the old ones already occupy.
        final = list(seed_records)
        final += assign_names(accepted[len(seed_records):],
                              taken={n for n, _ in final})
        write_fasta(a.out, final)

        # --- report: best surviving reference hit for each accepted sequence
        report_path = a.report or (a.out + ".tsv")
        write_report(report_path, final, dbs, tmpdir, a.ref_word_size, a.threads,
                     verbose, hamming_rc=not a.no_hamming_rc)

        if verbose:
            print(f"[orthoseq] wrote {len(final)} sequences -> {a.out}", file=sys.stderr)
            print(f"[orthoseq] report -> {report_path}", file=sys.stderr)
            if stats:
                summary = ", ".join(f"{k}={v}" for k, v in sorted(stats.items()))
                print(f"[orthoseq] rejections: {summary}", file=sys.stderr)
        return 0 if len(accepted) >= a.num else 1
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
