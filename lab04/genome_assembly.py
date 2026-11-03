"""
Lab 04 - Genome assembly from random reads.

Task: take a genome, cut 20,000 random pieces (reads) of 100-150 nucleotides
from it, hope the pieces cover the whole genome, then forget the genome and
rebuild it only from the pieces by looking at what neighbouring reads have in
common.

How it works:
  1. Reads: pick a random length (100-150) and a random start, copy that piece.
     Repeat 20,000 times. This is what a sequencing machine does.
  2. Coverage: count how many reads cover each position of the genome. A
     position with coverage 0 was never read and can not be rebuilt.
  3. Overlap: two reads that sit next to each other in the genome share a
     stretch of letters: the END of the first equals the START of the second.

         contig: ...ACGTTAGGCATCCGA
         read  :          GGCATCCGATTACG     overlap = GGCATCCGA
         result: ...ACGTTAGGCATCCGATTACG     only the new letters are added

  4. Greedy assembly: start from one read (the contig). Look for reads whose
     start equals the end of the contig (overlap of at least 30 letters), glue
     the one that adds the most new letters, and repeat until nothing fits.
     Then do the same to the left (by reversing all the text).
  5. To find matching reads fast, reads are stored in a dictionary under their
     first 30 letters, so we never compare every read with every other read.
  6. If part of the genome was never read, the contig gets stuck at the gap and
     a new contig is started from a read that is not used yet.
  7. The longest contig is compared with the original genome.

Run:
  python genome_assembly.py
  python genome_assembly.py genomes/NC_005816.fasta
  python genome_assembly.py --reads 60        (too few reads: gaps and several contigs)
Saves the coverage chart as coverage.png (add --no-plot to skip it).
"""

import argparse
import gzip
import random
import time

DEFAULT_GENOME = "genomes/GCF_000907355.1_ViralProj202889_genomic.fna.gz"


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    name, parts = "", []
    with opener(path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = name or line[1:]
            else:
                parts.append(line.upper())
    return name, "".join(parts)


def take_reads(genome, count, shortest, longest, rng):
    """Random windows of random length. Also returns where each one started."""
    reads, places = [], []
    for _ in range(count):
        length = rng.randint(shortest, longest)
        start = rng.randint(0, len(genome) - length)
        reads.append(genome[start:start + length])
        places.append((start, length))
    return reads, places


def coverage_per_position(genome_length, places):
    """How many reads cover each position (difference array, one pass)."""
    change = [0] * (genome_length + 1)
    for start, length in places:
        change[start] += 1
        change[start + length] -= 1
    cover, running = [], 0
    for i in range(genome_length):
        running += change[i]
        cover.append(running)
    return cover


def index_by_start(reads, size):
    index = {}
    for read in reads:
        index.setdefault(read[:size], []).append(read)
    return index


def extend_right(contig, index, size, longest):
    """Glue reads onto the right end of the contig until none fits."""
    steps = 0
    while True:
        best_read, best_overlap, best_gain = None, 0, 0
        # the overlap can start at most `longest` letters before the end
        for p in range(max(0, len(contig) - longest + 1), len(contig) - size + 1):
            tail = contig[p:]
            for read in index.get(tail[:size], ()):
                gain = len(read) - len(tail)
                if gain > best_gain and read.startswith(tail):
                    best_read, best_overlap, best_gain = read, len(tail), gain
        if best_read is None:
            return contig, steps
        contig += best_read[best_overlap:]
        steps += 1


def build_contig(seed, right_index, left_index, size, longest):
    contig, steps_right = extend_right(seed, right_index, size, longest)
    # the left side is the right side of the reversed text
    reversed_contig, steps_left = extend_right(contig[::-1], left_index, size, longest)
    return reversed_contig[::-1], steps_right + steps_left


def assemble(reads, size, longest):
    unique = sorted(set(reads), key=len, reverse=True)
    right_index = index_by_start(unique, size)
    left_index = index_by_start([r[::-1] for r in unique], size)

    contigs, total_steps = [], 0
    unused = unique
    while unused:
        contig, steps = build_contig(unused[0], right_index, left_index, size, longest)
        contigs.append(contig)
        total_steps += steps
        unused = [r for r in unused if r not in contig]
    contigs.sort(key=len, reverse=True)
    return contigs, len(unique), total_steps


def compare(original, rebuilt):
    """Percent of the original covered by exactly matching letters."""
    if rebuilt == original:
        return 100.0
    at = original.find(rebuilt[:30])
    if at == -1:
        return 0.0
    same = sum(1 for a, b in zip(original[at:], rebuilt) if a == b)
    return same / len(original) * 100


def draw_coverage(cover, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.figure(figsize=(11, 3.5))
    plt.plot(cover, color="steelblue", linewidth=0.8)
    plt.xlabel("Position in the genome")
    plt.ylabel("Reads covering it")
    plt.title("Coverage along the genome")
    plt.tight_layout()
    plt.savefig(path, dpi=120)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Cut a genome into random reads and rebuild it.")
    parser.add_argument("fasta", nargs="?", default=DEFAULT_GENOME)
    parser.add_argument("--reads", type=int, default=20000, help="number of reads (default 20000)")
    parser.add_argument("--min", dest="shortest", type=int, default=100, help="shortest read")
    parser.add_argument("--max", dest="longest", type=int, default=150, help="longest read")
    parser.add_argument("--overlap", type=int, default=30, help="shortest overlap to accept")
    parser.add_argument("--seed", type=int, default=42, help="random seed (same seed = same reads)")
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    name, genome = read_fasta(args.fasta)
    rng = random.Random(args.seed)

    print("Genome :", name)
    print("Length :", f"{len(genome):,} nt")

    started = time.time()
    reads, places = take_reads(genome, args.reads, args.shortest, args.longest, rng)
    cover = coverage_per_position(len(genome), places)
    average = sum(len(r) for r in reads) / len(reads)

    print("\n1. Reads")
    print(f"   Taken            : {args.reads:,} reads of {args.shortest}-{args.longest} nt (average {average:.1f})")
    print(f"   Expected coverage: {args.reads * average / len(genome):.0f}x  (reads x average length / genome length)")
    print(f"   Real coverage    : min {min(cover)}x, max {max(cover)}x, mean {sum(cover) / len(cover):.0f}x")
    print(f"   Never read       : {cover.count(0)} positions")

    contigs, unique_count, steps = assemble(reads, args.overlap, args.longest)
    rebuilt = contigs[0]

    print("\n2. Assembly")
    print(f"   Different reads  : {unique_count:,} (identical ones merged)")
    print(f"   Contigs          : {len(contigs)}" + ("  (one complete piece)" if len(contigs) == 1 else "  (gaps: some parts were never read)"))
    if len(contigs) > 1:
        print("   Contig lengths   :", ", ".join(f"{len(c):,}" for c in contigs[:10]))
    print(f"   Glue steps       : {steps}")
    print(f"   Time             : {time.time() - started:.2f} s")

    print("\n3. Result")
    print(f"   Original length  : {len(genome):,} nt")
    print(f"   Rebuilt length   : {len(rebuilt):,} nt" + (" (longest contig)" if len(contigs) > 1 else ""))
    if rebuilt == genome:
        print("   Match            : 100% identical, the genome was rebuilt perfectly")
    else:
        print(f"   Match            : {compare(genome, rebuilt):.2f}% of the original")
    print("   Original start   :", genome[:50])
    print("   Rebuilt start    :", rebuilt[:50])

    if not args.no_plot:
        draw_coverage(cover, "coverage.png")
        print("\nCoverage chart saved: coverage.png")


if __name__ == "__main__":
    main()
