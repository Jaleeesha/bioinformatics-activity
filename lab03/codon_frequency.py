"""
Lab 03 - codon frequency of a genome.

Reads a genome from a FASTA file (plain or .gz), counts how often each of the
64 codons appears and shows the most frequent ones.

How it works:
  1. Read the FASTA file: the line starting with ">" is the name, all other
     lines are the sequence and are joined into one string.
  2. Make a counter for all 64 codons (4 x 4 x 4), each starting at 0.
  3. Slide a window of 3 letters along the genome one step at a time and add 1
     to the codon inside the window. A genome of n letters has n - 2 windows.
     (With --frame 1, 2 or 3 only the codons of that one reading frame are
     counted: the window jumps 3 letters at a time.)
  4. frequency = count / number of codons counted.
  5. Sort from the most to the least frequent, print the top ones and save a
     bar chart of all 64 codons.

Run:
  python codon_frequency.py genomes/NC_005816.fasta
  python codon_frequency.py genomes/GCF_000907355.1_ViralProj202889_genomic.fna.gz --top 5
  python codon_frequency.py genomes/NC_005816.fasta --frame 1
Saves the chart as <file name>_codons.png (add --no-plot to skip it).
"""

import argparse
import gzip
import os
from itertools import product

AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = ["".join(c) for c in product("TCAG", repeat=3)]
AMINO_ACID_OF = dict(zip(CODONS, AMINO_ACIDS))


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


def count_codons(genome, frame=0):
    """frame = 0: every window; frame = 1, 2, 3: only that reading frame."""
    counts = {codon: 0 for codon in CODONS}
    step = 1 if frame == 0 else 3
    first = 0 if frame == 0 else frame - 1
    for i in range(first, len(genome) - 2, step):
        codon = genome[i:i + 3]
        if codon in counts:  # skips windows with N or other unknown letters
            counts[codon] += 1
    return counts


def gc_content(genome):
    return (genome.count("G") + genome.count("C")) / len(genome) * 100


def draw_chart(ranked, total, title, out_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.figure(figsize=(13, 5))
    plt.bar([c for c, _ in ranked], [n / total * 100 for _, n in ranked], color="seagreen")
    plt.xticks(rotation=90, fontsize=7)
    plt.ylabel("Frequency (%)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Codon frequencies of a genome (FASTA).")
    parser.add_argument("fasta", help="FASTA file, plain or .gz")
    parser.add_argument("--top", type=int, default=10, help="how many codons to list (default 10)")
    parser.add_argument("--frame", type=int, choices=[0, 1, 2, 3], default=0,
                        help="0 = all windows (default), 1-3 = one reading frame")
    parser.add_argument("--no-plot", action="store_true", help="do not save the chart")
    args = parser.parse_args()

    name, genome = read_fasta(args.fasta)
    if len(genome) < 3:
        print("The sequence is too short to hold a codon.")
        return

    counts = count_codons(genome, args.frame)
    total = sum(counts.values())
    ranked = sorted(counts.items(), key=lambda item: item[1], reverse=True)

    print("Genome :", name)
    print("Length :", f"{len(genome):,} nt")
    print("GC     :", f"{gc_content(genome):.2f}%")
    mode = "all windows" if args.frame == 0 else f"reading frame {args.frame}"
    print("Codons :", f"{total:,} counted ({mode})")
    if total < len(genome) - 2 and args.frame == 0:
        print("Skipped:", len(genome) - 2 - total, "windows with letters other than A, C, G, T")

    print(f"\nTop {args.top} codons")
    print(f"{'Rank':<6}{'Codon':<7}{'AA':<4}{'Count':>6}{'Freq %':>9}  Bar")
    print("-" * 52)
    biggest = ranked[0][1]
    for rank, (codon, n) in enumerate(ranked[:args.top], 1):
        bar = "#" * round(n / biggest * 25)
        print(f"{rank:<6}{codon:<7}{AMINO_ACID_OF[codon]:<4}{n:>6}{n / total * 100:>8.2f}%  {bar}")

    never = [codon for codon, n in counts.items() if n == 0]
    print(f"\nCodons that never appear: {len(never)}" + (f" ({', '.join(never)})" if never else ""))

    if not args.no_plot:
        base = os.path.basename(args.fasta).split(".")[0]
        out = f"{base}_codons.png"
        draw_chart(ranked, total, f"Codon frequencies - {name[:70]}", out)
        print("Chart saved:", out)


if __name__ == "__main__":
    main()
