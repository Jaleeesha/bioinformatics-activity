"""
Lab 01 - DNA sequence analysis

1. Generate a random DNA sequence S.
2. Count the four bases A, C, G, T and their frequencies.
3. Count all 2-letter (16) and 3-letter (64) combinations and their frequencies.
4. Show the results as bar charts.

Run:
  python lab01.py            (random sequence of 20 letters)
  python lab01.py 100        (random sequence of 100 letters)
  python lab01.py ATGCGATAC  (your own sequence)
"""

import random
import sys
from itertools import product

import matplotlib.pyplot as plt

ALPHABET = "ACGT"


def make_sequence(length):
    return "".join(random.choice(ALPHABET) for _ in range(length))


def count_kmers(seq, k):
    """Slide a window of k letters over seq and count every combination."""
    counts = {"".join(p): 0 for p in product(ALPHABET, repeat=k)}
    for i in range(len(seq) - k + 1):
        counts[seq[i:i + k]] += 1
    return counts


def frequencies(counts):
    total = sum(counts.values())
    return {key: (c / total if total else 0) for key, c in counts.items()}


def print_table(title, counts, freqs):
    print(f"\n{title}")
    for key in counts:
        print(f"  {key}: count={counts[key]:>2}  freq={freqs[key]:.3f} ({freqs[key] * 100:.1f}%)")


def main():
    arg = sys.argv[1].upper() if len(sys.argv) > 1 else "20"
    seq = make_sequence(int(arg)) if arg.isdigit() else arg
    if not set(seq) <= set(ALPHABET) or len(seq) < 3:
        print("Sequence must have at least 3 letters, all from A, C, G, T.")
        return

    print("S =", seq)
    print("Length:", len(seq))

    results = {}
    for k, name in ((1, "Single bases"), (2, "2-letter combinations"), (3, "3-letter combinations")):
        counts = count_kmers(seq, k)
        freqs = frequencies(counts)
        results[k] = (name, freqs)
        print_table(name, counts, freqs)

    fig, axes = plt.subplots(3, 1, figsize=(11, 9))
    colors = {1: "tab:orange", 2: "steelblue", 3: "seagreen"}
    for ax, (k, (name, freqs)) in zip(axes, results.items()):
        ax.bar(freqs.keys(), freqs.values(), color=colors[k])
        ax.set_title(name)
        ax.set_ylabel("Frequency")
        if k == 3:
            ax.tick_params(axis="x", rotation=90, labelsize=7)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
