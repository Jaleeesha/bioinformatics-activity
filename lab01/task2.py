"""
Task 2 - count 2-letter and 3-letter combinations in a DNA sequence.

How it works:
  1. Make a random sequence S of 20 letters (A, C, G, T).
  2. List every possible combination first: 4^2 = 16 pairs, 4^3 = 64 triples,
     each starting at 0, so combinations that never appear still show up.
  3. Slide a window of k letters along S one step at a time and add 1 to
     the counter of the combination inside the window.
     A sequence of length n gives n - k + 1 windows.
  4. frequency = count / number of windows.
  5. Draw one bar chart for the pairs and one for the triples.

Run:
  python task2.py            (different sequence every time)
  python task2.py 7          (seed 7: same sequence every time)
Saves the chart as task2_chart.png and also shows it on screen.
"""

import random
import sys
from itertools import product

import matplotlib.pyplot as plt

BASES = "ACGT"
LENGTH = 20


def count_combinations(seq, k):
    counts = {"".join(p): 0 for p in product(BASES, repeat=k)}
    windows = len(seq) - k + 1
    for start in range(windows):
        counts[seq[start:start + k]] += 1
    return counts, windows


def print_result(title, counts, windows):
    print(f"\n{title} ({len(counts)} combinations, {windows} windows)")
    for combo, count in counts.items():
        print(f"  {combo}: count={count}  freq={count / windows:.3f}")


if len(sys.argv) > 1:
    random.seed(int(sys.argv[1]))

S = "".join(random.choices(BASES, k=LENGTH))
print("S =", S)

pairs, pair_windows = count_combinations(S, 2)
triples, triple_windows = count_combinations(S, 3)

print_result("2-letter combinations", pairs, pair_windows)
print_result("3-letter combinations", triples, triple_windows)

fig, (top, bottom) = plt.subplots(2, 1, figsize=(10, 7))

top.bar(pairs.keys(), [c / pair_windows for c in pairs.values()], color="steelblue")
top.set_title(f"2-letter frequencies   S = {S}")
top.set_ylabel("Frequency")

bottom.bar(triples.keys(), [c / triple_windows for c in triples.values()], color="seagreen")
bottom.set_title("3-letter frequencies")
bottom.set_ylabel("Frequency")
bottom.tick_params(axis="x", rotation=90, labelsize=7)

plt.tight_layout()
plt.savefig("task2_chart.png", dpi=120)
plt.show()
