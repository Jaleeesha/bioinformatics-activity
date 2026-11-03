"""
Task 1-1 - the same count as task1.py, but on a random sequence.

Builds a random DNA sequence of 20 letters, then counts A, C, G and T
and prints the count and frequency of each.

Run:
  python task1-1.py          (different sequence every time)
  python task1-1.py 7        (seed 7: same sequence every time)
"""

import random
import sys

BASES = "ACGT"
LENGTH = 20

if len(sys.argv) > 1:
    random.seed(int(sys.argv[1]))

S = "".join(random.choices(BASES, k=LENGTH))

print("S =", S)
print("Length:", len(S))
for base in BASES:
    count = S.count(base)
    freq = count / len(S)
    print(f"{base}: count={count:>2}  freq={freq:.2f} ({freq * 100:.0f}%)")
