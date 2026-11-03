"""
Task 1 - count the bases of a given DNA sequence.

Takes a fixed sequence S, counts how many times A, C, G and T appear,
and prints the count and the frequency (count / length of S) of each.

Run:
  python task1.py
"""

S = "ATGCGATACGCTTAGGCTAA"
BASES = "ACGT"

counts = {}
for base in BASES:
    counts[base] = S.count(base)

print("S =", S)
print("Length:", len(S))
for base in BASES:
    freq = counts[base] / len(S)
    print(f"{base}: count={counts[base]:>2}  freq={freq:.2f} ({freq * 100:.0f}%)")

# The four frequencies must add up to 100%
print("Total:", sum(counts.values()), "letters =", f"{sum(counts.values()) / len(S) * 100:.0f}%")
