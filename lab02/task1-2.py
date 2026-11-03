"""
Task 1-2 - translation table: DNA codon, mRNA codon and amino acid side by side.

For every codon from the first AUG to the first stop it prints the DNA codon,
the mRNA codon, the amino acid (3-letter and 1-letter), then the full protein.

Run:
  python task1-2.py
  python task1-2.py ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
"""

import sys

from genetic_code import CODON_TABLE, THREE_LETTER, get_sequence, to_mrna

dna = get_sequence(sys.argv[1:])
mrna = to_mrna(dna)
start = mrna.find("AUG")

if start == -1:
    print("No start codon (AUG) found.")
    sys.exit()

print(f"{'DNA':<6}{'mRNA':<6}{'Amino acid':<12}{'Letter'}")
print("-" * 30)

names, letters = [], ""
for i in range(start, len(mrna) - 2, 3):
    letter = CODON_TABLE.get(mrna[i:i + 3], "?")
    name = THREE_LETTER.get(letter, "?")
    print(f"{dna[i:i + 3]:<6}{mrna[i:i + 3]:<6}{name:<12}{letter}")
    if letter == "*":
        break
    names.append(name)
    letters += letter

print("\nProtein  :", "-".join(names))
print("1-letter :", letters)
