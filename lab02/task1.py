"""
Task 1 - translate a DNA / mRNA sequence into a protein.

Steps:
  1. Transcription: DNA -> mRNA (T becomes U).
  2. Find the first start codon AUG.
  3. Read the mRNA three letters (one codon) at a time from there.
  4. Look every codon up in the genetic code and print its amino acid.
  5. Stop at the first stop codon (UAA, UAG or UGA).

Run:
  python task1.py
  python task1.py ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
"""

import sys

from genetic_code import CODON_TABLE, THREE_LETTER, get_sequence, to_mrna

dna = get_sequence(sys.argv[1:])
mrna = to_mrna(dna)
start = mrna.find("AUG")

print("DNA :", dna)
print("mRNA:", mrna)

if start == -1:
    print("No start codon (AUG) found.")
    sys.exit()

print(f"Start codon AUG found at position {start + 1}\n")

protein = []
for i in range(start, len(mrna) - 2, 3):
    codon = mrna[i:i + 3]
    amino_acid = THREE_LETTER.get(CODON_TABLE.get(codon, "?"), "?")
    print(f"{codon} -> {amino_acid}")
    if amino_acid == "Stop":
        break
    protein.append(amino_acid)
else:
    print("(no stop codon: the sequence ended first)")

print("\nProtein:", "-".join(protein))
print("Length :", len(protein), "amino acids")
