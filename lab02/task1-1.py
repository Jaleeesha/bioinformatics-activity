"""
Task 1-1 - a closer look at the sequence: reading frames, genes and mutations.

  Option 1  Read-through   translate from the first AUG, mark stops with *
  Option 2  Fragments      cut the read-through at every stop
  Option 3  All genes      every AUG ... stop stretch (a new gene may start
                           after each stop)
  Option 4  Stop count     how many UAA / UAG / UGA are in frame
  Option 5  3 frames       translate the 3 reading frames, find the longest ORF
  Option 6  Mutation       compare with a mutated sequence (optional)

Run:
  python task1-1.py
  python task1-1.py ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
"""

import sys

from genetic_code import first_protein, get_sequence, to_mrna, translate


def heading(text):
    print(f"\n--- {text} ---")


def all_genes(mrna):
    """Every protein that starts at an AUG and ends at the next in-frame stop."""
    genes = []
    pos = mrna.find("AUG")
    while pos != -1:
        end = pos
        protein = ""
        while end + 3 <= len(mrna):
            aa = translate(mrna[end:end + 3])
            if aa == "*":
                break
            protein += aa
            end += 3
        genes.append(protein)
        pos = mrna.find("AUG", end + 3)
    return genes


def longest_orf(letters):
    best = ""
    for piece in letters.split("*"):
        m = piece.find("M")
        if m != -1 and len(piece) - m > len(best):
            best = piece[m:]
    return best


mrna = to_mrna(get_sequence(sys.argv[1:]))
start = mrna.find("AUG")
print("mRNA:", mrna)

if start == -1:
    print("No start codon (AUG) found.")
    sys.exit()

reading = translate(mrna, start)

heading("Option 1: read-through (stops shown as *)")
print(reading)

heading("Option 2: fragments between stops")
for n, piece in enumerate([p for p in reading.split("*") if p], 1):
    print(f"Fragment {n}: {piece}")

heading("Option 3: all genes (AUG ... stop)")
for n, protein in enumerate(all_genes(mrna), 1):
    print(f"Gene {n}: {protein}  ({len(protein)} aa)")

heading("Option 4: stop codons in frame")
stops = {"UAA": 0, "UAG": 0, "UGA": 0}
for i in range(start, len(mrna) - 2, 3):
    if mrna[i:i + 3] in stops:
        stops[mrna[i:i + 3]] += 1
for codon, n in stops.items():
    print(f"{codon}: {n}")
print("Total:", sum(stops.values()))

heading("Option 5: the 3 reading frames")
best, best_frame = "", 0
for frame in range(3):
    letters = translate(mrna, frame)
    orf = longest_orf(letters)
    print(f"Frame {frame + 1}: {letters}   longest ORF: {orf or '-'}")
    if len(orf) > len(best):
        best, best_frame = orf, frame + 1
print(f"Most likely gene: frame {best_frame}, {best} ({len(best)} aa)" if best else "No ORF found.")

heading("Option 6: mutation check")
try:
    other = input("Mutated sequence (Enter to skip): ")
except EOFError:
    other = ""
if not other.strip():
    print("Skipped.")
else:
    normal, mutated = first_protein(mrna), first_protein(to_mrna(other))
    print(f"Normal : {normal} ({len(normal)} aa)")
    print(f"Mutated: {mutated} ({len(mutated)} aa)")
    if normal == mutated:
        print("No change in the protein (silent mutation).")
    elif len(mutated) < len(normal) and normal.startswith(mutated):
        print(f"Nonsense mutation: early stop after amino acid {len(mutated)}.")
    else:
        diff = next((k for k in range(min(len(normal), len(mutated))) if normal[k] != mutated[k]),
                    min(len(normal), len(mutated)))
        print(f"The protein changes from amino acid {diff + 1}.")
