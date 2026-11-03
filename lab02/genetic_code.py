"""
Shared helpers for Lab 02: the genetic code and small translation functions.
Imported by task1.py, task1-1.py and task1-2.py.
"""

BASES = "UCAG"
# Amino acid for every codon, in the order UUU, UUC, UUA, UUG, UCU ... GGG
# ("*" = stop codon).
AMINO_ACIDS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"

CODON_TABLE = {
    a + b + c: AMINO_ACIDS[i]
    for i, (a, b, c) in enumerate((x, y, z) for x in BASES for y in BASES for z in BASES)
}

THREE_LETTER = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys",
    "Q": "Gln", "E": "Glu", "G": "Gly", "H": "His", "I": "Ile",
    "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe", "P": "Pro",
    "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val",
    "*": "Stop",
}


def clean(seq):
    """Upper-case the text and drop spaces and line breaks."""
    return "".join(seq.split()).upper()


def to_mrna(seq):
    """Transcription: DNA -> mRNA (every T becomes U)."""
    return clean(seq).replace("T", "U")


def translate(mrna, start=0):
    """Read codons from `start` and return the one-letter protein string
    (stop codons are kept as '*', unknown codons as '?')."""
    protein = ""
    for i in range(start, len(mrna) - 2, 3):
        protein += CODON_TABLE.get(mrna[i:i + 3], "?")
    return protein


def first_protein(mrna):
    """Protein from the first AUG up to the first stop ('' if no AUG)."""
    start = mrna.find("AUG")
    if start == -1:
        return ""
    return translate(mrna, start).split("*")[0]


def get_sequence(args):
    """Sequence from the command line, otherwise ask the user."""
    text = " ".join(args) if args else input("Enter DNA or mRNA sequence: ")
    return clean(text)
