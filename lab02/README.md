# Lab 02 - DNA to Protein Translation

Turning a DNA sequence into mRNA and then into a protein, codon by codon.

| File | What it does |
|---|---|
| `genetic_code.py` | The genetic code (64 codons) and small helper functions used by the tasks |
| `task1.py` | DNA -> mRNA -> protein: finds the first AUG and translates until a stop codon |
| `task1-1.py` | Extra analysis: read-through, fragments, all genes, stop codons, the 3 reading frames, mutation check |
| `task1-2.py` | Translation table with DNA codon, mRNA codon and amino acid in columns |

## How it works

- **Transcription**: DNA is copied into mRNA by changing every T into U.
- **Start codon**: reading starts at the first `AUG` (Met).
- **Codons**: the mRNA is read three letters at a time. 4^3 = 64 codons code for 20 amino acids and 3 stops, so several codons give the same amino acid.
- **Stop codons**: `UAA`, `UAG` and `UGA` end the protein.
- **Reading frames**: starting at position 1, 2 or 3 gives three different proteins; the real gene is usually the frame with the longest stretch from `M` to a stop (the longest ORF).
- **Nonsense mutation**: a change that creates an early stop codon and makes the protein shorter.

## Run

Give the sequence as an argument, or run without one and type it in.

```
python task1.py ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
python task1-1.py ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
python task1-2.py ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
```

## Output

Sequence used: `ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG`

### task1.py

```
DNA : ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG
mRNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG
Start codon AUG found at position 1

AUG -> Met
GCC -> Ala
AUU -> Ile
GUA -> Val
AUG -> Met
GGC -> Gly
CGC -> Arg
UGA -> Stop

Protein: Met-Ala-Ile-Val-Met-Gly-Arg
Length : 7 amino acids
```

### task1-2.py

```
DNA   mRNA  Amino acid  Letter
------------------------------
ATG   AUG   Met         M
GCC   GCC   Ala         A
ATT   AUU   Ile         I
GTA   GUA   Val         V
ATG   AUG   Met         M
GGC   GGC   Gly         G
CGC   CGC   Arg         R
TGA   UGA   Stop        *

Protein  : Met-Ala-Ile-Val-Met-Gly-Arg
1-letter : MAIVMGR
```

### task1-1.py

For the mutation check, `ATGGCCATTTAATGGGCCGCTGAAAGGGTGCCCGATAG` was typed in (the codon GTA became TAA).

```
mRNA: AUGGCCAUUGUAAUGGGCCGCUGAAAGGGUGCCCGAUAG

--- Option 1: read-through (stops shown as *) ---
MAIVMGR*KGAR*

--- Option 2: fragments between stops ---
Fragment 1: MAIVMGR
Fragment 2: KGAR

--- Option 3: all genes (AUG ... stop) ---
Gene 1: MAIVMGR  (7 aa)

--- Option 4: stop codons in frame ---
UAA: 0
UAG: 1
UGA: 1
Total: 2

--- Option 5: the 3 reading frames ---
Frame 1: MAIVMGR*KGAR*   longest ORF: MAIVMGR
Frame 2: WPL*WAAERVPD   longest ORF: -
Frame 3: GHCNGPLKGCPI   longest ORF: -
Most likely gene: frame 1, MAIVMGR (7 aa)

--- Option 6: mutation check ---
Mutated sequence (Enter to skip): ATGGCCATTTAATGGGCCGCTGAAAGGGTGCCCGATAG
Normal : MAIVMGR (7 aa)
Mutated: MAI (3 aa)
Nonsense mutation: early stop after amino acid 3.
```
