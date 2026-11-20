# Lab 5 - Detecting repetitions in a DNA sequence (NCBI: M77789.2, pUC19, 2686 bp)
#
# A repetition = a unit of 3-6 letters written several times in a row,
# with a total length between 6 and 10 letters. Examples:
#   unit "ACG"   (3 letters) x2 -> ACGACG       (6 letters)
#   unit "ACG"   (3 letters) x3 -> ACGACGACG    (9 letters)
#   unit "ACGT"  (4 letters) x2 -> ACGTACGT     (8 letters)
#   unit "ACGTA" (5 letters) x2 -> ACGTAACGTA   (10 letters)
# The units are produced by brute force: every combination of A, C, G, T with 3 to 6 letters.

import itertools

# Reading the DNA sequence downloaded from NCBI
S = ""
with open("sequence.fasta") as file:
    for line in file:
        if not line.startswith(">"):
            S += line.strip()

print("Sequence length:", len(S))

# Brute force engine: try every unit of 3 to 6 letters
found = []
for k in range(3, 7):
    for combo in itertools.product("ACGT", repeat=k):
        unit = "".join(combo)
        times = 2
        while k * times <= 10:
            if k * times >= 6:
                repetition = unit * times
                # search the repetition at every position of the sequence
                positions = []
                for i in range(len(S) - len(repetition) + 1):
                    if S[i:i + len(repetition)] == repetition:
                        positions.append(i + 1)  # positions start from 1
                if positions:
                    found.append((repetition, unit, times, positions))
            times += 1

print("Repetitions found:", len(found), "\n")
for repetition, unit, times, positions in found:
    print(repetition, "=", unit, "x", times, "| found", len(positions), "times at positions", positions)
