# Lab 6 - Detecting transposable elements in an artificial DNA sequence
#
# Transposon structure:
#   TSD | IR ...middle... reverse_complement(IR) | TSD
#   IR  = inverted repeat (the right end is the left end read backwards with A<->T, C<->G)
#   TSD = target site duplication (the same short piece copied on both sides)

import random


def random_dna(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def reverse_complement(s):
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(pairs[letter] for letter in reversed(s))


IR_LEN = 10

# Every transposon has its own inverted repeat
IR = {name: random_dna(IR_LEN) for name in ["T1", "T2", "T3", "T4"]}

# Building the artificial sequence and remembering the real positions (start, end)
S = ""
real = {}

S += random_dna(30)

# T1 and T2 overlap: T2 landed on the end of T1 and overwrote its last 15 letters
# (so the right inverted repeat of T1 is destroyed)
T1 = IR["T1"] + random_dna(30) + reverse_complement(IR["T1"])  # 50 letters
OVERWRITTEN = 15
tsd1 = random_dna(5)
S += tsd1
start1 = len(S) + 1
S += T1[:-OVERWRITTEN]  # only the first part of T1 survives
real["T1"] = (start1, start1 + len(T1) - 1)  # where T1 would have ended
start2 = len(S) + 1
S += IR["T2"] + random_dna(25) + reverse_complement(IR["T2"])
real["T2"] = (start2, len(S))
S += random_dna(5)

S += random_dna(30)

# T3 with T4 nested inside it (T4 jumped into the middle of T3)
tsd3 = random_dna(5)
tsd4 = random_dna(5)
S += tsd3
start3 = len(S) + 1
S += IR["T3"] + random_dna(10)
S += tsd4
start4 = len(S) + 1
S += IR["T4"] + random_dna(20) + reverse_complement(IR["T4"])
real["T4"] = (start4, len(S))
S += tsd4 + random_dna(10) + reverse_complement(IR["T3"])
real["T3"] = (start3, len(S))
S += tsd3

S += random_dna(30)

print("Sequence length:", len(S))
print(S, "\n")

# Detector: it does NOT know the IRs or the positions.
# For every piece of IR_LEN letters, look for its reverse complement further on
# (at most MAX_LEN letters away). If found -> a transposon from that piece to the mirrored piece.
MAX_LEN = 150

detected = []
i = 0
while i <= len(S) - IR_LEN:
    left = S[i:i + IR_LEN]
    j = S.find(reverse_complement(left), i + IR_LEN, i + MAX_LEN)
    if j != -1:
        detected.append((i + 1, j + IR_LEN))  # positions start from 1
        i += IR_LEN  # skip the rest of this inverted repeat
    else:
        i += 1

print("Real transposons:")
for name, (start, end) in real.items():
    print(" ", name, "start =", start, "end =", end)

print("\nDetected transposons:")
for start, end in detected:
    print("  start =", start, "end =", end)

# Comparing (a detection counts if its start and end are at most 2 letters away from the real ones)
print("\nComparison:")
for name, (start, end) in real.items():
    result = "NOT found"
    for d_start, d_end in detected:
        if abs(d_start - start) <= 2 and abs(d_end - end) <= 2:
            if (d_start, d_end) == (start, end):
                result = "found exactly"
            else:
                result = "found, but the ends are shifted by " + str(abs(d_start - start)) + " letter(s)"
    print(" ", name, "->", result)
