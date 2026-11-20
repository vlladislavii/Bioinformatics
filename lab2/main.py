# Lab 2 - Melting temperature (Tm) of a DNA sequence
# Tm = 81.5 + 16.6 * log10([Na+]) + 41 * (%GC) - 600 / length

import math

S = "ACGCGTGCCA"
Na = 0.05  # Na+ concentration (mol/L)

gc = 0
for letter in S:
    match letter:
        case "G" | "C":
            gc += 1

length = len(S)
percent_gc = gc / length

Tm = 81.5 + 16.6 * math.log10(Na) + 41 * percent_gc - 600 / length

print("Sequence:", S)
print("Tm =", round(Tm, 2))
