# Lab 5 - Gel electrophoresis simulation (NCBI: M77789.2, pUC19, 2686 bp)

import math
import os
import random
import matplotlib.pyplot as plt

# Reading the DNA sequence downloaded from NCBI
S = ""
with open("sequence.fasta") as file:
    for line in file:
        if not line.startswith(">"):
            S += line.strip()

print("Sequence length:", len(S))

# Taking 10 random samples
samples = []
for i in range(10):
    length = random.randint(100, min(3000, len(S)))
    start = random.randint(0, len(S) - length)
    samples.append(S[start:start + length])

# 4. Migration on the gel
# Molecular weight: MW = length * 650 Da
# Distance: short fragments travel far, long fragments stay near the well
#   d = (log10(3000) - log10(L)) / (log10(3000) - log10(100))   -> 0 = top, 1 = bottom
distances = []
for i in range(10):
    L = len(samples[i])
    mw = L * 650
    d = (math.log10(3000) - math.log10(L)) / (math.log10(3000) - math.log10(100))
    distances.append(d)
    print("Sample", i + 1, ":", L, "bp, MW =", mw, "Da, distance =", round(d, 2))

plt.figure(figsize=(10, 6))
ax = plt.gca()
ax.set_facecolor("black")

for i in range(10):
    lane = i + 1
    plt.plot([lane - 0.3, lane + 0.3], [-0.05, -0.05], color="gray", linewidth=8)             # well
    plt.plot([lane - 0.3, lane + 0.3], [distances[i], distances[i]], color="orange", linewidth=6)  # band

plt.xlim(0.3, 10.7)
plt.ylim(1.05, -0.1)  # inverted: wells at the top, + electrode at the bottom
plt.xticks(range(1, 11), [str(len(s)) + " bp" for s in samples])
plt.yticks([])
plt.xlabel("Samples (fragment length)")
plt.title("Gel electrophoresis simulation (wells on top, DNA moves down)")
plt.tight_layout()

os.makedirs("out", exist_ok=True)
plt.savefig("out/gel.png")
plt.show()
