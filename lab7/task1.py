# Lab 7 - Restriction enzymes + gel electrophoresis (NCBI: M77789.2, pUC19, 2686 bp)

import math
import os
import matplotlib.pyplot as plt

# Restriction enzymes: name -> (recognition sequence, cleavage site = letters before the cut)
#   EcoRI   G^AATTC
#   BamHI   G^GATCC
#   HindIII A^AGCTT
ENZYMES = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
}

# Reading the DNA sequence downloaded from NCBI
S = ""
with open("sequence.fasta") as file:
    for line in file:
        if not line.startswith(">"):
            S += line.strip()

print("Sequence length:", len(S))


def cut(sequence, site, cleavage):
    # positions where the enzyme cuts
    cuts = []
    for i in range(len(sequence) - len(site) + 1):
        if sequence[i:i + len(site)] == site:
            cuts.append(i + cleavage)
    # fragments between the cuts
    fragments = []
    start = 0
    for c in cuts:
        fragments.append(sequence[start:c])
        start = c
    fragments.append(sequence[start:])
    return cuts, fragments


# Cutting the sequence with each enzyme
lanes = {}
for name, (site, cleavage) in ENZYMES.items():
    cuts, fragments = cut(S, site, cleavage)
    lanes[name] = [len(f) for f in fragments]
    print("\n" + name, "(" + site + ") cuts after positions:", cuts)
    print("  fragment lengths:", lanes[name])

# Gel electrophoresis: long fragments stay at the top, short ones travel down (log10 of length)
LADDER = [3000, 2000, 1500, 1000, 500, 250, 100]

plt.figure(figsize=(6, 7))
ax = plt.gca()
ax.set_facecolor("black")

all_lanes = {"Ladder": LADDER, **lanes}
for x, (name, sizes) in enumerate(all_lanes.items()):
    for size in sizes:
        y = math.log10(size)
        plt.plot([x - 0.3, x + 0.3], [y, y], color="orange", linewidth=5)
        if name == "Ladder":
            plt.text(x - 0.4, y, str(size), color="white", ha="right", va="center", fontsize=8)
        else:
            plt.text(x + 0.35, y, str(size), color="cyan", ha="left", va="center", fontsize=8)

plt.xlim(-1, len(all_lanes))
plt.ylim(1.8, 3.6)
plt.xticks(range(len(all_lanes)), list(all_lanes.keys()))
plt.yticks([])
plt.title("pUC19 cut with EcoRI, BamHI, HindIII")
plt.tight_layout()

os.makedirs("out", exist_ok=True)
plt.savefig("out/gel_task1.png")
plt.show()
