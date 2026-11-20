# Lab 7 - Comparing two genomes side by side on the gel
#   Genome 1: SARS coronavirus (NC_004718.3, 29751 bp)
#   Genome 2: SARS coronavirus 2 (NC_045512.2, 29903 bp)

import math
import os
import matplotlib.pyplot as plt

ENZYMES = {
    "EcoRI": ("GAATTC", 1),
    "BamHI": ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
}


def read_fasta(filename):
    sequence = ""
    with open(filename) as file:
        for line in file:
            if not line.startswith(">"):
                sequence += line.strip()
    return sequence


def cut(sequence, site, cleavage):
    cuts = []
    for i in range(len(sequence) - len(site) + 1):
        if sequence[i:i + len(site)] == site:
            cuts.append(i + cleavage)
    fragments = []
    start = 0
    for c in cuts:
        fragments.append(sequence[start:c])
        start = c
    fragments.append(sequence[start:])
    return cuts, fragments


genomes = {
    "SARS-CoV": read_fasta("sars_cov.fasta"),
    "SARS-CoV-2": read_fasta("sars_cov_2.fasta"),
}

# Cutting both genomes with each enzyme; the same enzyme is put side by side for both genomes
LADDER = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 500, 250, 100]
lanes = {"Ladder": LADDER}

for enzyme, (site, cleavage) in ENZYMES.items():
    for genome_name, sequence in genomes.items():
        cuts, fragments = cut(sequence, site, cleavage)
        lanes[genome_name + "\n" + enzyme] = [len(f) for f in fragments]
        print(genome_name, enzyme, "->", len(fragments), "fragments:", [len(f) for f in fragments])
    print()

# Gel electrophoresis
plt.figure(figsize=(12, 8))
ax = plt.gca()
ax.set_facecolor("black")

for x, (name, sizes) in enumerate(lanes.items()):
    color = "orange" if "SARS-CoV-2" in name else "deepskyblue"
    if name == "Ladder":
        color = "white"
    for size in sizes:
        y = math.log10(size)
        plt.plot([x - 0.3, x + 0.3], [y, y], color=color, linewidth=4)
        if name == "Ladder":
            plt.text(x - 0.4, y, str(size), color="white", ha="right", va="center", fontsize=8)

plt.xlim(-1, len(lanes))
plt.ylim(1.8, 4.6)
plt.xticks(range(len(lanes)), list(lanes.keys()), fontsize=8)
plt.yticks([])
plt.title("SARS-CoV (blue) vs SARS-CoV-2 (orange) cut with EcoRI, BamHI, HindIII")
plt.tight_layout()

os.makedirs("out", exist_ok=True)
plt.savefig("out/gel_task2.png")
plt.show()
