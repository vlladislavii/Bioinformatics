# Lab 8 - DNA patterns (Objective Digital Stains) from gene promoters
#
# For each sliding window: x = C+G %, y = information content (IC %).
# The functions follow the teacher's code (github.com/Gagniuc - PromKappa, Self-sequence-alignment).

import os
import matplotlib.pyplot as plt

S = "CGGACTGATCTATCT" + "A" * 27 + "CGTAGCATCTATCGATCTATCTAGCGATCTA"
WINDOW = 30


def cg_content(s):
    return (s.count("C") + s.count("G")) / len(s) * 100


def information_content(s):
    # Self-sequence alignment: compare s with itself shifted by u letters,
    # take the % of matching letters, then the average over all shifts.
    total = 0
    for u in range(1, len(s)):
        matches = 0
        for i in range(len(s) - u):
            if s[i] == s[u + i]:  # f(x_i, x_u+i) = 1 if equal, else 0
                matches += 1
        total += matches / (len(s) - u) * 100
    return total / (len(s) - 1)


def pattern(sequence):
    points = []
    for i in range(len(sequence) - WINDOW + 1):
        window = sequence[i:i + WINDOW]
        points.append((cg_content(window), information_content(window)))
    return points


def center_of_weight(points):
    x = sum(p[0] for p in points) / len(points)
    y = sum(p[1] for p in points) / len(points)
    return x, y


def plot_pattern(points, center, title, filename, ylabel="Information content (%)"):
    plt.figure(figsize=(6, 6))
    plt.scatter([p[0] for p in points], [p[1] for p in points], s=15, label="sliding windows")
    plt.scatter(center[0], center[1], color="red", s=80, label="center of weight")
    plt.xlim(0, 100)
    plt.ylim(0, 100)
    plt.xlabel("C+G (%)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig("out/" + filename)


def read_fasta(filename):
    sequences = {}
    name = ""
    with open(filename) as file:
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                name = line.split(";")[0].split()[-1]  # e.g. ">EP 0 (+) INS_1; range..." -> INS_1
                sequences[name] = ""
            else:
                sequences[name] += line
    return sequences


os.makedirs("out", exist_ok=True)

# test sequence S
print("Test sequence S:", S, "(" + str(len(S)) + " letters)")
print("C+G of S =", round(cg_content(S), 2), "(board expects ~29.27)")
print("IC  of S =", round(information_content(S), 2), "(board expects ~27.53)")
print("First window:", S[:WINDOW], "-> C+G =", round(cg_content(S[:WINDOW]), 2),
      ", IC =", round(information_content(S[:WINDOW]), 2))

# pattern of S and its center of weight
points_S = pattern(S)
center_S = center_of_weight(points_S)
print("\nS: windows =", len(points_S), ", center of weight = (C+G", round(center_S[0], 2),
      ", IC", round(center_S[1], 2), ")")
plot_pattern(points_S, center_S, "Pattern of the test sequence S", "pattern_S.png")

# the same pattern with sigma(S) = 100 - IC (the formula written on the board)
points_sigma = [(x, 100 - y) for x, y in points_S]
center_sigma = center_of_weight(points_sigma)
print("S (sigma): center of weight = (C+G", round(center_sigma[0], 2), ", sigma", round(center_sigma[1], 2), ")")
plot_pattern(points_sigma, center_sigma, "Pattern of S with sigma(S) = 100 - IC", "pattern_S_sigma.png",
             ylabel="sigma(S) = 100 - IC (%)")

# patterns of real gene promoters (from the teacher's PromKappa promoter file)
centers = {"S": center_S}
promoters = read_fasta("promoters.fasta")
for name, sequence in promoters.items():
    points = pattern(sequence)
    centers[name] = center_of_weight(points)
    print(name + ": windows =", len(points), ", center of weight = (C+G", round(centers[name][0], 2),
          ", IC", round(centers[name][1], 2), ")")
    plot_pattern(points, centers[name], "Pattern of the " + name + " promoter", "pattern_" + name + ".png")

# the center of each pattern on a second chart
plt.figure(figsize=(6, 6))
for name, (x, y) in centers.items():
    plt.scatter(x, y, s=60)
    plt.text(x + 1, y + 1, name, fontsize=9)
plt.xlim(0, 100)
plt.ylim(0, 100)
plt.xlabel("C+G (%)")
plt.ylabel("Information content (%)")
plt.title("Centers of weight of all patterns")
plt.tight_layout()
plt.savefig("out/centers.png")

plt.show()
