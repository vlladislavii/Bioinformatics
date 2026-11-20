# Lab 3 - Codon frequencies of the SARS-CoV-2 genome (NCBI: NC_045512.2)

import os
import matplotlib.pyplot as plt
import ollama

# Genetic Code table
GENETIC_CODE = {
    "TTT": "Phenylalanine", "TTC": "Phenylalanine", "TTA": "Leucine", "TTG": "Leucine",
    "CTT": "Leucine", "CTC": "Leucine", "CTA": "Leucine", "CTG": "Leucine",
    "ATT": "Isoleucine", "ATC": "Isoleucine", "ATA": "Isoleucine", "ATG": "Methionine",
    "GTT": "Valine", "GTC": "Valine", "GTA": "Valine", "GTG": "Valine",
    "TCT": "Serine", "TCC": "Serine", "TCA": "Serine", "TCG": "Serine",
    "CCT": "Proline", "CCC": "Proline", "CCA": "Proline", "CCG": "Proline",
    "ACT": "Threonine", "ACC": "Threonine", "ACA": "Threonine", "ACG": "Threonine",
    "GCT": "Alanine", "GCC": "Alanine", "GCA": "Alanine", "GCG": "Alanine",
    "TAT": "Tyrosine", "TAC": "Tyrosine", "TAA": "Stop", "TAG": "Stop",
    "CAT": "Histidine", "CAC": "Histidine", "CAA": "Glutamine", "CAG": "Glutamine",
    "AAT": "Asparagine", "AAC": "Asparagine", "AAA": "Lysine", "AAG": "Lysine",
    "GAT": "Aspartic acid", "GAC": "Aspartic acid", "GAA": "Glutamic acid", "GAG": "Glutamic acid",
    "TGT": "Cysteine", "TGC": "Cysteine", "TGA": "Stop", "TGG": "Tryptophan",
    "CGT": "Arginine", "CGC": "Arginine", "CGA": "Arginine", "CGG": "Arginine",
    "AGT": "Serine", "AGC": "Serine", "AGA": "Arginine", "AGG": "Arginine",
    "GGT": "Glycine", "GGC": "Glycine", "GGA": "Glycine", "GGG": "Glycine",
}

def codon_to_amino_acid(codon):
    return GENETIC_CODE[codon]


# Reading the genome downloaded from NCBI
S = ""
with open("sequence.fasta") as file:
    for line in file:
        if not line.startswith(">"):  # skip the header line
            S += line.strip()

print("Genome length:", len(S))

# Codon frequencies in all 3 frames
counts = {}
for frame in range(3):
    for i in range(frame, len(S) - 2, 3):
        codon = S[i:i + 3]
        counts[codon] = counts.get(codon, 0) + 1

sorted_codons = sorted(counts.items(), key=lambda x: x[1], reverse=True)
top10 = sorted_codons[:10]

print("\nTop 10 codons:")
for codon, count in top10:
    print(codon, count, codon_to_amino_acid(codon))

labels = [codon + "\n" + codon_to_amino_acid(codon) for codon, count in top10]
values = [count for codon, count in top10]

plt.figure(figsize=(12, 5))
plt.bar(labels, values)
plt.title("Top 10 most frequent codons (all 3 frames)")
plt.xlabel("Codon")
plt.ylabel("Frequency")
plt.tight_layout()
os.makedirs("out", exist_ok=True)
plt.savefig("out/top10_codons.png")
plt.show()

# Asking AI
prompt = "These are the 10 most frequent codons in the SARS-CoV-2 genome and their amino acids:\n"
for codon, count in top10:
    prompt += codon + " - " + codon_to_amino_acid(codon) + "\n"
prompt += "Make a list with the foods that contain the lowest amount of these amino acids."

print("\nPrompt for the AI:\n")
print(prompt)

# AI can be wrong - its answer should be checked before being trusted (since the model is small)
answer = ollama.generate(model="llama3.2:3b", prompt=prompt)
print("\nAI answer:\n")
print(answer["response"])
