# Lab 4 - Rebuilding a genome from random samples (NCBI: NC_001422.1, phage phiX174, 5386 bp)

import random

# Reading the genome downloaded from NCBI
S = ""
with open("sequence.fasta") as file:
    for line in file:
        if not line.startswith(">"):
            S += line.strip()

print("Genome length:", len(S))

# Taking 2000 random samples of 100-120 bases and storing them in a list
samples = []
for i in range(2000):
    length = random.randint(100, 120)
    start = random.randint(0, len(S) - length)
    samples.append(S[start:start + length])

print("Number of samples:", len(samples))

# Rebuilding the sequence
# Start from one sample and keep gluing the sample that overlaps the most with the end.
# The overlap must be at least MIN_OVERLAP letters.
MIN_OVERLAP = 20

def extend_right(contig, reads):
    while True:
        tail = contig[-120:]
        best_read = None
        best_overlap = 0
        for r in reads:
            pos = tail.find(r[:MIN_OVERLAP])
            if pos == -1:
                continue
            overlap = len(tail) - pos
            # the read must match the end of the contig and add new letters
            if overlap < len(r) and tail[pos:] == r[:overlap] and overlap > best_overlap:
                best_read = r
                best_overlap = overlap
        if best_read is None:
            return contig
        contig += best_read[best_overlap:]


rebuilt = extend_right(samples[0], samples)

# To extend to the left: reverse everything, extend to the right, reverse back
reversed_samples = [r[::-1] for r in samples]
rebuilt = extend_right(rebuilt[::-1], reversed_samples)[::-1]

# Results
print("Rebuilt length:", len(rebuilt))
print("Rebuilt is part of the original genome:", rebuilt in S)
print("Rebuilt is identical to the original genome:", rebuilt == S)
print("Missing letters:", len(S) - len(rebuilt))

# Main problem of this approach:
# The samples have no position - we only know their letters. We rebuild by guessing
# the order from overlaps. If a part of the genome is repeated (a repeat longer than
# the overlap), several samples fit equally well and the algorithm cannot know which
# one is correct, so it can glue the wrong pieces, skip or collapse parts of the genome.
# Also: parts of the genome that no sample covers (especially the two ends) cannot be
# rebuilt, and comparing every sample with every other one is slow for big genomes.
