# Lab 9 - finding a DNA motif with PFM -> PWM -> LLM and a sliding window

import os
import math
import matplotlib.pyplot as plt

MOTIFS = [
    "ATTCCGAAGTG",
    "TACGCAATCGA",
    "CGTGCTTACAG",
    "GTAATCTGAGC",
    "TGATCAGGCTA",
    "GACCTAAGGAT",
    "TCGATATGCCG",
]
S = "AGTCAATGGCTAAGCCATTA"
FASTA_FILE = "promoters.fasta"

LETTERS = "ACGT"
NULL_MODEL = 0.25
ZERO_SCORE = -20
L = len(MOTIFS[0])
N = len(MOTIFS)


def build_pfm():
    pfm = {letter: [0] * L for letter in LETTERS}
    for motif in MOTIFS:
        for i in range(L):
            pfm[motif[i]][i] += 1
    return pfm


def build_pwm(pfm):
    return {letter: [count / N for count in pfm[letter]] for letter in LETTERS}


def build_llm(pwm):
    llm = {}
    for letter in LETTERS:
        llm[letter] = []
        for p in pwm[letter]:
            if p == 0:
                llm[letter].append(ZERO_SCORE)
            else:
                llm[letter].append(math.log(p / NULL_MODEL))
    return llm


def print_matrix(name, matrix, decimals=2):
    print(name)
    print("    " + "".join(f"{i + 1:>7}" for i in range(L)))
    for letter in LETTERS:
        print(f"{letter}   " + "".join(f"{value:>7.{decimals}f}" for value in matrix[letter]))
    print()


def score_windows(sequence, llm):
    scores = []
    for start in range(len(sequence) - L + 1):
        window = sequence[start:start + L]
        score = 0
        for i in range(L):
            score += llm[window[i]][i]
        scores.append(score)
    return scores


def read_first_fasta(filename):
    sequence = ""
    with open(filename) as file:
        file.readline()  # skip the header of the first record
        for line in file:
            if line.startswith(">"):
                break
            sequence += line.strip().upper()
    return sequence


def plot_scores(scores, title, filename):
    plt.figure(figsize=(10, 4))
    plt.bar(range(1, len(scores) + 1), scores)
    plt.xlabel("Window start position")
    plt.ylabel("Score (sum of LLM)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig("out/" + filename)


def report(name, sequence, llm, filename, show_all):
    scores = score_windows(sequence, llm)
    print(name, "-", len(sequence), "letters,", len(scores), "windows")
    if show_all:
        for i in range(len(scores)):
            print(f"  S{i + 1} = {sequence[i:i + L]} = {scores[i]:.2f}")
    best = scores.index(max(scores))
    print(f"  Best window: S{best + 1} at position {best + 1} = "
          f"{sequence[best:best + L]}, score {scores[best]:.2f}")
    print()
    plot_scores(scores, "Window scores - " + name, filename)


os.makedirs("out", exist_ok=True)

pfm = build_pfm()
pwm = build_pwm(pfm)
llm = build_llm(pwm)
print_matrix("PFM (counts)", pfm, 0)
print_matrix("PWM (probabilities)", pwm)
print_matrix("LLM = ln(P / 0.25), 0 -> -20", llm)

report("Default S", S, llm, "scores_default.png", True)
report(FASTA_FILE + " (first record)", read_first_fasta(FASTA_FILE), llm, "scores_fasta.png", False)

plt.show()
