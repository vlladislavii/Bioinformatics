# Lab 8 - Aligning two DNA sequences with the Needleman-Wunsch algorithm

S1 = "ACCGTGAAGCCAATAC"
S2 = "AGCGTGCAGCCAATAC"

MATCH = 1
MISMATCH = -1
GAP = -1

n = len(S1)
m = len(S2)

# Score matrix: (n+1) x (m+1), first row and column are gaps
score = [[0] * (m + 1) for _ in range(n + 1)]
for i in range(1, n + 1):
    score[i][0] = i * GAP
for j in range(1, m + 1):
    score[0][j] = j * GAP

# Filling the matrix: each cell takes the best of diagonal / up / left
for i in range(1, n + 1):
    for j in range(1, m + 1):
        if S1[i - 1] == S2[j - 1]:
            diagonal = score[i - 1][j - 1] + MATCH
        else:
            diagonal = score[i - 1][j - 1] + MISMATCH
        up = score[i - 1][j] + GAP
        left = score[i][j - 1] + GAP
        score[i][j] = max(diagonal, up, left)

# Traceback: from the bottom-right corner back to the top-left corner
aligned1 = ""
aligned2 = ""
i = n
j = m
while i > 0 or j > 0:
    if i > 0 and j > 0:
        if S1[i - 1] == S2[j - 1]:
            diagonal = score[i - 1][j - 1] + MATCH
        else:
            diagonal = score[i - 1][j - 1] + MISMATCH
    else:
        diagonal = None
    if diagonal is not None and score[i][j] == diagonal:
        aligned1 = S1[i - 1] + aligned1
        aligned2 = S2[j - 1] + aligned2
        i -= 1
        j -= 1
    elif i > 0 and score[i][j] == score[i - 1][j] + GAP:
        aligned1 = S1[i - 1] + aligned1
        aligned2 = "-" + aligned2
        i -= 1
    else:
        aligned1 = "-" + aligned1
        aligned2 = S2[j - 1] + aligned2
        j -= 1

marks = ""
matches = 0
for a, b in zip(aligned1, aligned2):
    if a == b:
        marks += "|"
        matches += 1
    elif a == "-" or b == "-":
        marks += " "
    else:
        marks += "."

print("Score matrix:")
print("      " + "   ".join("-" + S2))
for i in range(n + 1):
    letter = "-" if i == 0 else S1[i - 1]
    print(letter, " ".join(f"{v:3d}" for v in score[i]))

print("\nAlignment:")
print(aligned1)
print(marks)
print(aligned2)
print("\nScore:", score[n][m])
print("Matches:", matches, "of", len(aligned1), "->", round(matches / len(aligned1) * 100, 2), "% similarity")
