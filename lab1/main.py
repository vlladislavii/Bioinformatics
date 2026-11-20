# Lab 1 - Frequency of the letters A, C, G, T in a DNA sequence S

S = "TGCATTACTTGCAGT"

a = c = g = t = 0

for letter in S:
    match letter:
        case "A":
            a += 1
        case "C":
            c += 1
        case "G":
            g += 1
        case "T":
            t += 1

n = len(S)

print("Sequence:", S)
print("Length:", n)
print("A:", a, "->", round(a / n * 100, 2), "%")
print("C:", c, "->", round(c / n * 100, 2), "%")
print("G:", g, "->", round(g / n * 100, 2), "%")
print("T:", t, "->", round(t / n * 100, 2), "%")
