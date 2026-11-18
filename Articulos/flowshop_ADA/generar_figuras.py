"""Genera las figuras del informe flowshop.tex y verifica los ejemplos numéricos.

Uso:
    python3 generar_figuras.py

Crea los PNG en la carpeta fig/ relativa a este script.
"""

from itertools import permutations
from math import factorial
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


BASE = Path(__file__).resolve().parent
FIG_DIR = BASE / "fig"
FIG_DIR.mkdir(exist_ok=True)


def makespan(seq, p):
    """Calcula Cmax para una secuencia mediante la recurrencia de flow shop."""
    m = len(p[0])
    c = [0] * m
    for j in seq:
        c[0] += p[j][0]
        for i in range(1, m):
            c[i] = max(c[i], c[i - 1]) + p[j][i]
    return c[-1]


def completion_table(seq, p):
    n = len(seq)
    m = len(p[0])
    table = [[0] * m for _ in range(n)]
    for k, j in enumerate(seq):
        for i in range(m):
            previous_job = table[k - 1][i] if k else 0
            previous_machine = table[k][i - 1] if i else 0
            table[k][i] = max(previous_job, previous_machine) + p[j][i]
    return table


def johnson(a, b):
    """Regla de Johnson para dos máquinas."""
    first = sorted((j for j in range(len(a)) if a[j] <= b[j]), key=lambda j: a[j])
    second = sorted((j for j in range(len(a)) if a[j] > b[j]), key=lambda j: -b[j])
    return first + second


def palmer(p):
    m = len(p[0])
    weights = [2 * i - (m + 1) for i in range(1, m + 1)]
    scores = [sum(weights[i] * p[j][i] for i in range(m)) for j in range(len(p))]
    return sorted(range(len(p)), key=lambda j: -scores[j]), scores


def cds(p):
    """Devuelve la mejor secuencia de CDS y su makespan."""
    m = len(p[0])
    best = None
    for k in range(1, m):
        a = [sum(row[:k]) for row in p]
        b = [sum(row[k:]) for row in p]
        seq = johnson(a, b)
        cmax = makespan(seq, p)
        if best is None or cmax < best[0]:
            best = (cmax, seq)
    return best


def neh(p):
    order = sorted(range(len(p)), key=lambda j: -sum(p[j]))
    seq = [order[0]]
    for j in order[1:]:
        candidates = [(makespan(seq[:pos] + [j] + seq[pos:], p), pos)
                      for pos in range(len(seq) + 1)]
        best_cmax = min(c for c, _ in candidates)
        pos = next(pos for c, pos in candidates if c == best_cmax)
        seq = seq[:pos] + [j] + seq[pos:]
    return seq


def lower_bound(p):
    n = len(p)
    m = len(p[0])
    bounds = []
    for k in range(m):
        before = min(sum(row[i] for i in range(k)) for row in p)
        after = min(sum(row[i] for i in range(k + 1, m)) for row in p)
        bounds.append(before + sum(row[k] for row in p) + after)
    return max(bounds)


def gantt(seq, p, title, filename):
    table = completion_table(seq, p)
    m = len(p[0])
    colors = plt.cm.tab20.colors
    fig, ax = plt.subplots(figsize=(9.5, 3.8))
    for k, j in enumerate(seq):
        start_previous = table[k][0] - p[j][0]
        for i in range(m):
            start = table[k][i] - p[j][i]
            ax.barh(i, p[j][i], left=start, height=0.65,
                    color=colors[j % len(colors)], edgecolor="black")
            ax.text(start + p[j][i] / 2, i, f"J{j + 1}",
                    ha="center", va="center", fontsize=8)
    ax.set_yticks(range(m), [f"M{i + 1}" for i in range(m)])
    ax.set_xlabel("Tiempo")
    ax.set_title(title)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(FIG_DIR / filename, dpi=180, bbox_inches="tight")
    plt.close(fig)


# Ejemplo de dos máquinas
johnson_a = [3, 5, 1, 6, 7]
johnson_b = [6, 2, 2, 6, 5]
johnson_seq = johnson(johnson_a, johnson_b)
johnson_p = [[a, b] for a, b in zip(johnson_a, johnson_b)]
print("Johnson:", [chr(65 + j) for j in johnson_seq], "Cmax =",
      makespan(johnson_seq, johnson_p))

# Diagrama de Gantt de Johnson
table = completion_table(johnson_seq, johnson_p)
fig, ax = plt.subplots(figsize=(9.5, 3.2))
colors = plt.cm.tab20.colors
for k, j in enumerate(johnson_seq):
    for i in range(2):
        start = table[k][i] - johnson_p[j][i]
        ax.barh(i, johnson_p[j][i], left=start, height=0.65,
                color=colors[j % len(colors)], edgecolor="black")
        ax.text(start + johnson_p[j][i] / 2, i, chr(65 + j),
                ha="center", va="center", fontsize=9)
ax.set_yticks([0, 1], ["M1", "M2"])
ax.set_xlabel("Tiempo")
ax.set_title("Secuencia óptima de Johnson: C--A--D--E--B")
ax.axvline(24, color="red", linestyle="--", label="Cmax = 24")
ax.legend(loc="lower right")
ax.invert_yaxis()
ax.grid(axis="x", linestyle="--", alpha=0.35)
fig.tight_layout()
fig.savefig(FIG_DIR / "gantt_johnson.png", dpi=180, bbox_inches="tight")
plt.close(fig)

# Ejemplo de cuatro máquinas
P = [[5, 9, 8, 10], [9, 3, 10, 1], [9, 4, 5, 8],
     [4, 8, 8, 7], [1, 5, 6, 8]]
palmer_seq, palmer_scores = palmer(P)
cds_cmax, cds_seq = cds(P)
neh_seq = neh(P)
original = list(range(len(P)))
print("Palmer:", [f"J{j + 1}" for j in palmer_seq], palmer_scores,
      "Cmax =", makespan(palmer_seq, P))
print("CDS:", [f"J{j + 1}" for j in cds_seq], "Cmax =", cds_cmax)
print("NEH:", [f"J{j + 1}" for j in neh_seq], "Cmax =", makespan(neh_seq, P))
print("Original: Cmax =", makespan(original, P))
print("Cota inferior =", lower_bound(P))

best_cmax = None
best_seq = None
for seq in permutations(range(len(P))):
    cmax = makespan(seq, P)
    if best_cmax is None or cmax < best_cmax:
        best_cmax = cmax
        best_seq = seq
print("Enumeración:", [f"J{j + 1}" for j in best_seq],
      "Cmax =", best_cmax, "permits =", factorial(len(P)))

# Figuras de Gantt y espacio de búsqueda
gantt(original, P, "Orden original: J1--J2--J3--J4--J5",
      "gantt_comparacion.png")
# Comparación: sustituye temporalmente con dos paneles para reflejar ambos órdenes
# Se rehace la figura como dos paneles, igual que en el documento.
fig, axes = plt.subplots(2, 1, figsize=(9.5, 7), sharex=True)
for ax, seq, title in zip(axes,
                          [original, neh_seq],
                          ["Orden original", "Secuencia NEH: J5--J4--J3--J1--J2"]):
    table = completion_table(seq, P)
    for k, j in enumerate(seq):
        for i in range(4):
            start = table[k][i] - P[j][i]
            ax.barh(i, P[j][i], left=start, height=0.65,
                    color=colors[j % len(colors)], edgecolor="black")
            ax.text(start + P[j][i] / 2, i, f"J{j + 1}",
                    ha="center", va="center", fontsize=8)
    ax.set_yticks(range(4), [f"M{i + 1}" for i in range(1, 5)])
    ax.set_title(title)
    ax.invert_yaxis()
    ax.grid(axis="x", linestyle="--", alpha=0.35)
axes[-1].set_xlabel("Tiempo")
fig.tight_layout()
fig.savefig(FIG_DIR / "gantt_comparacion.png", dpi=180, bbox_inches="tight")
plt.close(fig)

n_values = list(range(2, 21))
log_permutations = [factorial(n).bit_length() for n in n_values]
# Usa log10 directamente sin depender de la precisión factorial.
import math
log_permutations = [math.log10(math.factorial(n)) for n in n_values]
fig, ax = plt.subplots(figsize=(8.5, 5))
ax.plot(n_values, log_permutations, marker="o", linewidth=1.8,
        label="Permutacional: $n!$")
ax.plot(n_values, [n * math.log10(n) for n in n_values], linestyle="--",
        label="Caso general: $(n!)^m$ (ejemplo $m=n$)")
ax.set_xlabel("Número de trabajos $n$")
ax.set_ylabel("Logaritmo base 10 del número de secuencias")
ax.set_title("Crecimiento del espacio de búsqueda")
ax.grid(True, linestyle="--", alpha=0.35)
ax.legend()
fig.tight_layout()
fig.savefig(FIG_DIR / "espacio_busqueda.png", dpi=180, bbox_inches="tight")
plt.close(fig)

print("Figuras creadas en:", FIG_DIR)
