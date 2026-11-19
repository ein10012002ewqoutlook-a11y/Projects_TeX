#!/usr/bin/env python3
"""
generar_figuras.py
Genera las figuras del informe flowshop.tex y verifica los calculos del
ejemplo numerico (Johnson, Palmer, CDS, NEH y busqueda exhaustiva).

Uso:
    pip install matplotlib
    python3 generar_figuras.py

Crea la carpeta fig/ junto al script con:
    fig/gantt_johnson.png
    fig/gantt_comparacion.png
    fig/espacio_busqueda.png
"""
import itertools
import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fig")
os.makedirs(OUT, exist_ok=True)

COLORES = ["#8ecae6", "#ffb703", "#b7e4c7", "#f4a6a6", "#cdb4db",
           "#f1e3a1", "#a8dadc", "#e5b8a0"]

# ---------------------------------------------------------------------------
# Datos del ejemplo (4 maquinas, 5 trabajos): P[j][i] = tiempo del trabajo j
# en la maquina i
# ---------------------------------------------------------------------------
P = [
    [5, 9, 8, 10],   # J1
    [9, 3, 10, 1],   # J2
    [9, 4, 5, 8],    # J3
    [4, 8, 8, 7],    # J4
    [1, 5, 6, 8],    # J5
]
# Ejemplo de 2 maquinas para la regla de Johnson
J2 = [(3, 6), (5, 2), (1, 2), (6, 6), (7, 5)]   # A, B, C, D, E
NOM_J2 = ["A", "B", "C", "D", "E"]


def matriz_c(seq, p):
    """Tiempos de finalizacion C[k][i] del k-esimo trabajo de la secuencia."""
    m = len(p[0])
    c = [[0] * m for _ in seq]
    for k, j in enumerate(seq):
        for i in range(m):
            previo_trabajo = c[k - 1][i] if k > 0 else 0
            previo_maquina = c[k][i - 1] if i > 0 else 0
            c[k][i] = max(previo_trabajo, previo_maquina) + p[j][i]
    return c


def cmax(seq, p):
    return matriz_c(seq, p)[-1][-1] if seq else 0


def johnson(a, b):
    s1 = sorted([j for j in range(len(a)) if a[j] <= b[j]], key=lambda j: a[j])
    s2 = sorted([j for j in range(len(a)) if a[j] > b[j]], key=lambda j: -b[j])
    return s1 + s2


def palmer(p):
    m = len(p[0])
    pend = [sum((2 * (i + 1) - m - 1) * p[j][i] for i in range(m)) for j in range(len(p))]
    return sorted(range(len(p)), key=lambda j: -pend[j]), pend


def cds(p):
    m = len(p[0])
    mejor = None
    detalle = []
    for k in range(1, m):
        a = [sum(p[j][:k]) for j in range(len(p))]
        b = [sum(p[j][m - k:]) for j in range(len(p))]
        s = johnson(a, b)
        c = cmax(s, p)
        detalle.append((k, s, c))
        if mejor is None or c < mejor[1]:
            mejor = (s, c)
    return mejor[0], detalle


def neh(p, log=None):
    n = len(p)
    orden = sorted(range(n), key=lambda j: -sum(p[j]))
    seq = [orden[0]]
    if log is not None:
        log.append(("orden inicial", orden, None))
    for j in orden[1:]:
        evals = []
        mejor = None
        for pos in range(len(seq) + 1):
            cand = seq[:pos] + [j] + seq[pos:]
            c = cmax(cand, p)
            evals.append((pos + 1, cand, c))
            if mejor is None or c < mejor[0]:   # empate: primera posicion
                mejor = (c, cand)
        seq = mejor[1]
        if log is not None:
            log.append((j, evals, mejor[0]))
    return seq


def optimo(p):
    mejor = None
    for perm in itertools.permutations(range(len(p))):
        c = cmax(list(perm), p)
        if mejor is None or c < mejor[0]:
            mejor = (c, list(perm))
    return mejor


def cota_inferior(p):
    n, m = len(p), len(p[0])
    lb = 0
    for k in range(m):
        cab = min(sum(p[j][:k]) for j in range(n))
        cola = min(sum(p[j][k + 1:]) for j in range(n))
        lb = max(lb, cab + sum(p[j][k] for j in range(n)) + cola)
    return lb


def gantt(ax, seq, p, nombres, titulo):
    c = matriz_c(seq, p)
    m = len(p[0])
    for k, j in enumerate(seq):
        for i in range(m):
            fin = c[k][i]
            ini = fin - p[j][i]
            ax.barh(i, p[j][i], left=ini, height=0.62, color=COLORES[j % len(COLORES)],
                    edgecolor="black", linewidth=0.6)
            etiqueta = nombres[j] if p[j][i] >= 2 else nombres[j][-1]
            ax.text(ini + p[j][i] / 2, i, etiqueta, ha="center", va="center",
                    fontsize=8 if p[j][i] >= 2 else 6)
    cm = c[-1][-1]
    ax.axvline(cm, linestyle="--", color="black", linewidth=0.8)
    ax.set_yticks(range(m))
    ax.set_yticklabels([f"$M_{i + 1}$" for i in range(m)])
    ax.invert_yaxis()
    ax.set_xlim(0, cm * 1.04)
    ax.set_xlabel("Tiempo")
    ax.set_title(f"{titulo}  ($C_{{max}}={cm}$)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)


# ---------------------------------------------------------------------------
# Calculos y verificacion
# ---------------------------------------------------------------------------
nombres = [f"J{j + 1}" for j in range(len(P))]
etq = lambda s: " - ".join(nombres[j] for j in s)

print("=== Ejemplo 2 maquinas (Johnson) ===")
a = [x[0] for x in J2]
b = [x[1] for x in J2]
sj = johnson(a, b)
print("Johnson:", " - ".join(NOM_J2[j] for j in sj), "Cmax =", cmax(sj, [list(x) for x in J2]))
print("Optimo exhaustivo:", optimo([list(x) for x in J2])[0],
      "| LB:", cota_inferior([list(x) for x in J2]))

print("\n=== Ejemplo 4 maquinas, 5 trabajos ===")
print("Totales T_j:", [sum(r) for r in P])
print("Cota inferior:", cota_inferior(P))
orig = list(range(len(P)))
print("Orden original:", etq(orig), "Cmax =", cmax(orig, P))
sp, pend = palmer(P)
print("Palmer: pendientes", pend, "->", etq(sp), "Cmax =", cmax(sp, P))
sc, det = cds(P)
for k, s, c in det:
    print(f"  CDS k={k}: {etq(s)} Cmax = {c}")
print("CDS:", etq(sc), "Cmax =", cmax(sc, P))
log = []
sn = neh(P, log)
for entrada in log:
    if entrada[0] == "orden inicial":
        print("NEH orden inicial:", etq(entrada[1]))
    else:
        j, evals, c = entrada
        print(f"  insertar {nombres[j]}:")
        for pos, cand, cc in evals:
            print(f"     pos {pos}: {etq(cand)}  Cmax = {cc}")
print("NEH:", etq(sn), "Cmax =", cmax(sn, P))
copt, sopt = optimo(P)
print("Optimo exhaustivo:", etq(sopt), "Cmax =", copt)
print("Matriz C (NEH):")
for k, j in enumerate(sn):
    print("  ", nombres[j], matriz_c(sn, P)[k])

# ---------------------------------------------------------------------------
# Figura 1: Gantt Johnson (2 maquinas)
# ---------------------------------------------------------------------------
pj = [list(x) for x in J2]
fig, ax = plt.subplots(figsize=(7.2, 2.4))
gantt(ax, sj, pj, NOM_J2, "Secuencia de Johnson")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "gantt_johnson.png"), dpi=200)
plt.close(fig)

# ---------------------------------------------------------------------------
# Figura 2: Gantt comparativo (orden original vs NEH)
# ---------------------------------------------------------------------------
fig, axs = plt.subplots(2, 1, figsize=(7.2, 5.0), sharex=True)
gantt(axs[0], orig, P, nombres, "Orden original " + etq(orig))
gantt(axs[1], sn, P, nombres, "Secuencia NEH " + etq(sn))
xmax = max(cmax(orig, P), cmax(sn, P)) * 1.04
for ax in axs:
    ax.set_xlim(0, xmax)
axs[0].set_xlabel("")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "gantt_comparacion.png"), dpi=200)
plt.close(fig)

# ---------------------------------------------------------------------------
# Figura 3: tamano del espacio de busqueda
# ---------------------------------------------------------------------------
ns = list(range(2, 21))
fig, ax = plt.subplots(figsize=(6.4, 3.6))
ax.plot(ns, [math.log10(math.factorial(n)) for n in ns], marker="o", ms=3,
        label=r"Permutacional: $n!$")
ax.plot(ns, [5 * math.log10(math.factorial(n)) for n in ns], marker="s", ms=3,
        label=r"General, $m=5$: $(n!)^5$")
ax.plot(ns, [10 * math.log10(math.factorial(n)) for n in ns], marker="^", ms=3,
        label=r"General, $m=10$: $(n!)^{10}$")
ax.set_xlabel("Numero de trabajos $n$")
ax.set_ylabel(r"$\log_{10}$ (numero de secuencias)")
ax.set_xticks(range(2, 21, 2))
ax.grid(alpha=0.3)
ax.legend(frameon=False, fontsize=9)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "espacio_busqueda.png"), dpi=200)
plt.close(fig)

print("\nFiguras guardadas en", OUT)