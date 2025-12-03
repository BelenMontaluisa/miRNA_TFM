#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 10:34:40 2025

@author: belenmontaluisa
"""

from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt

seqs = [str(rec.seq) for rec in SeqIO.parse("hairpin_tg_1.fasta", "fasta")]

# --------- Longitudes ---------
lengths = [len(s) for s in seqs]

# --------- Conteo de bases ---------
counts = Counter("".join(seqs))
total = sum(counts[b] for b in "AUGC")
percent = {b: counts[b] / total * 100 for b in "AUGC"}
percent_rounded = {b: round(percent[b], 1) for b in percent}
print(percent_rounded)


# --------- HISTOGRAMA DE LONGITUD ---------
palette = [
    "steelblue", "orange", "purple", "tomato", "mediumseagreen",
    "goldenrod", "orchid", "slateblue", "darkcyan", "salmon"
]
plt.figure(figsize=(10, 6))
n, bins, patches = plt.hist(lengths, bins=10)

for patch, c in zip(patches, palette * 6):
    patch.set_facecolor(c)

plt.xlabel("Longitud (nt)")
plt.ylabel("Frecuencia")
plt.title("Distribución de longitudes de pre-miRNA en Taeniopygia guttata")
plt.show()


# =============================================================================
# Filtrado de secuencias
# =============================================================================
# Archivo de entrada
rna_inicial = "hairpin_tg_1.fasta"

# Archivo de salida
rna_final = "hairpin_tg_filtrado.fasta"
filtered_records = []
removed_records = []

# Procesar secuencias
for record in SeqIO.parse(rna_inicial, "fasta"):
    if len(record.seq) <= 120:
        filtered_records.append(record)
    else:
        removed_records.append(record)

# Guardar archivo filtrado
SeqIO.write(filtered_records, rna_final, "fasta")

# Reporte
print("Total de secuencias originales:", len(filtered_records) + len(removed_records))
print("Secuencias conservadas (<=120 nt):", len(filtered_records))
print("Secuencias eliminadas (>120 nt):", len(removed_records))

# (Opcional) Mostrar IDs eliminados
print("\nIDs eliminados (longitud > 120 nt):")
for r in removed_records:
    print(r.id, len(r.seq))
