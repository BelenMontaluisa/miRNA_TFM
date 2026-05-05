#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 10:34:40 2025

@author: belenmontaluisa
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import matplotlib.pyplot as plt
from Bio import SeqIO
import RNA
import pandas as pd
# =============================================================================
#  Revisión de calidad y filtrado de las secuencias positivas
# =============================================================================

archivos_input = ["hairpin_high_conf_gga_1.fasta", 
                  "hairpin_tg_1.fasta"]

LIMITE_NT = 150
archivos_filtrados = []  # <- para guardar nombres generados

for archivo in archivos_input:
    
    prefix = os.path.splitext(archivo)[0]
    rna_final = f"{prefix}_filtrado.fasta"
    archivos_filtrados.append(rna_final)

    lengths = []
    filtered_records = []
    removed_count = 0

    for record in SeqIO.parse(archivo, "fasta"):
        seq_len = len(record.seq)
        lengths.append(seq_len)
        
        if seq_len <= LIMITE_NT:
            filtered_records.append(record)
        else:
            removed_count += 1

    # --------- HISTOGRAMA DE LONGITUD ---------
    plt.figure(figsize=(8, 5))
    plt.hist(lengths, bins=15)
    plt.xlabel("Longitud (nt)")
    plt.ylabel("Frecuencia")
    plt.title(f"Distribución de longitudes - {prefix}")
    plt.tight_layout()
    plt.show()

    # Guardar archivo filtrado
    SeqIO.write(filtered_records, rna_final, "fasta")

    print(f"--- Reporte para: {archivo} ---")
    print(f"Secuencias conservadas (<= {LIMITE_NT} nt): {len(filtered_records)}")
    print(f"Secuencias eliminadas: {removed_count}")
    print(f"Archivo guardado como: {rna_final}\n")

# =============================================================================
# ARCHIVOS SEPARADOS POR ESPECIE 
# =============================================================================

os.rename("hairpin_high_conf_gga_1_filtrado.fasta", "positivos_gallus.fa")
os.rename("hairpin_tg_1_filtrado.fasta", "positivos_taeniopygia.fa")

# =============================================================================
#  GRÁFICA DE MFE 
# =============================================================================


resultados = []

for record in SeqIO.parse("positivos.fa", "fasta"):
    seq = str(record.seq).replace("T", "U")  # convertir a RNA
    estructura, mfe = RNA.fold(seq)
    
    resultados.append({
        "id": record.id,
        "longitud": len(seq),
        "estructura": estructura,
        "MFE": mfe,
        "MFE_normalizado": mfe/len(seq) 
    })

df = pd.DataFrame(resultados)
df.to_csv("positivos_mfe.csv", index=False)

print("Archivo con MFE generado.")

plt.figure(figsize=(6,6))
plt.boxplot(df["MFE_normalizado"])
plt.ylabel("MFE normalizado")
plt.title("Distribución de MFE normalizado")
plt.axhline(0)
plt.show()