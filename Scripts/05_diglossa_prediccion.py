import pandas as pd
import numpy as np
import joblib
import RNA
import random
import re
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import product

# ==========================================================
# 1. CARGAR MODELO Y CONFIGURACIÓN
# ==========================================================

svm_model = joblib.load('modelo_SVM_final.pkl')
scaler = joblib.load('scaler_final.pkl')

with open('variables_finales.txt', 'r') as f:
    selected_features = f.read().split(',')

# ==========================================================
# 2. FUNCIONES DE EXTRACCIÓN 
# ==========================================================
SEED = 7
random.seed(SEED)
np.random.seed(SEED)


def get_triplet_map():
    bases = ['A', 'U', 'G', 'C']
    struct_combos = [''.join(p) for p in product(['(', '.'], repeat=3)]
    return [b + s for b in bases for s in struct_combos]

TRIPLET_ORDER = get_triplet_map()
def get_triplet_features(seq, structure):
    """Calcula la frecuencia de tripletes base-estructura."""
    struct_simple = structure.replace(')', '(')
    counts = {k: 0 for k in TRIPLET_ORDER}
    total = 0

    for i in range(1, len(seq) - 1):
        base = seq[i]
        triplet_s = struct_simple[i-1:i+2]
        key = base + triplet_s
        if key in counts:
            counts[key] += 1
            total += 1

    if total == 0:
        return counts

    return {k: v / total for k, v in counts.items()}

def get_kmer_features(seq, k=3):
    kmers = [''.join(p) for p in product("AUGC", repeat=k)]
    counts = {kmer: 0 for kmer in kmers}
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in counts:
            counts[kmer] += 1
    total = sum(counts.values())
    return {k: v / total for k, v in counts.items()} if total > 0 else counts

def get_sequence_features(seq):
    length = len(seq)
    counts = {nt: seq.count(nt) for nt in "AUGC"}
    gc = (counts["G"] + counts["C"]) / length if length > 0 else 0

    return {
        "length": length,
        "A_freq": counts["A"] / length,
        "U_freq": counts["U"] / length,
        "G_freq": counts["G"] / length,
        "C_freq": counts["C"] / length,
        "GC_content": gc,
        "AU_GC_ratio": ((counts["A"] + counts["U"]) / (counts["G"] + counts["C"])) if (counts["G"] + counts["C"]) > 0 else 0
    }
def get_structure_features(structure):
    stems = re.findall(r"\(+", structure)
    loops = re.findall(r"\.+", structure)

    stem_total_length = sum(len(s) for s in stems)
    longest_pairing = max((len(s) for s in stems), default=0)
    terminal_loop_size = max((len(l) for l in loops), default=0)
    bulges = sum(1 for l in loops if 0 < len(l) <= 3)

    left = structure.count("(")
    right = structure.count(")")
    asymmetry = abs(left - right)

    total_dots = structure.count(".")
    stem_loop_ratio = stem_total_length / total_dots if total_dots > 0 else 0

    return {
        "num_base_pairs": left,
        "stem_total_length": stem_total_length,
        "longest_pairing": longest_pairing,
        "terminal_loop_size": terminal_loop_size,
        "num_bulges": bulges,
        "asymmetry": asymmetry,
        "stem_loop_ratio": stem_loop_ratio
    }

def compute_zscore(sequence, real_mfe, n=10):
    shuffled_mfes = []
    for _ in range(n):
        shuffled_seq = "".join(random.sample(sequence, len(sequence)))
        _, mfe = RNA.fold(shuffled_seq)
        shuffled_mfes.append(mfe)
    std_mfe = np.std(shuffled_mfes)
    if std_mfe == 0:
        return 0
    return (real_mfe - np.mean(shuffled_mfes)) / std_mfe

def extract_for_prediction(id_seq, seq, structure):
    seq = seq.upper().replace("T", "U")
    length = len(seq)
    gc = (seq.count("G") + seq.count("C")) / length
    _, mfe = RNA.fold(seq)
    
    # Cálculos básicos
    mfe_norm = mfe / length
    mfei = (mfe_norm / gc) if gc > 0 else 0
    
    # Diccionario base
    row = {
        "mfe": mfe,
        "mfe_norm": mfe_norm,
        "MFEI": mfei,
        "zscore": compute_zscore(seq, mfe, n=20) 
    }
    
    # Unir todas las features (Kmers, Triplets, etc)
    row.update(get_sequence_features(seq))
    row.update(get_structure_features(structure))
    row.update(get_kmer_features(seq, k=3))
    row.update(get_triplet_features(seq, structure))
    
    return row

# ==========================================================
# 3. PROCESAMIENTO DEL ARCHIVO .FOLD
# ==========================================================
input_file = "candidatos_pre.fold"
data_list = []
ids = []

print(f"Abriendo {input_file}...")
try:
    with open(input_file, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    
    print(f"Total de líneas leídas: {len(lines)}")
    
    # Un archivo .fold tiene: 1. ID, 2. Secuencia, 3. Estructura (MFE)
    # Por eso saltamos de 3 en 3
    for i in range(0, len(lines), 3):
        if i + 2 >= len(lines):
            break
            
        header = lines[i]
        sequence = lines[i+1]
        # Ejemplo: "(((.((...))) (-15.0)" -> "(((.((...)))"
        structure = lines[i+2].split()[0]
        
        # Solo procesamos si parece una secuencia de nucleótidos
        if all(c in "ATGCU" for c in sequence.upper()):
            features_dict = extract_for_prediction(header, sequence, structure)
            data_list.append(features_dict)
            ids.append(header)
            if len(ids) % 100 == 0:
                print(f"Procesados {len(ids)} fragmentos...")
        else:
            print(f"Saltando bloque {i} por secuencia no válida")

    if not data_list:
        print("¡ERROR! No se pudo extraer ninguna característica. Revisa el formato de tu archivo .fold")
        exit()

    df_genoma = pd.DataFrame(data_list)
    print(f"DataFrame creado con éxito. Filas: {len(df_genoma)}")

except FileNotFoundError:
    print(f"Error: El archivo {input_file} no existe en esta carpeta.")
    exit()
# ==========================================================
# 4. PREDICCIÓN
# ==========================================================
# 1. Seleccionar las variables que el SVM espera
X_genoma = df_genoma[selected_features]

# 2. Escalar (Muy importante para SVM)
X_scaled = scaler.transform(X_genoma)

# 3. Predecir
predicciones = svm_model.predict(X_scaled)
probabilidades = svm_model.predict_proba(X_scaled)[:, 1]

# 4. Guardar resultados
df_final = pd.DataFrame({
    'ID': ids,
    'SVM_Score': probabilidades,
    'Is_miRNA': predicciones
})

df_final_filtrado = df_final[df_final["SVM_Score"] > 0.8]
df_final_filtrado.to_csv("predicciones_diglossa.csv", index=False)
print("Proceso terminado. Resultados en predicciones_diglossa.csv")

