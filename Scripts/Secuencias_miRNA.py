import pandas as pd

# Cargar archivo
df = pd.read_csv("predicciones_diglossa.csv")

# Separar partes del ID
df[["miRNA_full","genomic_info"]] = df["ID"].str.split("::", expand=True)

df[["chr","coords_strand"]] = df["genomic_info"].str.split(":", expand=True)

df[["coords","strand"]] = df["coords_strand"].str.extract(r'([0-9\-]+)\(([\+\-])\)')

df[["start","end"]] = df["coords"].str.split("-", expand=True)

# Limpiar >
df["miRNA_full"] = df["miRNA_full"].str.replace(">", "")

# Convertir a numérico
df["start"] = df["start"].astype(int)
df["end"] = df["end"].astype(int)

# Reordenar columnas
df_final = df[[
    "miRNA_full",
    "chr",
    "start",
    "end",
    "strand",
    "SVM_Score",
    "Is_miRNA"
]]

# Guardar
df_final.to_csv("predicciones_diglossa_separado.csv", index=False)

print("Tabla guardada como predicciones_diglossa_separado.csv")