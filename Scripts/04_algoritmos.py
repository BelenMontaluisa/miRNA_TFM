# ==========================================================
# TFM - Clasificación de pre-miRNA
# Random Forest vs SVM
# ==========================================================

import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    confusion_matrix,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    roc_curve,
    balanced_accuracy_score
)

# ==========================================================
# 1. Selector Mann-Whitney integrado para Nested CV
# ==========================================================
class MannWhitneySelector(BaseEstimator, TransformerMixin):
    def __init__(self, alpha=0.05, exclude_features=None):
        self.alpha = alpha
        self.exclude_features = exclude_features if exclude_features else []

    def fit(self, X, y):
        X = pd.DataFrame(X)
        y = pd.Series(y).reset_index(drop=True)
        X = X.reset_index(drop=True)

        # Eliminar features constantes
        X = X.loc[:, X.var() > 0]

        features = X.columns
        p_values = []

        for feat in features:
            pos = X.loc[y == 1, feat]
            neg = X.loc[y == 0, feat]

            if len(pos) < 5 or len(neg) < 5:
                p = 1.0
            else:
                _, p = mannwhitneyu(pos, neg, alternative="two-sided")

            p_values.append(p)

        _, p_adj, _, _ = multipletests(p_values, method="fdr_bh")

        self.selected_features_ = [
            feat for feat, p in zip(features, p_adj)
            if p < self.alpha and feat not in self.exclude_features
        ]

        if len(self.selected_features_) == 0:
            self.selected_features_ = list(features)

        return self

    def transform(self, X):
        X = pd.DataFrame(X)
        return X[self.selected_features_]

# ==========================================================
# 2. Estructura de carpetas
# ==========================================================
base_dir = "resultados_alg_mirna"
tablas_dir = os.path.join(base_dir, "tablas")
graficos_dir = os.path.join(base_dir, "graficos")
modelos_dir = os.path.join(base_dir, "modelos")
validacion_dir = os.path.join(base_dir, "validacion")

for d in [base_dir, tablas_dir, graficos_dir, modelos_dir, validacion_dir]:
    os.makedirs(d, exist_ok=True)

# ==========================================================
# 3. Carga y limpieza de datos
# ==========================================================
df_train_balanced = pd.read_csv("dataset_train_gallus_balanced.csv")
df_train_balanced["label"] = df_train_balanced["label"].astype(int)

for col in df_train_balanced.columns:
    if col not in ["id", "label"]:
        df_train_balanced[col] = pd.to_numeric(df_train_balanced[col], errors="coerce")

# Definir X y y 
X = df_train_balanced.drop(columns=["id", "label"])
y = df_train_balanced["label"]

# Asignar variables de entrenamiento
X_train, y_train = X, y

print("Total secuencias finales:", len(df_train_balanced))

# ==========================================================
# 4. Carga del test set
# ==========================================================
df_test = pd.read_csv("dataset_test_taeniopygia.csv")
y_test = df_test["label"].astype(int)
X_test = df_test.drop(columns=["id", "label"]).copy()

# ==========================================================
# 5. ANÁLISIS ESTRUCTURAL PREVIO (MFE, MFEI, Z-SCORE)
# ==========================================================
struct_features = ["mfe_norm", "MFEI", "zscore"]
analysis_results = []

for feat in struct_features:
    if feat not in X_train.columns:
        print(f"{feat} no está en el dataset")
        continue

    pos = X_train.loc[y_train == 1, feat]
    neg = X_train.loc[y_train == 0, feat]

    stat, p = mannwhitneyu(pos, neg, alternative="two-sided")
    mean_diff = pos.mean() - neg.mean()
    pooled_std = np.sqrt((pos.var(ddof=1) + neg.var(ddof=1)) / 2)
    cohen_d = 0 if pooled_std == 0 else mean_diff / pooled_std

    analysis_results.append({
        "feature": feat,
        "mean_pos": pos.mean(),
        "mean_neg": neg.mean(),
        "p_value": p,
        "cohen_d": cohen_d
    })

    plt.figure(figsize=(6,4))
    plt.hist(pos, bins=30, alpha=0.6, label="Positivos")
    plt.hist(neg, bins=30, alpha=0.6, label="Negativos")
    plt.title(f"Distribución {feat}")
    plt.legend()
    plt.show()

analysis_df = pd.DataFrame(analysis_results)
analysis_df["p_adj"] = multipletests(analysis_df["p_value"], method="fdr_bh")[1]
analysis_df.to_csv(os.path.join(tablas_dir, "Analisis_Estructural_MFE_MFEI_Zscore.csv"), index=False)
print("\nRESULTADOS ESTRUCTURALES:")
print(analysis_df)

# ==========================================================
# 6. Selección de variables 
# ==========================================================
features = X_train.columns.tolist()
p_values = []

for feat in features:
    pos = X_train.loc[y_train == 1, feat]
    neg = X_train.loc[y_train == 0, feat]
    _, p = mannwhitneyu(pos, neg, alternative="two-sided")
    p_values.append(p)

p_adj = multipletests(p_values, method="fdr_bh")[1]
mw_results = pd.DataFrame({"feature": features, "p_value": p_values, "p_adj": p_adj})
mw_results.to_csv(os.path.join(tablas_dir, "MannWhitney_results.csv"), index=False)

selected_features = mw_results[mw_results["p_adj"] < 0.05]["feature"].tolist()
selected_features = [f for f in selected_features if f not in ["mfe", "length"]]

if len(selected_features) == 0:
    print("No se encontraron variables significativas. Seleccionando las 10 mejores.")
    selected_features = mw_results.nsmallest(10, "p_adj")["feature"].tolist()

X_train_sel = X_train[selected_features].copy()
X_test_sel = X_test[selected_features].copy()
print("Número de variables seleccionadas:", len(selected_features))

# ==========================================================
# 7. RANDOM FOREST
# ==========================================================
rf = RandomForestClassifier(n_estimators=500, class_weight="balanced", random_state=42)
rf.fit(X_train_sel, y_train)
y_pred_rf = rf.predict(X_test_sel)
y_prob_rf = rf.predict_proba(X_test_sel)[:, 1]

# ==========================================================
# 8. SVM
# ==========================================================
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train_sel)
X_test_scaled = scaler.transform(X_test_sel)

svm = SVC(kernel="rbf", probability=True, class_weight="balanced", random_state=42)
svm.fit(X_train_scaled, y_train)
y_pred_svm = svm.predict(X_test_scaled)
y_prob_svm = svm.predict_proba(X_test_scaled)[:, 1]

# ==========================================================
# 9. Evaluación
# ==========================================================
def evaluar(y_true, y_pred, y_prob):
    cm = confusion_matrix(y_true, y_pred)
    return {
        "Accuracy": accuracy_score(y_true, y_pred),
        "Balanced_Accuracy": balanced_accuracy_score(y_true, y_pred),
        "Precision": precision_score(y_true, y_pred),
        "Recall": recall_score(y_true, y_pred),
        "F1": f1_score(y_true, y_pred),
        "AUC": roc_auc_score(y_true, y_prob),
        "Confusion_Matrix": cm
    }

rf_metrics = evaluar(y_test, y_pred_rf, y_prob_rf)
svm_metrics = evaluar(y_test, y_pred_svm, y_prob_svm)

# ==========================================================
# 10. Nested Cross-Validation
# ==========================================================
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

rf_pipeline = Pipeline([
    ("selector", MannWhitneySelector(alpha=0.05, exclude_features=["mfe","length"])),
    ("model", RandomForestClassifier(n_estimators=500, class_weight="balanced", random_state=42))
])

svm_pipeline = Pipeline([
    ("selector", MannWhitneySelector(alpha=0.05, exclude_features=["mfe","length"])),
    ("scaler", StandardScaler()),
    ("model", SVC(kernel="rbf", probability=True, class_weight="balanced", random_state=42))
])

cv_results_rf = cross_validate(rf_pipeline, X, y, cv=cv,
                               scoring=["roc_auc", "accuracy", "f1", "balanced_accuracy"])
cv_results_svm = cross_validate(svm_pipeline, X, y, cv=cv,
                                scoring=["roc_auc", "accuracy", "f1", "balanced_accuracy"])

cv_summary = pd.DataFrame({
    "Modelo": ["RandomForest", "SVM"],
    "AUC_mean": [cv_results_rf["test_roc_auc"].mean(), cv_results_svm["test_roc_auc"].mean()],
    "AUC_std": [cv_results_rf["test_roc_auc"].std(), cv_results_svm["test_roc_auc"].std()]
})
cv_summary.to_csv(os.path.join(validacion_dir, "Validacion_Cruzada_Resumen.csv"), index=False)

# ==========================================================
# 11. Guardar modelos finales
# ==========================================================
# Crear df para gráficos
df_plot = X_test_sel.copy()
df_plot["Clase"] = y_test.map({1:"Positivo", 0:"Negativo"})

# Tabla comparativa de métricas
comparison = pd.DataFrame({
    "Metric": list(rf_metrics.keys())[:-1],
    "RandomForest": list(rf_metrics.values())[:-1],
    "SVM": list(svm_metrics.values())[:-1]
})
comparison.to_csv(os.path.join(tablas_dir, "Comparacion_Modelos.csv"), index=False)

# Histogramas de MFE normalizado
if "mfe_norm" in df_plot.columns:
    plt.figure(figsize=(6,5))
    sns.histplot(data=df_plot, x="mfe_norm", hue="Clase", bins=30, kde=True, palette="Set2", alpha=0.6)
    plt.title("Distribución de mfe_norm")
    plt.xlabel("mfe_norm")
    plt.ylabel("Frecuencia")
    plt.tight_layout()
    plt.savefig(os.path.join(graficos_dir, "Distribucion_mfe_norm.png"), dpi=300)
    plt.close()

# Boxplots de features seleccionadas
features_visualizar = ["mfe", "mfe_norm", "GC_content"]

for feat in features_visualizar:
    if feat in df_plot.columns:
        plt.figure(figsize=(6,5))
        sns.boxplot(data=df_plot, x="Clase", y=feat, palette="Set2")
        plt.title(f"Diferencia estadística en {feat}")
        plt.tight_layout()
        plt.savefig(os.path.join(graficos_dir, f"Boxplot_{feat}.png"), dpi=300)
        plt.close()

# Guardar matrices de confusión
pd.DataFrame(rf_metrics["Confusion_Matrix"]).to_csv(os.path.join(tablas_dir, "CM_RF.csv"), index=False)
pd.DataFrame(svm_metrics["Confusion_Matrix"]).to_csv(os.path.join(tablas_dir, "CM_SVM.csv"), index=False)

# Guardar modelos
joblib.dump(rf, os.path.join(modelos_dir, "modelo_RF_final.pkl"))
joblib.dump(svm, os.path.join(modelos_dir, "modelo_SVM_final.pkl"))
joblib.dump(scaler, os.path.join(modelos_dir, "scaler_final.pkl"))

# Guardar variables seleccionadas
with open(os.path.join(modelos_dir, "variables_finales.txt"), "w") as f:
    f.write(",".join(selected_features))

# ==========================================================
# 12. Importancia de Variables
# ==========================================================
result_rf = permutation_importance(rf, X_test_sel, y_test, n_repeats=10, random_state=42, n_jobs=-1)
result_svm = permutation_importance(svm, X_test_scaled, y_test, n_repeats=10, random_state=42, n_jobs=-1)

importancia_df = pd.DataFrame({
    'Feature': selected_features,
    'RF_Importance': result_rf.importances_mean,
    'SVM_Importance': result_svm.importances_mean
}).sort_values(by='RF_Importance', ascending=False)

importancia_df.to_csv(os.path.join(tablas_dir, "importancia_variables.csv"), index=False)

plt.figure(figsize=(10, 8))
sns.barplot(x='RF_Importance', y='Feature', data=importancia_df, palette='viridis')
plt.title('Importancia de las Variables - Random Forest (Permutation Importance)')
plt.xlabel('Caída en el desempeño del modelo al permutar')
plt.ylabel('Variables (Features)')
plt.tight_layout()
plt.savefig(os.path.join(graficos_dir, "Importancia_Variables_RF.png"))
plt.close()

# ==========================================================
# 13. Gráficas ROC comparativa
# ==========================================================
plt.figure(figsize=(6,6))
fpr_rf, tpr_rf, _ = roc_curve(y_test, y_prob_rf)
fpr_svm, tpr_svm, _ = roc_curve(y_test, y_prob_svm)

plt.plot(fpr_rf, tpr_rf, label=f"RF (AUC={rf_metrics['AUC']:.2f})")
plt.plot(fpr_svm, tpr_svm, label=f"SVM (AUC={svm_metrics['AUC']:.2f})")
plt.plot([0,1],[0,1],'--', color='gray')
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("Curva ROC Comparativa")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(graficos_dir, "ROC_Comparativa.png"), dpi=300)
plt.close()

# Importancia de variables basada en Random Forest
importances = rf.feature_importances_
indices = np.argsort(importances)[::-1]
features_sorted = [selected_features[i] for i in indices]

plt.figure(figsize=(8,6))
sns.barplot(x=importances[indices], y=features_sorted)
plt.title("Importancia de Variables - Random Forest")
plt.tight_layout()
plt.savefig(os.path.join(graficos_dir, "Importancia_RF.png"), dpi=300)
plt.close()