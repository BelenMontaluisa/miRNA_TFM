#!/bin/bash
# ------------------------------
    ### LIMPIEZA DEL GENOMA
# ------------------------------
# Descripción
# Este script procesa los genomas y anotaciones GTF:
# Filtra sólo los cromosomas nucleares, quitando el genoma mitocondrial y los scaffolds no alineados a ningún cromosoma. 
# Indexar el genoma
# Generar archivos BED de genes, exones, regiones codificantes, 
# Después se filtra 
# Entrada: Genoma en archivos fasta (dna.toplevel.fa) de las especies Gallus gallus y Taeniopygia gutatta
# Salida: Carpeta de limpieza_resultados en donde se encuentran los archivos 
# 1. Archivo de texto plano con la lista de cromosomas ({prefix}_croms.txt )
# 2. Cromosomas nucleares (.nucleares.fa)
# 3. Indice del fasta (.nucleares.fa.fai y genome.sizes) para acceder a regiones específicas, 
# indicando cuántos cromosomas existen y la longitud. 

# ------------------------------
    ### VARIABLES
# ------------------------------

OUTDIR="limpieza_resultados"
mkdir -p "$OUTDIR"
echo -e "Especie\tIntrones\tExones_nc\tIntergenicas" > resumen_global.tsv
MINLEN=60
MAXLEN=140
GC_MIN=0.40
GC_MAX=0.60
MAX_N=0.05

# Loop por todos los GTF
for gtf in *.gtf; do

    prefix="${gtf%%.*}"
    echo "Procesando $prefix ..."

    fasta=$(ls ${prefix}*.dna.toplevel.fa 2>/dev/null | head -n 1)

    if [[ -z "$fasta" ]]; then
        echo "No se encontró FASTA para $prefix"
        continue
    fi

    # ------------------------------
    # 1. Cromosomas válidos
    # ------------------------------
    # Se extrae el encabezado del archivo y se filtra por parámetros (0-9,abecedario y los cromosomas sexuales, Z y W)
    grep "^>" "$fasta" | \
        sed 's/^>//' | \
        cut -d' ' -f1 | \
        grep -xE '([0-9]+[A-Z]?|Z|W)' | \
        sort -u \
        > "$OUTDIR/${prefix}_croms.txt"

    # ------------------------------
    # 2. FASTA nuclear 
    # ------------------------------
    # se compara y extrae secuencias del archivo fasta con los nombres de cromosomas guardados anteriormente
    # Se obtiene archivos fasta solo con el genoma nuclear 
    seqkit grep -f "$OUTDIR/${prefix}_croms.txt" "$fasta" \
        -o "$OUTDIR/${prefix}.cromosomas.fa"

    seqkit grep -v -p "MT" \
        "$OUTDIR/${prefix}.cromosomas.fa" \
        -o "$OUTDIR/${prefix}.nucleares.fa"

    # ------------------------------
    # 3. Indexación del genoma
    # ------------------------------
    # Se indexa para obtener un archivo con la indexación del genoma, con el nombre del cromosoma, la longitud. 
    samtools faidx "$OUTDIR/${prefix}.nucleares.fa"

    cut -f1,2 "$OUTDIR/${prefix}.nucleares.fa.fai" | \
        sort -k1,1 \
        > "$OUTDIR/${prefix}.genome.sizes"

    # ------------------------------
    # 4. Extracción de genes
    # ------------------------------
    # Se extrae los genes desde el archivo GTF, con la información del cromosoma ($1) e inicio y final ($4 y $5). 
    # Debido a que los archivos GTF son 1-based se pasan a 0-based para evitar errores de desplazamiento cuando se utilice Bedtools. 
    grep -v "^#" "$gtf" | \
        awk 'BEGIN{OFS="\t"} $3=="gene" {print $1,$4-1,$5}' \
        > "$OUTDIR/${prefix}.genes.bed"

    awk 'NR==FNR{a[$1];next} $1 in a' \
        "$OUTDIR/${prefix}.genome.sizes" \
        "$OUTDIR/${prefix}.genes.bed" \
        > "$OUTDIR/${prefix}.genes.filtered.bed"
    # Filtra cromosomas a través de genome.sizes y genes.bed, eliminando lo restante para mantener coherencia con el archivo FASTA.  
    bedtools sort \
        -i "$OUTDIR/${prefix}.genes.filtered.bed" \
        -g "$OUTDIR/${prefix}.genome.sizes" \
        > "$OUTDIR/${prefix}.genes.sorted.bed"

    # ------------------------------
    # 5. Intergénicas
    # ------------------------------
    # Se extraen las regiones intergénicas utilizando el bedtools complement para obtener las regiones donde no hay genes. 
    bedtools complement \
        -i "$OUTDIR/${prefix}.genes.sorted.bed" \
        -g "$OUTDIR/${prefix}.genome.sizes" \
        > "$OUTDIR/${prefix}.intergenic.bed"

    # ------------------------------
    # 6. EXONES
    # ------------------------------
    # Se filtra por exones y se agrega la información de los atributos del gen al que pertenece ($9)
    grep -v "^#" "$gtf" | \
        awk 'BEGIN{OFS="\t"} $3=="exon" {print $1,$4-1,$5,$9}' \
        > "$OUTDIR/${prefix}.exons.bed"

    awk 'NR==FNR{a[$1];next} $1 in a' \
        "$OUTDIR/${prefix}.genome.sizes" \
        "$OUTDIR/${prefix}.exons.bed" \
        > "$OUTDIR/${prefix}.exons.filtered.bed"

    bedtools sort \
        -i "$OUTDIR/${prefix}.exons.filtered.bed" \
        -g "$OUTDIR/${prefix}.genome.sizes" \
        > "$OUTDIR/${prefix}.exons.sorted.bed"

    # ------------------------------
    # 7. CDS
    # ------------------------------
    # Se filtra por regiones codificantes, CDS. 

    grep -v "^#" "$gtf" | \
        awk 'BEGIN{OFS="\t"} $3=="CDS" {print $1,$4-1,$5}' \
        > "$OUTDIR/${prefix}.cds.bed"

    awk 'NR==FNR{a[$1];next} $1 in a' \
        "$OUTDIR/${prefix}.genome.sizes" \
        "$OUTDIR/${prefix}.cds.bed" \
        > "$OUTDIR/${prefix}.cds.filtered.bed"

    bedtools sort \
        -i "$OUTDIR/${prefix}.cds.filtered.bed" \
        -g "$OUTDIR/${prefix}.genome.sizes" \
        > "$OUTDIR/${prefix}.cds.sorted.bed"

    # ------------------------------
    # 8. Exones no codificantes
    # ------------------------------
    bedtools subtract \
        -a "$OUTDIR/${prefix}.exons.sorted.bed" \
        -b "$OUTDIR/${prefix}.cds.sorted.bed" \
        > "$OUTDIR/${prefix}.exons_nc.bed"


    # ------------------------------
    # 9. Intrones 
    # ------------------------------

    # Extraer transcripts
    grep -v "^#" "$gtf" | \
        awk 'BEGIN{OFS="\t"} $3=="transcript" {print $1,$4-1,$5,$9}' \
        > "$OUTDIR/${prefix}.transcripts.bed"

    awk 'NR==FNR{a[$1];next} $1 in a' \
        "$OUTDIR/${prefix}.genome.sizes" \
        "$OUTDIR/${prefix}.transcripts.bed" \
        > "$OUTDIR/${prefix}.transcripts.filtered.bed"

    bedtools sort \
        -i "$OUTDIR/${prefix}.transcripts.filtered.bed" \
        -g "$OUTDIR/${prefix}.genome.sizes" \
        > "$OUTDIR/${prefix}.transcripts.sorted.bed"

    bedtools subtract \
    -a "$OUTDIR/${prefix}.transcripts.sorted.bed" \
    -b "$OUTDIR/${prefix}.exons.sorted.bed" \
    > "$OUTDIR/${prefix}.introns.bed"

    # ------------------------------
    # 10. Resumen por especie
    # ------------------------------
    introns=$(wc -l < "$OUTDIR/${prefix}.introns.bed")
    exons_nc=$(wc -l < "$OUTDIR/${prefix}.exons_nc.bed")
    intergenic=$(wc -l < "$OUTDIR/${prefix}.intergenic.bed")

    echo -e "${prefix}\t${introns}\t${exons_nc}\t${intergenic}" >> resumen_global.tsv

done

echo "Proceso completado."
