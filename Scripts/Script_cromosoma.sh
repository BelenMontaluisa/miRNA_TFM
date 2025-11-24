#!/bin/bash

URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/023/105/GCA_019023105.1_LSU_DiBr_2.0/GCA_019023105.1_LSU_DiBr_2.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/"
# La variable URL se cambiará dependiendo del genoma del ave de interés. 
# Obtener listado de archivos
curl -s $URL | grep -Eo 'href="[^"]+\.gz"' | cut -d'"' -f2 > lista_archivos.txt

# Descargar cada archivo
while read -r file; do
  curl -O "${URL}${file}"
done < lista_archivos.txt

# Descomprimir
gzip -d *.gz

