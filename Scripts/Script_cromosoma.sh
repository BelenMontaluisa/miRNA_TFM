#!/bin/bash

URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/023/105/GCA_019023105.1_LSU_DiBr_2.0/GCA_019023105.1_LSU_DiBr_2.0_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/"
# La variable URL se cambiará dependiendo del genoma del ave de interés. 

# Obtener listado de archivos
curl -s $URL | grep -Eo 'href="[^"]+\.gz"' | cut -d'"' -f2 > lista_archivos_cromosomas.txt
# Se obtiene el filtrado de archivos con grep -E que permite usar expresiones regulares extendidas (para buscar patrones más complejos) y -o extrae sólo la parte del texto que
#coincide con el patrón. 'href="[^"]+\.gz"' 
## se utiliza pipe | para agregar cut -d "', que corta usando como delimitador (") y con -f2 se queda con el segundo campo. 
## Todos los nombres de los archivos se guardan en la lista de archivos de cromosomas.

# Descargar cada archivo
while read -r file; do
  curl -O "${URL}${file}"
done < lista_archivos.txt
# Bucle para descargar cada archivo
## while read -r file permits leer el archivo lista_archivos_cromosomas.txt línea por línea y guarda la linea en la variable file. 
## Se descarga el archivo con curl -O que guarda el archivo con su nombre original. 
## "${URL}${file}" une la URL con el nombre del archivo. 
 

# Descomprimir
gzip -d *.gz

