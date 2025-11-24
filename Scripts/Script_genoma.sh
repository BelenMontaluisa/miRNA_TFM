#!/bin/bash
# Variable para guardar el link de descarga del genoma de Diglossa b.

URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/019/023/105/GCA_019023105.1_LSU_DiBr_2.0/"

# Se obtiene el filtrado de archivos con grep -E permite usar expresiones regulares extendidas (para buscar patrones más complejos) y -o extrae sólo la parte del texto 
# que coincide con el patrón. 'href="[^"]+\.gz"' 
# se agrega cut -d "', que corta usando como delimitador (") y con -f2 se queda con el segundo campo. 
# Todos los nombres de los archivos se guardan en la lista de archivos de cromosomas.

curl -s $URL | grep -Eo 'href="[^"]+\.gz"' | cut -d'"' -f2 > lista_archivos_cromosomas.txt


# Bucle para descargar cada archivo
## while read -r file permits leer el archive lista_archivos_cromosomas.txt línea por línea y guarda la linea en la variable file. 
## Se descarga el archivo con curl -O guarda el archivo con su nombre original. 
## "${URL}${file}" une la URL con el nombre del archivo. 
 
while read -r file; do
  curl -O "${URL}${file}"
done < lista_archivos.txt


# Descomprimir con gzip para todos los archivos con extensión .gz
gzip -d *.gz

