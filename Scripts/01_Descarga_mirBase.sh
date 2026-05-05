# Script para la descarga de secuencias de mirBase

#!/bin/bash
# =============================================================================
#  Descarga de secuencias de miRBase
# =============================================================================

# La variable URL base del FTP de miRBase. (HTML con links)
URL="https://www.mirbase.org/download/CURRENT/"

# Obtener listado de archivos fasta
curl -s $URL | grep -Eo 'href="[^"]+\.fa"' | sed 's/href="//; s/$//' > lista_archivos_fasta.txt

# Descargar cada archivo
while read -r file; do  
  if [[ "$file" != http* ]]; then
    full_url="${URL}${file}"
  else
    full_url="$file"
  fi


  curl -O "$full_URL"
done < lista_archivos_fasta.txt

#Filtrado secuencias de Gallus gallus del archivo de miRNA de alta confianza 
awk '/^>/ { f = match($0, /Gallus gallus/) } f' hairpin_high_conf.fasta > hairpin_high_conf_gga_1.fasta\n
grep -o '>' hairpin_high_conf_gga.fasta | wc -l
head hairpin_high_conf_gga_1.fasta


#Filtrado de secuencias Taeniogypia guttata el archivo de miRNA hairpin. 
awk '/^>/ { f = match($0, /Taeniopygia guttata/) } f' hairpin.fasta > hairpin_tg_1.fasta\n
grep -o '>' hairpin_tg_1.fasta | wc -l
head hairpin_tg_1.fasta

