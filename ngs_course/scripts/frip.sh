#!/bin/bash

OPTIND=1 # Reset in case getopts has been previously used in the cell
verbose=0

sample=""
bam_file=""
peak_file=""

while getopts "h:?:s:b:p:o:" opt; do
        case "$opt" in
        h|\?)
          	echo "Usage: frip.sh -s <sample> -b <bam_file> -p <peak_file>"
		echo "Example: frip.sh -s sample1 -b sample1.bam -p sample1.narrowPeak"
                exit 0
                ;;
        s) sample=$OPTARG
           ;;
	b) bam_file=$OPTARG
           ;;
	p) peak_file=$OPTARG
           ;;
	esac
done

shift $((OPTIND-1))

[ "${1:-}" = "--" ] && shift

if [[ -z "$sample" || -z "$bam_file" || -z "$peak_file" ]]; then
    echo "Error: missing required arguments"
    echo "Usage: $0 -s <sample> -b <bam_file> -p <peak_file>"
    exit 1
fi


echo "Sample: ${sample}"
echo "BAM file: ${bam_file}"
echo "Peak file: ${peak_file}"
echo ""
echo ""

calculate_frip() {

  echo "Processing sample: ${sample}"

  # total reads
  total_reads=$(samtools view -c "${bam_file}")

  # reads in peaks
  reads_in_peaks=$(bedtools sort -i "${peak_file}" \
    | bedtools merge -i stdin | bedtools intersect -u -nonamecheck \
    -a "${bam_file}" -b stdin -ubam | samtools view -c)

  # FRiP score
  FRiP=$(awk "BEGIN {print ${reads_in_peaks}/${total_reads}}")

  # Imprimir el resultado
  echo "FRiP score for ${sample}: ${FRiP}"
  echo "--------------------------"
  
  # Salvar en un archivo
  echo "FRiP score for ${sample}: ${FRiP}" >> "FRiP_summary.txt"
  echo "--------------------------" >> "FRiP_summary.txt"

}
  
calculate_frip
