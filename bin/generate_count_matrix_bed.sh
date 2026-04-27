#!/usr/bin/env bash
set -euo pipefail
export LC_ALL=C

# Map positional arguments
SAMPLE_ID=$1
MODALITY=$2
FILE_PATH=$3

mkdir -p "${SAMPLE_ID}_${MODALITY}"
OUTPUT_FILE="${SAMPLE_ID}_${MODALITY}/${SAMPLE_ID}_${MODALITY}.tab"

READ_CMD="cat"
if [[ "${FILE_PATH}" == *.gz ]]; then
    READ_CMD="zcat"
fi

{
    echo -e "chr\tstart\tend\tCov\tMeth"
    
    join -1 1 -2 1 -a 1 -a 2 \
        <( ${READ_CMD} "${FILE_PATH}" \
            | awk -v OFS='\t' '$6=="+" {print $1"_"$2, $1, $2, $6, $5, ($11*$5)/100}' \
            | sort -k1,1 ) \
        <( ${READ_CMD} "${FILE_PATH}" \
            | awk -v OFS='\t' '$6=="-" {print $1"_"($2-1), $1, ($2-1), $6, $5, ($11*$5)/100}' \
            | sort -k1,1 ) \
    | sort -k2,2 -k3,3n \
    | awk -v OFS='\t' '{
        sub(/^chr/, "", $2);
        printf ("%s\t%d\t%d\t%d\t%d\n", $2, $3, $3+1, $5+$10, $6+$11+0.5)
    }' \
    | grep -E -v "random|GL|NC|M|hs|hap|Un|J|EBV|ph|L"

} > "${OUTPUT_FILE}"