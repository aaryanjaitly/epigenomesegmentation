#!/usr/bin/env bash
set -euo pipefail

# Map positional arguments
GENOME=$1

DOWNLOAD_FILE="downloaded.chrom.sizes"
TMP_FILE="temp.chrom.sizes"
FINAL_FILE="${GENOME}.chrom.sizes"

URL="http://hgdownload.soe.ucsc.edu/goldenPath/${GENOME}/bigZips/${FINAL_FILE}"

curl -s -o "${DOWNLOAD_FILE}" "${URL}"

if [[ "${GENOME}" == hg* ]]; then
    awk 'BEGIN{OFS="\t"} {
        sub(/^chr/, "", $1);
        if ($1 ~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|X|Y)$/) print
    }' "${DOWNLOAD_FILE}" > "${TMP_FILE}"

elif [[ "${GENOME}" == mm* ]]; then
    awk 'BEGIN{OFS="\t"} {
        sub(/^chr/, "", $1);
        if ($1 ~ /^(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|X|Y)$/) print
    }' "${DOWNLOAD_FILE}" > "${TMP_FILE}"

else
    awk 'BEGIN{OFS="\t"} { 
        sub(/^chr/, "", $1); 
        print 
    }' "${DOWNLOAD_FILE}" > "${TMP_FILE}"
fi

sort -k1,1V "${TMP_FILE}" > "${FINAL_FILE}"

rm "${DOWNLOAD_FILE}" "${TMP_FILE}"