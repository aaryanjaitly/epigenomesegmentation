#!/usr/bin/env bash
set -euo pipefail

# Map positional arguments to variables
SAMPLE_ID=$1
SEGMENT_COORDS=$2
WGBS_PATH=$3
NOME_PATH=$4
HISTONE_PATH=$5

SORTED_SEG="${SAMPLE_ID}_segments_sorted.bed"
bedtools sort -i "${SEGMENT_COORDS}" > "$SORTED_SEG"

CURRENT_BED="$SORTED_SEG"
HAS_WGBS="false"
HAS_NOME="false"

# 1. Process WGBS Map (with Rounding)
if [[ -n "${WGBS_PATH}" ]]; then
    bedtools map \
        -a "$CURRENT_BED" \
        -b <(tail -n +2 "${WGBS_PATH}") \
        -c 4,5 -o median -null 0 \
    | awk -v OFS='\t' '{$4=int($4+0.5); $5=int($5+0.5); print}' \
    > temp_wgbs.bed
    mv temp_wgbs.bed "$CURRENT_BED"
    HAS_WGBS="true"
fi

# 2. Process NOMe Map (with Rounding)
if [[ -n "${NOME_PATH}" ]]; then
    if [[ "$HAS_WGBS" == "true" ]]; then
        # If WGBS ran, NOMe columns are appended as 6 and 7
        bedtools map \
            -a "$CURRENT_BED" \
            -b <(tail -n +2 "${NOME_PATH}") \
            -c 4,5 -o median -null 0 \
        | awk -v OFS='\t' '{$6=int($6+0.5); $7=int($7+0.5); print}' \
        > temp_nome.bed
    else
        # If no WGBS, NOMe columns are 4 and 5
        bedtools map \
            -a "$CURRENT_BED" \
            -b <(tail -n +2 "${NOME_PATH}") \
            -c 4,5 -o median -null 0 \
        | awk -v OFS='\t' '{$4=int($4+0.5); $5=int($5+0.5); print}' \
        > temp_nome.bed
    fi
    mv temp_nome.bed "$CURRENT_BED"
    HAS_NOME="true"
fi

FINAL_OUT="${SAMPLE_ID}_seg_tabs_counts.bed"
INTERSECT_TEMP="intersect_temp.bed"
HISTONE_HEADER=$(head -n 1 "${HISTONE_PATH}")

# 3. Intersect with Histone data
bedtools intersect \
    -a "$CURRENT_BED" \
    -b <(tail -n+2 "${HISTONE_PATH}") \
    -wa -wb \
    > "$INTERSECT_TEMP"

# 4. Construct Table based on available modalities
if [[ "$HAS_WGBS" == "true" && "$HAS_NOME" == "true" ]]; then
    echo -e "${HISTONE_HEADER}\tCov_WGBS\tMeth_WGBS\tCov_NOME\tMeth_NOME" > "$FINAL_OUT"
    awk -v OFS='\t' '{
        for(i=8; i<=NF; i++) printf "%s%s", $i, OFS;
        printf "%s\t%s\t%s\t%s\n", $4, $5, $6, $7;
    }' "$INTERSECT_TEMP" >> "$FINAL_OUT"

elif [[ "$HAS_WGBS" == "true" ]]; then
    echo -e "${HISTONE_HEADER}\tCov_WGBS\tMeth_WGBS" > "$FINAL_OUT"
    awk -v OFS='\t' '{
        for(i=6; i<=NF; i++) printf "%s%s", $i, OFS;
        printf "%s\t%s\n", $4, $5;
    }' "$INTERSECT_TEMP" >> "$FINAL_OUT"

elif [[ "$HAS_NOME" == "true" ]]; then
    echo -e "${HISTONE_HEADER}\tCov_NOME\tMeth_NOME" > "$FINAL_OUT"
    awk -v OFS='\t' '{
        for(i=6; i<=NF; i++) printf "%s%s", $i, OFS;
        printf "%s\t%s\n", $4, $5;
    }' "$INTERSECT_TEMP" >> "$FINAL_OUT"
else
    exit 1
fi

# -----------------------------------------------------------------
# 5. SORT COLUMNS ALPHABETICALLY TO ENSURE CONSISTENCY
# -----------------------------------------------------------------
HEADER=$(head -n 1 "$FINAL_OUT")

NUM_METRICS=0
if [[ "$HAS_WGBS" == "true" ]]; then NUM_METRICS=$((NUM_METRICS + 2)); fi
if [[ "$HAS_NOME" == "true" ]]; then NUM_METRICS=$((NUM_METRICS + 2)); fi

HIST_COLS=$(echo "$HEADER" | awk -F'\t' -v m="$NUM_METRICS" '{
    for(i=4; i<=NF-m; i++) print i, $i
}' | sort -k2,2 | awk '{printf "$%s,", $1}' | sed 's/,$//')

PRINT_COLS='$1, $2, $3'
if [[ -n "$HIST_COLS" ]]; then
    PRINT_COLS="$PRINT_COLS, $HIST_COLS"
fi
if [[ $NUM_METRICS -gt 0 ]]; then
    METRIC_COLS=$(echo "$HEADER" | awk -F'\t' -v m="$NUM_METRICS" '{
        for(i=NF-m+1; i<=NF; i++) printf "$%s,", i
    }' | sed 's/,$//')
    PRINT_COLS="$PRINT_COLS, $METRIC_COLS"
fi

awk -F'\t' -v OFS='\t' '{print '"$PRINT_COLS"'}' "$FINAL_OUT" > temp_sorted.bed
mv temp_sorted.bed "$FINAL_OUT"