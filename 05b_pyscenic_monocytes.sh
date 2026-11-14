#!/usr/bin/env bash
# Fig. 4 companion analysis: monocyte regulons for Supplementary Fig. 5F.
# Input: Mono.loom containing raw RNA counts from all monocyte subsets.
# Infer the network jointly; select CM, NCM and Mono_C1QC for the final heatmap.
# Requires pySCENIC and the three reference files listed below.
# Usage: bash 05b_pyscenic_monocytes.sh /path/to/references [Mono.loom] [output_dir]

set -euo pipefail

if [[ $# -lt 1 || $# -gt 3 ]]; then
    echo "Usage: bash $0 REFERENCE_DIR [INPUT_LOOM] [OUTPUT_DIR]" >&2
    exit 1
fi

reference_dir="$1"
input_loom="${2:-Mono.loom}"
output_dir="${3:-scenic_monocytes}"
tfs="$reference_dir/allTFs_hg38.txt"
rankings="$reference_dir/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"
annotations="$reference_dir/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

command -v pyscenic >/dev/null 2>&1 || {
    echo "Activate the environment containing pySCENIC first." >&2
    exit 1
}
for input_file in "$input_loom" "$tfs" "$rankings" "$annotations"; do
    if [[ ! -s "$input_file" ]]; then
        echo "Missing or empty input: $input_file" >&2
        exit 1
    fi
done
mkdir -p "$output_dir"

# 1. Infer TF-target associations using GRNBoost2.
pyscenic grn \
    --num_workers 1 \
    --sparse \
    --method grnboost2 \
    --output "$output_dir/Mono_grn.csv" \
    "$input_loom" "$tfs"

# 2. Identify motif-supported regulons.
pyscenic ctx "$output_dir/Mono_grn.csv" "$rankings" \
    --annotations_fname "$annotations" \
    --expression_mtx_fname "$input_loom" \
    --output "$output_dir/Mono_regulons.csv" \
    --num_workers 3 \
    --mode dask_multiprocessing \
    --mask_dropouts

# 3. Score regulon activity in individual cells.
pyscenic aucell \
    --num_workers 3 \
    --output "$output_dir/Mono_SCENIC.loom" \
    "$input_loom" "$output_dir/Mono_regulons.csv"

# In R, load 05_C1QC_monocyte_analysis.R and use the matching Seurat object:
# sfig5F(sobj, "scenic_monocytes/Mono_SCENIC.loom")
# Cell IDs must match between sobj and the loom file.
echo "Completed: $output_dir/Mono_SCENIC.loom"
