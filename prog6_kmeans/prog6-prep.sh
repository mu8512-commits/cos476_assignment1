#!/usr/bin/env bash

# -----------------------------
# 0) Ensure local data dir exists + link dataset
# -----------------------------
mkdir -p data

DATA_SRC="/scratch/network/zl1111/ece476-public/pa1/data.dat"
DATA_DST="data/data.dat"

# Create symlink only if it doesn't already exist
if [ ! -e "$DATA_DST" ]; then
  ln -s "$DATA_SRC" "$DATA_DST"
fi

# -----------------------------
# 1) Load Anaconda module
# -----------------------------
module load anaconda3/2025.12

eval "$(conda shell.bash hook)"   # <-- makes `conda activate` work

conda activate kmeans_p310
