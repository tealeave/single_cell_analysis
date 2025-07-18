#!/bin/bash
# GSE235063 Dataset Download Script
# Run this as a background job on HPC

# Load required modules
module load python/3.10.2

# Set up environment
export DEBUG_MODE=${DEBUG_MODE:-false}
export TMPDIR=/tmp  # Use local temp for better performance

# Create directories
mkdir -p data/raw data/processed data/references

# Download strategy - use aria2c for parallel downloads if available
echo "Starting GSE235063 dataset download..."
echo "Estimated download time: 2-4 hours for 7GB dataset"

# Check if aria2c is available for faster downloads
if command -v aria2c &> /dev/null; then
    echo "Using aria2c for parallel downloads"
    
    # Download matrix file (largest file)
    aria2c -x 16 -s 16 -c \
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE235063\&format=file\&file=GSE235063%5Fmatrix%2Emtx%2Egz" \
        -o data/raw/GSE235063_matrix.mtx.gz &
    
    # Download barcodes
    aria2c -x 8 -s 8 -c \
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE235063\&format=file\&file=GSE235063%5Fbarcodes%2Etsv%2Egz" \
        -o data/raw/GSE235063_barcodes.tsv.gz &
    
    # Download features
    aria2c -x 8 -s 8 -c \
        "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE235063\&format=file\&file=GSE235063%5Ffeatures%2Etsv%2Egz" \
        -o data/raw/GSE235063_features.tsv.gz &
    
    wait  # Wait for all downloads to complete
    
else
    echo "Using uv python for sequential downloads"
    uv run python src/download_data.py --output data/raw/
fi

# Validate downloads
echo "Validating downloads..."
ls -lh data/raw/

# Create sample manifest
echo "Creating sample manifest..."
uv run python -c "
import os
from pathlib import Path
import pandas as pd

# Create basic manifest for 75 samples
samples = [f'S{i+1}' for i in range(75)]
manifest = pd.DataFrame({
    'sample_id': samples,
    'file_path': [f'data/raw/{s}' for s in samples],
    'status': 'downloading',
    'n_cells': None,
    'n_genes': None
})
manifest.to_csv('data/raw/sample_manifest.csv', index=False)
print('Created sample manifest for 75 samples')
"

echo "Download complete! Files saved to data/raw/"
echo "Next steps:"
echo "1. Run: uv run python src/download_data.py --validate-only"
echo "2. Start debug mode: export DEBUG_MODE=true && uv run nextflow run main.nf --debug"