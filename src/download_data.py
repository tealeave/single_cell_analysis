"""Download and prepare GSE235063 dataset from GEO."""

import argparse
import gzip
import logging
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple

import pandas as pd
import requests
from tqdm import tqdm

from config import RAW_DATA_DIR
from utils import setup_logging


def download_geo_study(study_id: str, output_dir: Path) -> List[Path]:
    """Download GSE235063 study from GEO."""
    
    # GSE235063 specific URLs and files
    base_url = "https://www.ncbi.nlm.nih.gov/geo/download/"
    study_files = [
        f"?acc={study_id}&format=file",
        f"?acc={study_id}&format=file&file=GSE235063%5Fmatrix%2Emtx%2Egz",
        f"?acc={study_id}&format=file&file=GSE235063%5Fbarcodes%2Etsv%2Egz",
        f"?acc={study_id}&format=file&file=GSE235063%5Ffeatures%2Etsv%2Egz"
    ]
    
    downloaded_files = []
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for file_path in study_files:
        url = base_url + file_path
        filename = file_path.split('=')[-1].replace('%2E', '.').replace('%5F', '_')
        local_path = output_dir / filename
        
        logging.info(f"Downloading {filename}...")
        
        try:
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            total_size = int(response.headers.get('content-length', 0))
            
            with open(local_path, 'wb') as f, tqdm(
                desc=filename,
                total=total_size,
                unit='B',
                unit_scale=True,
                unit_divisor=1024,
            ) as pbar:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))
            
            downloaded_files.append(local_path)
            logging.info(f"Downloaded {local_path} ({local_path.stat().st_size / 1024**2:.1f} MB)")
            
        except Exception as e:
            logging.error(f"Failed to download {url}: {e}")
    
    return downloaded_files


def extract_sample_data(geo_dir: Path, output_dir: Path) -> List[str]:
    """Extract individual sample data from GEO matrix files."""
    
    # Expected files from GSE235063
    matrix_file = geo_dir / "GSE235063_matrix.mtx.gz"
    barcodes_file = geo_dir / "GSE235063_barcodes.tsv.gz"
    features_file = geo_dir / "GSE235063_features.tsv.gz"
    
    if not all([matrix_file.exists(), barcodes_file.exists(), features_file.exists()]):
        raise FileNotFoundError("Required GEO files not found")
    
    # Create sample directories
    samples = []
    
    # Read barcodes to identify samples
    logging.info("Reading barcodes...")
    barcodes_df = pd.read_csv(barcodes_file, sep='\t', compression='gzip', header=None)
    
    # Extract sample names from barcodes
    # Format: AAACCTGAGAAACGCC-1-1_S1_AAACCTGAGAAACGCC-1-1
    sample_names = barcodes_df[0].str.extract(r'-([A-Z]\d+)_')[0].unique()
    
    logging.info(f"Found {len(sample_names)} samples: {sample_names}")
    
    # For now, we'll create symbolic links to the matrix files
    # Each sample will use the same matrix but filtered barcodes
    for sample in sample_names:
        sample_dir = output_dir / sample
        sample_dir.mkdir(exist_ok=True)
        
        # Create symbolic links
        (sample_dir / "matrix.mtx.gz").symlink_to(matrix_file.resolve())
        (sample_dir / "features.tsv.gz").symlink_to(features_file.resolve())
        
        # Filter barcodes for this sample
        sample_barcodes = barcodes_df[barcodes_df[0].str.contains(f"{sample}_")]
        sample_barcodes_file = sample_dir / "barcodes.tsv.gz"
        sample_barcodes.to_csv(sample_barcodes_file, sep='\t', compression='gzip', 
                              index=False, header=False)
        
        samples.append(sample)
        logging.info(f"Created sample directory: {sample}")
    
    return samples


def validate_downloads(data_dir: Path) -> bool:
    """Validate downloaded data files."""
    
    required_files = [
        "GSE235063_matrix.mtx.gz",
        "GSE235063_barcodes.tsv.gz", 
        "GSE235063_features.tsv.gz"
    ]
    
    for file_name in required_files:
        file_path = data_dir / file_name
        if not file_path.exists():
            logging.error(f"Missing required file: {file_path}")
            return False
        
        # Check file size
        size_mb = file_path.stat().st_size / 1024**2
        if size_mb < 1:
            logging.warning(f"File {file_name} seems too small ({size_mb:.1f} MB)")
    
    logging.info("All required files present")
    return True


def create_sample_metadata(samples: List[str], output_dir: Path) -> Path:
    """Create sample metadata file."""
    
    metadata = []
    for sample in samples:
        metadata.append({
            'sample_id': sample,
            'file_path': str(output_dir / sample),
            'n_cells': None,  # Will be populated later
            'status': 'raw'
        })
    
    metadata_df = pd.DataFrame(metadata)
    metadata_file = output_dir / "sample_metadata.csv"
    metadata_df.to_csv(metadata_file, index=False)
    
    logging.info(f"Created sample metadata: {metadata_file}")
    return metadata_file


def main():
    """Main download function."""
    
    parser = argparse.ArgumentParser(description="Download GSE235063 dataset")
    parser.add_argument("--study", default="GSE235063", help="GEO study ID")
    parser.add_argument("--output", default=RAW_DATA_DIR, help="Output directory")
    parser.add_argument("--validate-only", action="store_true", help="Only validate existing files")
    parser.add_argument("--skip-download", action="store_true", help="Skip download if files exist")
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging("download_data.log")
    
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Check if data already exists
    if args.validate_only:
        if validate_downloads(output_dir):
            logger.info("Data validation passed")
            sys.exit(0)
        else:
            logger.error("Data validation failed")
            sys.exit(1)
    
    # Skip download if files exist and skip flag is set
    if args.skip_download and validate_downloads(output_dir):
        logger.info("Data already exists, skipping download")
        samples = extract_sample_data(output_dir, output_dir)
        create_sample_metadata(samples, output_dir)
        return
    
    # Download study
    logger.info(f"Starting download of {args.study}")
    downloaded_files = download_geo_study(args.study, output_dir)
    
    # Extract individual samples
    logger.info("Extracting sample data...")
    samples = extract_sample_data(output_dir, output_dir)
    
    # Create metadata
    metadata_file = create_sample_metadata(samples, output_dir)
    
    # Validate downloads
    if validate_downloads(output_dir):
        logger.info("Download completed successfully")
        logger.info(f"Downloaded {len(downloaded_files)} files")
        logger.info(f"Extracted {len(samples)} samples")
        logger.info(f"Metadata saved to: {metadata_file}")
    else:
        logger.error("Download validation failed")
        sys.exit(1)


if __name__ == "__main__":
    main()