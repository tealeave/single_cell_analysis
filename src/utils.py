"""Utility functions for single-cell analysis pipeline."""

import logging
import random
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats


def setup_logging(log_file: str, level: int = logging.INFO) -> logging.Logger:
    """Set up logging configuration."""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)


def select_debug_samples(
    sample_files: List[Path], 
    n_samples: int = 3, 
    random_seed: int = 42
) -> List[Path]:
    """
    Select a representative subset of samples for debug mode.
    
    Args:
        sample_files: List of sample file paths
        n_samples: Number of samples to select
        random_seed: Random seed for reproducibility
    
    Returns:
        List of selected sample paths
    """
    if len(sample_files) <= n_samples:
        return sample_files
    
    random.seed(random_seed)
    
    # Select samples evenly distributed across the dataset
    indices = np.linspace(0, len(sample_files) - 1, n_samples, dtype=int)
    selected_samples = [sample_files[i] for i in indices]
    
    return selected_samples


def mad_outlier_detection(
    data: np.ndarray, 
    threshold: float = 3.0
) -> np.ndarray:
    """
    Detect outliers using Median Absolute Deviation (MAD).
    
    Args:
        data: Input data array
        threshold: MAD threshold for outlier detection
    
    Returns:
        Boolean array indicating outliers
    """
    median = np.median(data)
    mad = stats.median_abs_deviation(data)
    
    if mad == 0:
        return np.zeros_like(data, dtype=bool)
    
    modified_z_scores = 0.6745 * (data - median) / mad
    return np.abs(modified_z_scores) > threshold


def create_sample_manifest(
    data_dir: Path, 
    pattern: str = "*.h5ad"
) -> pd.DataFrame:
    """
    Create a manifest of samples with basic metadata.
    
    Args:
        data_dir: Directory containing sample files
        pattern: File pattern to match
    
    Returns:
        DataFrame with sample metadata
    """
    sample_files = list(data_dir.glob(pattern))
    
    manifest = []
    for file_path in sample_files:
        try:
            # Load minimal data to get basic info
            adata = sc.read_h5ad(file_path)
            manifest.append({
                'sample_id': file_path.stem,
                'file_path': str(file_path),
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'total_counts': adata.X.sum(),
            })
        except Exception as e:
            logging.warning(f"Could not process {file_path}: {e}")
    
    return pd.DataFrame(manifest)


def validate_sample_files(sample_files: List[Path]) -> Tuple[List[Path], List[str]]:
    """
    Validate sample files and return valid files with error messages.
    
    Args:
        sample_files: List of file paths to validate
    
    Returns:
        Tuple of (valid_files, error_messages)
    """
    valid_files = []
    errors = []
    
    for file_path in sample_files:
        if not file_path.exists():
            errors.append(f"File not found: {file_path}")
            continue
        
        try:
            # Quick validation - just check if file can be read
            sc.read_h5ad(file_path)
            valid_files.append(file_path)
        except Exception as e:
            errors.append(f"Invalid file {file_path}: {e}")
    
    return valid_files, errors


def get_resource_usage():
    """Get current resource usage for debugging."""
    import psutil
    
    process = psutil.Process()
    memory_info = process.memory_info()
    
    return {
        'memory_mb': memory_info.rss / 1024 / 1024,
        'cpu_percent': process.cpu_percent(),
        'num_threads': process.num_threads()
    }