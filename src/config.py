"""Configuration module for single-cell analysis pipeline."""

import logging
import os
from pathlib import Path


def setup_logging(log_file: str = "pipeline.log", level: int = logging.INFO) -> logging.Logger:
    """Set up logging configuration."""
    
    # Create logs directory if it doesn't exist
    LOGS_DIR.mkdir(exist_ok=True)
    
    # Configure logging
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(LOGS_DIR / log_file),
            logging.StreamHandler()
        ]
    )
    
    # Suppress verbose logging from some libraries
    logging.getLogger('scvi').setLevel(logging.WARNING)
    logging.getLogger('celltypist').setLevel(logging.WARNING)
    
    logger = logging.getLogger('pipeline')
    logger.info(f"Logging initialized: {log_file}")
    
    return logger

# Base paths
PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
LOGS_DIR = PROJECT_ROOT / "logs"

# Debug mode configuration
DEBUG_MODE = os.getenv("DEBUG_MODE", "false").lower() == "true"
DEBUG_SAMPLE_COUNT = 3
FULL_SAMPLE_COUNT = 75

# Resource configuration
if DEBUG_MODE:
    CELLBENDER_CPUS = 2
    CELLBENDER_MEMORY = "8GB"
    CELLBENDER_TIME = "2h"
    QC_CPUS = 2
    QC_MEMORY = "4GB"
    QC_TIME = "1h"
else:
    CELLBENDER_CPUS = 8
    CELLBENDER_MEMORY = "32GB"
    CELLBENDER_TIME = "8h"
    QC_CPUS = 4
    QC_MEMORY = "16GB"
    QC_TIME = "4h"

# Data paths
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
REFERENCES_DIR = DATA_DIR / "references"

# Results paths
QC_REPORTS_DIR = RESULTS_DIR / "qc_reports"
INTEGRATED_DIR = RESULTS_DIR / "integrated"
ANNOTATIONS_DIR = RESULTS_DIR / "annotations"

# CellBender settings
CELLBENDER_EXPECTED_CELLS = 5000
CELLBENDER_TOTAL_DROPLETS = 50000 if not DEBUG_MODE else 10000
CELLBENDER_FPR = 0.01

# QC thresholds
QC_MIN_GENES = 200
QC_MAX_GENES = 6000
QC_MAX_MT_PERCENT = 20
QC_MIN_CELLS = 3
QC_MAD_THRESHOLD = 3

# Integration settings
SCVI_LATENT_DIM = 30
SCVI_N_LAYERS = 2
SCVI_N_HIDDEN = 128
SCVI_MAX_EPOCHS = 400 if not DEBUG_MODE else 50
SCVI_BATCH_SIZE = 128
SCVI_LEARNING_RATE = 1e-3

# Debug sample selection criteria
DEBUG_SAMPLE_SELECTION = {
    "criteria": "representative",  # low, medium, high cell count
    "min_cells": 1000,
    "max_cells": 10000,
    "random_seed": 42
}