"""CellBender ambient RNA removal module for single-cell analysis."""

import argparse
import logging
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import AnnData

from config import (
    CELLBENDER_EXPECTED_CELLS,
    CELLBENDER_FPR,
    CELLBENDER_TOTAL_DROPLETS,
    DEBUG_MODE,
    setup_logging,
)
from utils import get_resource_usage


class CellBenderProcessor:
    """Handle CellBender ambient RNA removal."""
    
    def __init__(self, 
                 input_dir: Path,
                 output_dir: Path,
                 expected_cells: int = CELLBENDER_EXPECTED_CELLS,
                 total_droplets: int = CELLBENDER_TOTAL_DROPLETS,
                 fpr: float = CELLBENDER_FPR):
        """
        Initialize CellBender processor.
        
        Args:
            input_dir: Directory with raw 10x data
            output_dir: Directory for CellBender outputs
            expected_cells: Expected number of cells
            total_droplets: Total droplets to include
            fpr: False positive rate
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.expected_cells = expected_cells
        self.total_droplets = total_droplets
        self.fpr = fpr
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def check_cellbender_installation(self) -> bool:
        """Check if CellBender is available."""
        try:
            subprocess.run(["cellbender", "--help"], 
                         capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            logging.error("CellBender not found. Install with: pip install cellbender")
            return False
    
    def prepare_input_files(self, sample_id: str) -> Path:
        """Prepare input directory structure for CellBender."""
        sample_dir = self.output_dir / sample_id
        sample_dir.mkdir(exist_ok=True)
        
        # Check if 10x files exist
        matrix_file = self.input_dir / sample_id / "matrix.mtx.gz"
        barcodes_file = self.input_dir / sample_id / "barcodes.tsv.gz"
        features_file = self.input_dir / sample_id / "features.tsv.gz"
        
        if not all([matrix_file.exists(), barcodes_file.exists(), 
                   features_file.exists()]):
            raise FileNotFoundError(
                f"Missing 10x files for {sample_id}: {matrix_file}, "
                f"{barcodes_file}, {features_file}"
            )
        
        return sample_dir
    
    def run_cellbender(self, sample_id: str, force: bool = False) -> Dict[str, Path]:
        """
        Run CellBender remove-background for a single sample.
        
        Args:
            sample_id: Sample identifier
            force: Force re-run even if output exists
            
        Returns:
            Dictionary with output file paths
        """
        sample_dir = self.prepare_input_files(sample_id)
        
        # Output files
        output_prefix = sample_dir / f"{sample_id}_cellbender"
        output_h5 = f"{output_prefix}.h5"
        output_filtered = f"{output_prefix}_filtered.h5"
        metrics_file = f"{output_prefix}_metrics.csv"
        
        # Skip if already processed and not forced
        if Path(output_h5).exists() and not force:
            logging.info(f"CellBender already run for {sample_id}, skipping")
            return {
                'raw': Path(output_h5),
                'filtered': Path(output_filtered),
                'metrics': Path(metrics_file)
            }
        
        # Build CellBender command
        cmd = [
            "cellbender", "remove-background",
            "--input", str(self.input_dir / sample_id),
            "--output", str(output_h5),
            "--expected-cells", str(self.expected_cells),
            "--total-droplets-included", str(self.total_droplets),
            "--fpr", str(self.fpr),
            "--epochs", "150" if not DEBUG_MODE else "50",
            "--cuda", "False",  # Force CPU mode
            "--low-count-threshold", "15"
        ]
        
        if DEBUG_MODE:
            cmd.extend(["--epochs", "10"])
        
        logging.info(f"Running CellBender for {sample_id}")
        logging.info(f"Command: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logging.info(f"CellBender completed for {sample_id}")
            
            # Parse metrics from output
            metrics = self.parse_cellbender_output(result.stdout, result.stderr)
            pd.DataFrame([metrics]).to_csv(metrics_file, index=False)
            
            return {
                'raw': Path(output_h5),
                'filtered': Path(output_filtered),
                'metrics': Path(metrics_file)
            }
            
        except subprocess.CalledProcessError as e:
            logging.error(f"CellBender failed for {sample_id}: {e}")
            logging.error(f"STDOUT: {e.stdout}")
            logging.error(f"STDERR: {e.stderr}")
            raise
    
    def parse_cellbender_output(self, stdout: str, stderr: str) -> Dict[str, float]:
        """Parse CellBender output metrics."""
        metrics = {
            'fraction_counts_removed': None,
            'estimated_cells': None,
            'ambient_RNA_fraction': None
        }
        
        lines = (stdout + stderr).split('\n')
        for line in lines:
            if "Fraction of counts due to ambient RNA" in line:
                try:
                    metrics['fraction_counts_removed'] = float(
                        line.split(':')[-1].strip().replace('%', '')) / 100
                except ValueError:
                    pass
            elif "Estimated number of cells" in line:
                try:
                    metrics['estimated_cells'] = int(line.split(':')[-1].strip())
                except ValueError:
                    pass
            elif "Ambient RNA fraction" in line:
                try:
                    metrics['ambient_RNA_fraction'] = float(
                        line.split(':')[-1].strip())
                except ValueError:
                    pass
        
        return metrics
    
    def convert_to_anndata(self, cellbender_h5: Path, output_h5ad: Path) -> AnnData:
        """Convert CellBender output to AnnData format."""
        
        logging.info(f"Converting {cellbender_h5} to AnnData")
        
        # Read CellBender output
        adata = sc.read_10x_h5(str(cellbender_h5))
        
        # Add metadata
        adata.var_names_make_unique()
        adata.var['gene_ids'] = adata.var.index
        adata.var['feature_types'] = 'Gene Expression'
        
        # Add cellbender info
        adata.uns['cellbender'] = {
            'expected_cells': self.expected_cells,
            'total_droplets': self.total_droplets,
            'fpr': self.fpr,
            'processed_at': pd.Timestamp.now().isoformat()
        }
        
        # Save as AnnData
        adata.write_h5ad(output_h5ad)
        logging.info(f"Saved AnnData to {output_h5ad}")
        
        return adata
    
    def create_qc_plots(self, adata: AnnData, sample_id: str, output_dir: Path) -> None:
        """Create QC plots for CellBender results."""
        
        plots_dir = output_dir / "qc_plots"
        plots_dir.mkdir(exist_ok=True)
        
        # Basic QC metrics
        sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
        
        # Plot distributions
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        # Total counts per cell
        axes[0, 0].hist(adata.obs['n_counts'], bins=50, edgecolor='black')
        axes[0, 0].set_xlabel('Total counts per cell')
        axes[0, 0].set_ylabel('Frequency')
        axes[0, 0].set_title(f'{sample_id}: Counts Distribution')
        
        # Number of genes per cell
        axes[0, 1].hist(adata.obs['n_genes_by_counts'], bins=50, edgecolor='black')
        axes[0, 1].set_xlabel('Number of genes per cell')
        axes[0, 1].set_ylabel('Frequency')
        axes[0, 1].set_title(f'{sample_id}: Genes Distribution')
        
        # Mitochondrial gene percentage
        mito_genes = adata.var_names.str.startswith('MT-')
        if mito_genes.any():
            adata.obs['percent_mito'] = adata[:, mito_genes].X.sum(axis=1) / adata.X.sum(axis=1) * 100
            axes[1, 0].hist(adata.obs['percent_mito'], bins=50, edgecolor='black')
            axes[1, 0].set_xlabel('Mitochondrial gene %')
            axes[1, 0].set_ylabel('Frequency')
            axes[1, 0].set_title(f'{sample_id}: Mitochondrial %')
        
        # Scatter plot: n_genes vs n_counts
        axes[1, 1].scatter(adata.obs['n_counts'], adata.obs['n_genes_by_counts'], alpha=0.5)
        axes[1, 1].set_xlabel('Total counts')
        axes[1, 1].set_ylabel('Number of genes')
        axes[1, 1].set_title(f'{sample_id}: Genes vs Counts')
        
        plt.tight_layout()
        plt.savefig(plots_dir / f"{sample_id}_qc_plots.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"Created QC plots for {sample_id}")
    
    def process_sample(self, sample_id: str, force: bool = False) -> Dict[str, Path]:
        """Process a single sample through CellBender pipeline."""
        
        logging.info(f"Processing sample: {sample_id}")
        
        # Run CellBender
        outputs = self.run_cellbender(sample_id, force)
        
        # Convert to AnnData
        output_h5ad = self.output_dir / f"{sample_id}_denoised.h5ad"
        adata = self.convert_to_anndata(outputs['raw'], output_h5ad)
        
        # Create QC plots
        self.create_qc_plots(adata, sample_id, self.output_dir)
        
        # Log resource usage
        usage = get_resource_usage()
        logging.info(f"Resource usage: {usage}")
        
        return {
            'denoised': output_h5ad,
            'metrics': outputs['metrics'],
            'qc_plots': self.output_dir / "qc_plots" / f"{sample_id}_qc_plots.png"
        }
    
    def process_all_samples(self, sample_list: list = None, force: bool = False) -> Dict[str, Dict[str, Path]]:
        """Process all samples in the dataset."""
        
        if sample_list is None:
            # Discover samples from input directory
            sample_list = [d.name for d in self.input_dir.iterdir() if d.is_dir()]
        
        if DEBUG_MODE and len(sample_list) > 3:
            sample_list = sample_list[:3]
            logging.info(f"Debug mode: processing only {len(sample_list)} samples")
        
        results = {}
        
        for sample_id in sample_list:
            try:
                results[sample_id] = self.process_sample(sample_id, force)
                logging.info(f"Successfully processed {sample_id}")
            except Exception as e:
                logging.error(f"Failed to process {sample_id}: {e}")
                continue
        
        return results


def main():
    """Main function for CellBender processing."""
    
    parser = argparse.ArgumentParser(description="CellBender ambient RNA removal")
    parser.add_argument("--input", default="data/raw", help="Input directory")
    parser.add_argument("--output", default="data/processed", help="Output directory")
    parser.add_argument("--sample", help="Process specific sample")
    parser.add_argument("--force", action="store_true", help="Force reprocessing")
    parser.add_argument("--expected-cells", type=int, default=CELLBENDER_EXPECTED_CELLS)
    parser.add_argument("--total-droplets", type=int, default=CELLBENDER_TOTAL_DROPLETS)
    parser.add_argument("--fpr", type=float, default=CELLBENDER_FPR)
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging("preprocessing.log")
    
    processor = CellBenderProcessor(
        input_dir=args.input,
        output_dir=args.output,
        expected_cells=args.expected_cells,
        total_droplets=args.total_droplets,
        fpr=args.fpr
    )
    
    # Check CellBender installation
    if not processor.check_cellbender_installation():
        sys.exit(1)
    
    if args.sample:
        # Process single sample
        results = processor.process_sample(args.sample, args.force)
    else:
        # Process all samples
        results = processor.process_all_samples(force=args.force)
    
    logger.info(f"Processed {len(results)} samples successfully")
    
    # Create summary report
    summary_file = Path(args.output) / "cellbender_summary.csv"
    summary_data = []
    
    for sample_id, outputs in results.items():
        if outputs['metrics'].exists():
            metrics_df = pd.read_csv(outputs['metrics'])
            summary_data.append({
                'sample_id': sample_id,
                **metrics_df.iloc[0].to_dict()
            })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(summary_file, index=False)
        logger.info(f"Created summary: {summary_file}")


if __name__ == "__main__":
    main()