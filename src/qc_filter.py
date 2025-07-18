"""Quality control and filtering module with MAD outlier detection."""

import argparse
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import doubletdetection as dd
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import AnnData
from scipy import stats

from config import (
    DEBUG_MODE,
    QC_MAX_GENES,
    QC_MAX_MT_PERCENT,
    QC_MIN_CELLS,
    QC_MIN_GENES,
    QC_MAD_THRESHOLD,
    setup_logging,
)
from utils import mad_outlier_detection, select_debug_samples


class QCProcessor:
    """Handle quality control and filtering."""
    
    def __init__(self,
                 input_dir: Path,
                 output_dir: Path,
                 min_genes: int = QC_MIN_GENES,
                 max_genes: int = QC_MAX_GENES,
                 max_mt_percent: float = QC_MAX_MT_PERCENT,
                 min_cells: int = QC_MIN_CELLS,
                 mad_threshold: float = QC_MAD_THRESHOLD):
        """
        Initialize QC processor.
        
        Args:
            input_dir: Directory with CellBender outputs
            output_dir: Directory for filtered outputs
            min_genes: Minimum genes per cell
            max_genes: Maximum genes per cell
            max_mt_percent: Maximum mitochondrial gene percentage
            min_cells: Minimum cells per gene
            mad_threshold: MAD threshold for outlier detection
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.min_genes = min_genes
        self.max_genes = max_genes
        self.max_mt_percent = max_mt_percent
        self.min_cells = min_cells
        self.mad_threshold = mad_threshold
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.qc_reports_dir = self.output_dir / "qc_reports"
        self.qc_reports_dir.mkdir(exist_ok=True)
        
    def calculate_qc_metrics(self, adata: AnnData) -> AnnData:
        """Calculate QC metrics for the dataset."""
        
        # Basic QC metrics
        sc.pp.calculate_qc_metrics(
            adata, 
            percent_top=None, 
            log1p=False, 
            inplace=True
        )
        
        # Calculate mitochondrial gene percentage
        mito_genes = adata.var_names.str.startswith('MT-')
        if mito_genes.any():
            adata.obs['percent_mt'] = adata[:, mito_genes].X.sum(axis=1) / adata.X.sum(axis=1) * 100
        else:
            # Try lowercase MT
            mito_genes = adata.var_names.str.startswith('mt-')
            if mito_genes.any():
                adata.obs['percent_mt'] = adata[:, mito_genes].X.sum(axis=1) / adata.X.sum(axis=1) * 100
            else:
                adata.obs['percent_mt'] = 0
        
        # Calculate ribosomal gene percentage
        ribo_genes = adata.var_names.str.startswith(('RPS', 'RPL'))
        if ribo_genes.any():
            adata.obs['percent_ribo'] = adata[:, ribo_genes].X.sum(axis=1) / adata.X.sum(axis=1) * 100
        else:
            adata.obs['percent_ribo'] = 0
        
        return adata
    
    def detect_outliers_mad(self, adata: AnnData) -> pd.Series:
        """Detect outliers using MAD-based approach."""
        
        metrics_to_check = [
            'n_counts',
            'n_genes_by_counts', 
            'percent_mt'
        ]
        
        outlier_mask = pd.Series([False] * len(adata), index=adata.obs.index)
        
        for metric in metrics_to_check:
            if metric in adata.obs.columns:
                data = adata.obs[metric].values
                outliers = mad_outlier_detection(data, self.mad_threshold)
                outlier_mask |= outliers
        
        adata.obs['is_outlier'] = outlier_mask
        return outlier_mask
    
    def detect_doublets(self, adata: AnnData) -> pd.Series:
        """Detect doublets using DoubletDetection."""
        
        logging.info("Running doublet detection...")
        
        # Prepare data for doublet detection
        counts = adata.X.copy()
        
        # Doublet detection
        clf = dd.BoostClassifier(
            n_iters=25 if not DEBUG_MODE else 10,
            clustering_algorithm='phenograph',
            standard_scaling=True
        )
        
        try:
            doublets = clf.fit(counts).predict(p_thresh=1e-16, voter_thresh=0.5)
            doublet_score = clf.doublet_score()
            
            adata.obs['doublet'] = doublets.astype(bool)
            adata.obs['doublet_score'] = doublet_score
            
            logging.info(f"Detected {sum(doublets)} doublets out of {len(adata)} cells")
            
        except Exception as e:
            logging.warning(f"Doublet detection failed: {e}")
            adata.obs['doublet'] = False
            adata.obs['doublet_score'] = 0.0
        
        return adata.obs['doublet']
    
    def apply_filters(self, adata: AnnData) -> Tuple[AnnData, pd.DataFrame]:
        """Apply QC filters and return filtered data."""
        
        initial_cells = adata.n_obs
        initial_genes = adata.n_vars
        
        # Basic filtering
        sc.pp.filter_cells(adata, min_genes=self.min_genes)
        sc.pp.filter_cells(adata, max_genes=self.max_genes)
        sc.pp.filter_genes(adata, min_cells=self.min_cells)
        
        # Mitochondrial percentage filter
        mito_filter = adata.obs['percent_mt'] <= self.max_mt_percent
        adata = adata[mito_filter, :]
        
        # Remove outliers
        if 'is_outlier' in adata.obs.columns:
            outlier_filter = ~adata.obs['is_outlier']
            adata = adata[outlier_filter, :]
        
        # Remove doublets
        if 'doublet' in adata.obs.columns:
            doublet_filter = ~adata.obs['doublet']
            adata = adata[doublet_filter, :]
        
        # Create filtering report
        filtering_report = pd.DataFrame({
            'step': [
                'initial',
                'min_genes_filter',
                'max_genes_filter',
                'min_cells_filter',
                'mito_filter',
                'outlier_filter',
                'doublet_filter',
                'final'
            ],
            'n_cells': [
                initial_cells,
                None,  # Will be calculated
                None,
                None,
                None,
                None,
                None,
                adata.n_obs
            ],
            'n_genes': [
                initial_genes,
                None,
                None,
                adata.n_vars,
                adata.n_vars,
                adata.n_vars,
                adata.n_vars,
                adata.n_vars
            ]
        })
        
        return adata, filtering_report
    
    def create_qc_plots(self, adata: AnnData, sample_id: str) -> None:
        """Create comprehensive QC plots."""
        
        plots_dir = self.qc_reports_dir / sample_id
        plots_dir.mkdir(exist_ok=True)
        
        # Calculate metrics if not already done
        if 'n_counts' not in adata.obs.columns:
            adata = self.calculate_qc_metrics(adata)
        
        # Ridge plots for distributions
        fig, axes = plt.subplots(2, 3, figsize=(18, 12))
        
        metrics = ['n_counts', 'n_genes_by_counts', 'percent_mt']
        titles = ['Total Counts', 'Number of Genes', 'Mitochondrial %']
        
        for i, (metric, title) in enumerate(zip(metrics, titles)):
            if metric in adata.obs.columns:
                # Histogram
                axes[0, i].hist(adata.obs[metric], bins=50, edgecolor='black', alpha=0.7)
                axes[0, i].set_title(f'{sample_id}: {title}')
                axes[0, i].set_xlabel(title)
                axes[0, i].set_ylabel('Frequency')
                
                # Boxplot with outliers
                axes[1, i].boxplot(adata.obs[metric])
                axes[1, i].set_title(f'{sample_id}: {title} (Boxplot)')
                axes[1, i].set_ylabel(title)
        
        plt.tight_layout()
        plt.savefig(plots_dir / "qc_distributions.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Scatter plots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # n_genes vs n_counts
        axes[0, 0].scatter(adata.obs['n_counts'], adata.obs['n_genes_by_counts'], alpha=0.5)
        axes[0, 0].set_xlabel('Total counts')
        axes[0, 0].set_ylabel('Number of genes')
        axes[0, 0].set_title(f'{sample_id}: Genes vs Counts')
        
        # n_counts vs percent_mt
        axes[0, 1].scatter(adata.obs['n_counts'], adata.obs['percent_mt'], alpha=0.5)
        axes[0, 1].set_xlabel('Total counts')
        axes[0, 1].set_ylabel('Mitochondrial %')
        axes[0, 1].set_title(f'{sample_id}: Counts vs Mitochondrial %')
        
        # n_genes vs percent_mt
        axes[1, 0].scatter(adata.obs['n_genes_by_counts'], adata.obs['percent_mt'], alpha=0.5)
        axes[1, 0].set_xlabel('Number of genes')
        axes[1, 0].set_ylabel('Mitochondrial %')
        axes[1, 0].set_title(f'{sample_id}: Genes vs Mitochondrial %')
        
        # Doublet score distribution
        if 'doublet_score' in adata.obs.columns:
            axes[1, 1].hist(adata.obs['doublet_score'], bins=50, edgecolor='black')
            axes[1, 1].set_xlabel('Doublet score')
            axes[1, 1].set_ylabel('Frequency')
            axes[1, 1].set_title(f'{sample_id}: Doublet Score Distribution')
        
        plt.tight_layout()
        plt.savefig(plots_dir / "qc_scatter_plots.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"Created QC plots for {sample_id}")
    
    def process_sample(self, sample_id: str, force: bool = False) -> Dict[str, Path]:
        """Process a single sample through QC pipeline."""
        
        logging.info(f"Processing QC for sample: {sample_id}")
        
        # Input and output paths
        input_file = self.input_dir / f"{sample_id}_denoised.h5ad"
        output_file = self.output_dir / f"{sample_id}_filtered.h5ad"
        report_file = self.output_dir / "qc_reports" / f"{sample_id}_qc_report.html"
        
        # Skip if already processed and not forced
        if output_file.exists() and not force:
            logging.info(f"QC already run for {sample_id}, skipping")
            return {
                'filtered': output_file,
                'report': report_file
            }
        
        # Load data
        if not input_file.exists():
            logging.error(f"Input file not found: {input_file}")
            return {}
        
        adata = sc.read_h5ad(input_file)
        logging.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
        
        # Calculate QC metrics
        adata = self.calculate_qc_metrics(adata)
        
        # Detect outliers
        outlier_mask = self.detect_outliers_mad(adata)
        
        # Detect doublets
        doublet_mask = self.detect_doublets(adata)
        
        # Apply filters
        filtered_adata, filtering_report = self.apply_filters(adata)
        
        # Create QC plots
        self.create_qc_plots(adata, sample_id)
        
        # Save filtered data
        filtered_adata.write_h5ad(output_file)
        
        # Save filtering report
        report_file.parent.mkdir(exist_ok=True)
        filtering_report.to_csv(
            self.output_dir / "qc_reports" / f"{sample_id}_filtering_report.csv",
            index=False
        )
        
        # Create summary stats
        summary_stats = {
            'sample_id': sample_id,
            'initial_cells': adata.n_obs,
            'initial_genes': adata.n_vars,
            'final_cells': filtered_adata.n_obs,
            'final_genes': filtered_adata.n_vars,
            'cells_removed': adata.n_obs - filtered_adata.n_obs,
            'genes_removed': adata.n_vars - filtered_adata.n_vars,
            'outliers_detected': outlier_mask.sum() if 'is_outlier' in adata.obs else 0,
            'doublets_detected': doublet_mask.sum() if 'doublet' in adata.obs else 0,
            'percent_cells_kept': (filtered_adata.n_obs / adata.n_obs) * 100,
            'percent_genes_kept': (filtered_adata.n_vars / adata.n_vars) * 100
        }
        
        logging.info(f"QC completed for {sample_id}")
        logging.info(f"Kept {summary_stats['final_cells']} cells ({summary_stats['percent_cells_kept']:.1f}%)")
        logging.info(f"Kept {summary_stats['final_genes']} genes ({summary_stats['percent_genes_kept']:.1f}%)")
        
        return {
            'filtered': output_file,
            'report': report_file,
            'summary': summary_stats
        }
    
    def process_all_samples(self, sample_list: list = None, force: bool = False) -> Dict[str, Dict[str, Path]]:
        """Process all samples through QC pipeline."""
        
        if sample_list is None:
            # Discover samples from input directory
            sample_list = [f.stem.replace('_denoised', '') 
                          for f in self.input_dir.glob('*_denoised.h5ad')]
        
        if DEBUG_MODE and len(sample_list) > 3:
            sample_list = sample_list[:3]
            logging.info(f"Debug mode: processing only {len(sample_list)} samples")
        
        results = {}
        all_summaries = []
        
        for sample_id in sample_list:
            try:
                result = self.process_sample(sample_id, force)
                if result:
                    results[sample_id] = result
                    all_summaries.append(result['summary'])
                    logging.info(f"Successfully processed QC for {sample_id}")
            except Exception as e:
                logging.error(f"Failed to process QC for {sample_id}: {e}")
                continue
        
        # Create overall summary
        if all_summaries:
            summary_df = pd.DataFrame(all_summaries)
            summary_file = self.output_dir / "qc_summary.csv"
            summary_df.to_csv(summary_file, index=False)
            logging.info(f"Created QC summary: {summary_file}")
        
        return results


def main():
    """Main function for QC processing."""
    
    parser = argparse.ArgumentParser(description="Quality control and filtering")
    parser.add_argument("--input", default="data/processed", help="Input directory")
    parser.add_argument("--output", default="data/processed", help="Output directory")
    parser.add_argument("--sample", help="Process specific sample")
    parser.add_argument("--force", action="store_true", help="Force reprocessing")
    parser.add_argument("--min-genes", type=int, default=QC_MIN_GENES)
    parser.add_argument("--max-genes", type=int, default=QC_MAX_GENES)
    parser.add_argument("--max-mt-percent", type=float, default=QC_MAX_MT_PERCENT)
    parser.add_argument("--min-cells", type=int, default=QC_MIN_CELLS)
    parser.add_argument("--mad-threshold", type=float, default=QC_MAD_THRESHOLD)
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging("qc_filter.log")
    
    processor = QCProcessor(
        input_dir=args.input,
        output_dir=args.output,
        min_genes=args.min_genes,
        max_genes=args.max_genes,
        max_mt_percent=args.max_mt_percent,
        min_cells=args.min_cells,
        mad_threshold=args.mad_threshold
    )
    
    if args.sample:
        # Process single sample
        results = processor.process_sample(args.sample, args.force)
    else:
        # Process all samples
        results = processor.process_all_samples(force=args.force)
    
    logger.info(f"QC completed for {len(results)} samples")


if __name__ == "__main__":
    main()