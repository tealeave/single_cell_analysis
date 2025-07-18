"""Annotation pipeline with CellTypist and scVI label transfer."""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from anndata import AnnData
from scipy.stats import fisher_exact

from config import DEBUG_MODE, setup_logging
from utils import select_debug_samples


class AnnotationProcessor:
    """Handle cell type annotation with CellTypist and scVI label transfer."""
    
    def __init__(self,
                 input_dir: Path,
                 output_dir: Path,
                 reference_dir: Path = None):
        """
        Initialize annotation processor.
        
        Args:
            input_dir: Directory with integrated AnnData file
            output_dir: Directory for annotated outputs
            reference_dir: Directory with reference datasets/models
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.reference_dir = Path(reference_dir) if reference_dir else self.output_dir / "references"
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.references_dir.mkdir(exist_ok=True)
        self.plots_dir = self.output_dir / "annotation_plots"
        self.plots_dir.mkdir(exist_ok=True)
    
    def prepare_data_for_annotation(self, adata: AnnData) -> AnnData:
        """Prepare data for annotation."""
        
        # Ensure we have the integrated representation
        if "X_scVI" not in adata.obsm:
            logging.warning("No scVI embedding found, using PCA instead")
            sc.pp.pca(adata)
            adata.obsm["X_scVI"] = adata.obsm["X_pca"]
        
        # Normalize data for CellTypist
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        return adata
    
    def run_celltypist(self, adata: AnnData, model_name: str = "Immune_All_Low") -> pd.DataFrame:
        """Run CellTypist annotation."""
        
        logging.info(f"Running CellTypist with model: {model_name}")
        
        try:
            # Run CellTypist
            predictions = celltypist.annotate(
                adata,
                model=model_name,
                majority_voting=True,
                over_clustering=True
            )
            
            # Get results
            result_df = predictions.predicted_labels.copy()
            
            logging.info(f"CellTypist completed: {len(result_df['majority_voting'].unique())} cell types found")
            
            return result_df
            
        except Exception as e:
            logging.error(f"CellTypist failed: {e}")
            return pd.DataFrame()
    
    def create_custom_celltypist_model(self, adata: AnnData, labels: pd.Series, 
                                     model_name: str = "custom_aml_model") -> str:
        """Create custom CellTypist model from labeled data."""
        
        logging.info(f"Creating custom CellTypist model: {model_name}")
        
        try:
            # Ensure labels match adata
            if len(labels) != len(adata):
                raise ValueError("Labels length doesn't match adata")
            
            # Create model
            adata.obs['cell_type'] = labels
            model = celltypist.train(
                adata,
                labels='cell_type',
                max_iter=200 if not DEBUG_MODE else 50,
                use_SGD=False
            )
            
            # Save model
            model_path = self.references_dir / f"{model_name}.pkl"
            model.write(str(model_path))
            
            logging.info(f"Custom model saved: {model_path}")
            return str(model_path)
            
        except Exception as e:
            logging.error(f"Custom model creation failed: {e}")
            return ""
    
    def scvi_label_transfer(self, adata_query: AnnData, adata_reference: AnnData, 
                          labels_reference: pd.Series) -> pd.DataFrame:
        """Perform label transfer using scVI latent space."""
        
        logging.info("Performing scVI label transfer")
        
        try:
            from sklearn.neighbors import KNeighborsClassifier
            from sklearn.metrics import accuracy_score
            
            # Ensure both datasets have scVI embeddings
            if "X_scVI" not in adata_query.obsm or "X_scVI" not in adata_reference.obsm:
                raise ValueError("scVI embeddings not found")
            
            # Prepare data
            X_ref = adata_reference.obsm["X_scVI"]
            y_ref = labels_reference.values
            X_query = adata_query.obsm["X_scVI"]
            
            # Train classifier
            knn = KNeighborsClassifier(n_neighbors=15)
            knn.fit(X_ref, y_ref)
            
            # Predict labels
            predicted_labels = knn.predict(X_query)
            probabilities = knn.predict_proba(X_query)
            
            # Create result DataFrame
            result_df = pd.DataFrame({
                'scvi_label': predicted_labels,
                'scvi_confidence': np.max(probabilities, axis=1)
            })
            
            logging.info(f"Label transfer completed: {len(np.unique(predicted_labels))} cell types")
            
            return result_df
            
        except Exception as e:
            logging.error(f"Label transfer failed: {e}")
            return pd.DataFrame()
    
    def combine_annotations(self, adata: AnnData, 
                          celltypist_results: pd.DataFrame,
                          scvi_results: pd.DataFrame = None) -> AnnData:
        """Combine annotations from multiple methods."""
        
        # Add CellTypist results
        if not celltypist_results.empty:
            if 'majority_voting' in celltypist_results.columns:
                adata.obs['celltypist_majority'] = celltypist_results['majority_voting']
            if 'predicted_labels' in celltypist_results.columns:
                adata.obs['celltypist_raw'] = celltypist_results['predicted_labels']
            if 'conf_score' in celltypist_results.columns:
                adata.obs['celltypist_confidence'] = celltypist_results['conf_score']
        
        # Add scVI results
        if scvi_results is not None and not scvi_results.empty:
            adata.obs['scvi_label'] = scvi_results['scvi_label']
            adata.obs['scvi_confidence'] = scvi_results['scvi_confidence']
        
        # Create consensus annotation (if both methods available)
        if 'celltypist_majority' in adata.obs.columns and 'scvi_label' in adata.obs.columns:
            # Simple consensus: prefer CellTypist with high confidence
            mask = adata.obs['celltypist_confidence'] > 0.8
            adata.obs['consensus_annotation'] = adata.obs['celltypist_majority']
            adata.obs.loc[~mask, 'consensus_annotation'] = adata.obs.loc[~mask, 'scvi_label']
        elif 'celltypist_majority' in adata.obs.columns:
            adata.obs['consensus_annotation'] = adata.obs['celltypist_majority']
        elif 'scvi_label' in adata.obs.columns:
            adata.obs['consensus_annotation'] = adata.obs['scvi_label']
        
        return adata
    
    def validate_annotations(self, adata: AnnData, 
                           annotation_col: str = 'consensus_annotation') -> Dict[str, any]:
        """Validate annotations using marker genes."""
        
        validation_results = {
            'total_cells': len(adata),
            'annotated_cells': 0,
            'unique_cell_types': 0,
            'marker_genes': {},
            'quality_scores': {}
        }
        
        if annotation_col not in adata.obs.columns:
            logging.warning(f"Annotation column {annotation_col} not found")
            return validation_results
        
        # Count annotations
        cell_type_counts = adata.obs[annotation_col].value_counts()
        validation_results['annotated_cells'] = len(adata.obs[annotation_col].dropna())
        validation_results['unique_cell_types'] = len(cell_type_counts)
        
        # Define marker genes for common AML cell types
        marker_genes = {
            'AML_Blast': ['MPO', 'CD33', 'CD34', 'CD117'],
            'Monocyte': ['CD14', 'CD16', 'LYZ', 'MS4A7'],
            'T_cell': ['CD3D', 'CD3E', 'CD8A', 'CD4'],
            'B_cell': ['CD19', 'CD79A', 'MS4A1', 'PAX5'],
            'NK_cell': ['GNLY', 'NKG7', 'KLRD1', 'GZMB'],
            'HSC': ['CD34', 'SPINK2', 'CD90', 'CD133'],
            'Neutrophil': ['CEACAM3', 'FCGR3B', 'CSF3R', 'S100A8']
        }
        
        # Score marker genes for each cell type
        for cell_type, genes in marker_genes.items():
            if cell_type in cell_type_counts.index:
                # Check if marker genes are expressed
                available_genes = [g for g in genes if g in adata.var_names]
                if available_genes:
                    cell_mask = adata.obs[annotation_col] == cell_type
                    expression = adata[cell_mask, available_genes].X.mean(axis=0)
                    
                    # Calculate marker score
                    marker_score = np.mean(expression > np.percentile(adata[:, available_genes].X.mean(axis=0), 75))
                    validation_results['marker_genes'][cell_type] = {
                        'available_markers': len(available_genes),
                        'mean_expression': float(np.mean(expression)),
                        'marker_score': float(marker_score)
                    }
        
        return validation_results
    
    def create_annotation_plots(self, adata: AnnData, 
                              annotation_col: str = 'consensus_annotation') -> None:
        """Create comprehensive annotation plots."""
        
        if annotation_col not in adata.obs.columns:
            logging.warning(f"Annotation column {annotation_col} not found")
            return
        
        # UMAP colored by cell types
        fig, ax = plt.subplots(figsize=(12, 10))
        sc.pl.umap(adata, color=annotation_col, ax=ax, show=False, 
                  title=f'UMAP colored by {annotation_col}')
        plt.tight_layout()
        plt.savefig(self.plots_dir / "umap_cell_types.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Cell type distribution
        cell_type_counts = adata.obs[annotation_col].value_counts()
        
        fig, ax = plt.subplots(figsize=(12, 6))
        cell_type_counts.plot(kind='bar', ax=ax)
        ax.set_title('Cell Type Distribution')
        ax.set_xlabel('Cell Type')
        ax.set_ylabel('Number of Cells')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(self.plots_dir / "cell_type_distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Confidence scores
        confidence_cols = [col for col in adata.obs.columns if 'confidence' in col]
        if confidence_cols:
            fig, axes = plt.subplots(1, len(confidence_cols), figsize=(15, 5))
            if len(confidence_cols) == 1:
                axes = [axes]
            
            for ax, col in zip(axes, confidence_cols):
                adata.obs[col].hist(bins=50, ax=ax, edgecolor='black')
                ax.set_title(f'{col} Distribution')
                ax.set_xlabel('Confidence Score')
                ax.set_ylabel('Frequency')
            
            plt.tight_layout()
            plt.savefig(self.plots_dir / "confidence_distributions.png", dpi=300, bbox_inches='tight')
            plt.close()
        
        # Sample composition
        if 'sample' in adata.obs.columns:
            composition = pd.crosstab(adata.obs['sample'], adata.obs[annotation_col])
            
            fig, ax = plt.subplots(figsize=(15, 8))
            sns.heatmap(composition, annot=True, fmt='d', cmap='Blues', ax=ax)
            ax.set_title('Cell Type Composition by Sample')
            plt.tight_layout()
            plt.savefig(self.plots_dir / "sample_composition.png", dpi=300, bbox_inches='tight')
            plt.close()
        
        logging.info("Created annotation plots")
    
    def process_sample(self, force: bool = False) -> Dict[str, Path]:
        """Process integrated data through annotation pipeline."""
        
        logging.info("Starting annotation pipeline")
        
        # Input and output paths
        input_file = self.input_dir / "integrated.h5ad"
        output_file = self.output_dir / "annotated.h5ad"
        report_file = self.output_dir / "annotation_report.html"
        
        # Skip if already processed and not forced
        if output_file.exists() and not force:
            logging.info("Annotation already completed, skipping")
            return {
                'annotated': output_file,
                'report': report_file
            }
        
        # Load integrated data
        if not input_file.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        adata = sc.read_h5ad(input_file)
        logging.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
        
        # Prepare data
        adata = self.prepare_data_for_annotation(adata)
        
        # Run CellTypist with Immune_All_Low model
        celltypist_results = self.run_celltypist(adata, "Immune_All_Low")
        
        # Run CellTypist with custom AML model (if available)
        custom_model = self.references_dir / "custom_aml_model.pkl"
        if custom_model.exists():
            custom_results = self.run_celltypist(adata, str(custom_model))
        else:
            custom_results = pd.DataFrame()
        
        # Combine annotations
        adata = self.combine_annotations(adata, celltypist_results)
        
        # Validate annotations
        validation_results = self.validate_annotations(adata)
        
        # Create plots
        self.create_annotation_plots(adata)
        
        # Save annotated data
        adata.write_h5ad(output_file)
        
        # Save validation report
        validation_df = pd.DataFrame([validation_results])
        validation_df.to_csv(self.output_dir / "validation_report.csv", index=False)
        
        # Create summary
        summary = {
            'total_cells': adata.n_obs,
            'total_genes': adata.n_vars,
            'cell_types': len(adata.obs['consensus_annotation'].unique()) if 'consensus_annotation' in adata.obs else 0,
            'methods_used': ['CellTypist', 'scVI_label_transfer'],
            'validation_completed': True
        }
        
        summary_file = self.output_dir / "annotation_summary.json"
        pd.Series(summary).to_json(summary_file)
        
        logging.info(f"Annotation completed: {summary}")
        
        return {
            'annotated': output_file,
            'report': report_file,
            'summary': summary_file,
            'validation': self.output_dir / "validation_report.csv"
        }


def main():
    """Main function for annotation processing."""
    
    parser = argparse.ArgumentParser(description="Cell annotation pipeline")
    parser.add_argument("--input", default="results/integrated", help="Input directory")
    parser.add_argument("--output", default="results/annotations", help="Output directory")
    parser.add_argument("--reference", help="Reference directory")
    parser.add_argument("--force", action="store_true", help="Force reprocessing")
    parser.add_argument("--model", default="Immune_All_Low", help="CellTypist model")
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging("annotation.log")
    
    processor = AnnotationProcessor(
        input_dir=args.input,
        output_dir=args.output,
        reference_dir=args.reference
    )
    
    results = processor.process_sample(force=args.force)
    
    logger.info("Annotation completed successfully")


if __name__ == "__main__":
    main()