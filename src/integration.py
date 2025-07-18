"""scVI integration module with hyperparameter tuning using Ray Tune."""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
from anndata import AnnData
from ray import tune
from ray.tune import CLIReporter
from ray.tune.schedulers import ASHAScheduler
from ray.tune.search.hyperopt import HyperOptSearch

from config import (
    DEBUG_MODE,
    SCVI_BATCH_SIZE,
    SCVI_LATENT_DIM,
    SCVI_LEARNING_RATE,
    SCVI_MAX_EPOCHS,
    SCVI_N_HIDDEN,
    SCVI_N_LAYERS,
    setup_logging,
)


class IntegrationProcessor:
    """Handle scVI integration with hyperparameter tuning."""
    
    def __init__(self,
                 input_dir: Path,
                 output_dir: Path,
                 latent_dim: int = SCVI_LATENT_DIM,
                 max_epochs: int = SCVI_MAX_EPOCHS,
                 batch_size: int = SCVI_BATCH_SIZE,
                 learning_rate: float = SCVI_LEARNING_RATE):
        """
        Initialize integration processor.
        
        Args:
            input_dir: Directory with filtered AnnData files
            output_dir: Directory for integrated outputs
            latent_dim: scVI latent dimension
            max_epochs: Maximum training epochs
            batch_size: Training batch size
            learning_rate: Learning rate
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.latent_dim = latent_dim
        self.max_epochs = max_epochs
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.models_dir = self.output_dir / "models"
        self.models_dir.mkdir(exist_ok=True)
        
        # Configure scVI for CPU
        scvi.settings.seed = 42
        if not torch.cuda.is_available():
            scvi.settings.dl_num_workers = 0
    
    def load_and_concatenate_data(self, sample_list: List[str]) -> AnnData:
        """Load and concatenate all sample data."""
        
        adatas = []
        
        for sample_id in sample_list:
            input_file = self.input_dir / f"{sample_id}_filtered.h5ad"
            
            if not input_file.exists():
                logging.warning(f"File not found: {input_file}")
                continue
            
            adata = sc.read_h5ad(input_file)
            adata.obs['sample'] = sample_id
            adatas.append(adata)
        
        if not adatas:
            raise ValueError("No valid input files found")
        
        # Concatenate all samples
        concatenated = sc.concat(adatas, join='outer', label='sample')
        
        # Handle missing values
        concatenated.X = concatenated.X.fillna(0)
        concatenated.obs = concatenated.obs.fillna('')
        concatenated.var = concatenated.var.fillna('')
        
        logging.info(f"Concatenated {len(adatas)} samples: {concatenated.n_obs} cells, {concatenated.n_vars} genes")
        
        return concatenated
    
    def prepare_data_for_scvi(self, adata: AnnData) -> AnnData:
        """Prepare data for scVI training."""
        
        # Ensure we have raw counts
        if 'counts' not in adata.layers:
            adata.layers['counts'] = adata.X.copy()
        
        # Make variable names unique
        adata.var_names_make_unique()
        
        # Filter highly variable genes if requested
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=2000 if not DEBUG_MODE else 1000,
            subset=False,
            flavor="seurat_v3",
            layer="counts",
            batch_key="sample"
        )
        
        # Set up scVI
        scvi.model.SCVI.setup_anndata(
            adata,
            layer="counts",
            batch_key="sample"
        )
        
        return adata
    
    def train_scvi_model(self, adata: AnnData, hyperparams: Dict) -> scvi.model.SCVI:
        """Train scVI model with given hyperparameters."""
        
        model = scvi.model.SCVI(
            adata,
            n_layers=hyperparams.get('n_layers', 2),
            n_latent=hyperparams.get('n_latent', self.latent_dim),
            n_hidden=hyperparams.get('n_hidden', 128),
            dropout_rate=hyperparams.get('dropout_rate', 0.1),
            dispersion=hyperparams.get('dispersion', 'gene'),
            gene_likelihood=hyperparams.get('gene_likelihood', 'nb'),
            use_observed_lib_size=hyperparams.get('use_observed_lib_size', True)
        )
        
        early_stopping_kwargs = {
            "early_stopping": True,
            "early_stopping_monitor": "elbo_train",
            "early_stopping_patience": 15,
            "early_stopping_min_delta": 0.0
        }
        
        model.train(
            max_epochs=hyperparams.get('max_epochs', self.max_epochs),
            batch_size=hyperparams.get('batch_size', self.batch_size),
            lr=hyperparams.get('lr', self.learning_rate),
            use_gpu=torch.cuda.is_available(),
            **early_stopping_kwargs
        )
        
        return model
    
    def evaluate_model(self, model: scvi.model.SCVI, adata: AnnData) -> Dict[str, float]:
        """Evaluate scVI model performance."""
        
        metrics = {}
        
        try:
            # Get reconstruction loss
            reconstruction_loss = model.get_reconstruction_error()
            metrics['reconstruction_loss'] = float(reconstruction_loss)
            
            # Get latent representation
            latent = model.get_latent_representation()
            
            # Calculate variance explained
            total_variance = np.var(adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X)
            latent_variance = np.sum(np.var(latent, axis=0))
            metrics['variance_explained'] = float(latent_variance / total_variance)
            
            # Calculate batch mixing score (simplified)
            # This is a placeholder - in practice, you'd use more sophisticated metrics
            unique_batches = adata.obs['sample'].unique()
            if len(unique_batches) > 1:
                batch_entropy = self.calculate_batch_mixing_score(latent, adata.obs['sample'])
                metrics['batch_mixing_score'] = float(batch_entropy)
            
        except Exception as e:
            logging.warning(f"Error evaluating model: {e}")
            metrics['reconstruction_loss'] = float('inf')
            metrics['variance_explained'] = 0.0
            metrics['batch_mixing_score'] = 0.0
        
        return metrics
    
    def calculate_batch_mixing_score(self, latent: np.ndarray, batch_labels: pd.Series) -> float:
        """Calculate batch mixing score."""
        
        from sklearn.metrics import silhouette_score
        
        try:
            # Use negative silhouette score (lower is better for batch mixing)
            score = silhouette_score(latent, batch_labels)
            return max(0, 1 - abs(score))  # Normalize to [0,1]
        except:
            return 0.5  # Default score if calculation fails
    
    def hyperparameter_search(self, adata: AnnData, n_trials: int = 50) -> Dict:
        """Perform hyperparameter tuning for scVI."""
        
        if DEBUG_MODE:
            n_trials = 5
            logging.info(f"Debug mode: running {n_trials} trials")
        
        search_space = {
            'n_layers': tune.choice([1, 2, 3]),
            'n_latent': tune.choice([20, 30, 50]),
            'n_hidden': tune.choice([64, 128, 256]),
            'dropout_rate': tune.choice([0.1, 0.2, 0.3]),
            'lr': tune.loguniform(1e-4, 1e-2),
            'batch_size': tune.choice([64, 128, 256]),
            'max_epochs': tune.choice([200, 300, 400]) if not DEBUG_MODE else tune.choice([50, 100]),
            'dispersion': tune.choice(['gene', 'gene-batch']),
            'gene_likelihood': tune.choice(['nb', 'zinb'])
        }
        
        def trainable(config):
            """Training function for Ray Tune."""
            try:
                model = self.train_scvi_model(adata, config)
                metrics = self.evaluate_model(model, adata)
                
                # Tune minimizes the metric, so use negative reconstruction loss
                tune.report(
                    reconstruction_loss=metrics['reconstruction_loss'],
                    variance_explained=metrics['variance_explained'],
                    batch_mixing_score=metrics.get('batch_mixing_score', 0.0)
                )
                
                # Clean up model to save memory
                del model
                
            except Exception as e:
                logging.error(f"Training failed: {e}")
                tune.report(reconstruction_loss=float('inf'))
        
        # Set up search algorithm
        search_alg = HyperOptSearch(metric="reconstruction_loss", mode="min")
        
        # Set up scheduler
        scheduler = ASHAScheduler(
            max_t=400 if not DEBUG_MODE else 100,
            grace_period=10,
            reduction_factor=2
        )
        
        # Run hyperparameter search
        analysis = tune.run(
            trainable,
            config=search_space,
            num_samples=n_trials,
            search_alg=search_alg,
            scheduler=scheduler,
            resources_per_trial={"cpu": 4, "gpu": 0},
            progress_reporter=CLIReporter(
                metric_columns=["reconstruction_loss", "variance_explained"],
                max_report_frequency=30
            ),
            local_dir=str(self.output_dir / "ray_results"),
            name="scvi_hyperopt"
        )
        
        return analysis
    
    def train_best_model(self, adata: AnnData, analysis) -> scvi.model.SCVI:
        """Train final model with best hyperparameters."""
        
        best_config = analysis.get_best_config(metric="reconstruction_loss", mode="min")
        logging.info(f"Best hyperparameters: {best_config}")
        
        # Train final model
        model = self.train_scvi_model(adata, best_config)
        
        # Save model
        model_file = self.models_dir / "best_scvi_model"
        model.save(str(model_file))
        
        # Save best config
        config_file = self.models_dir / "best_config.json"
        pd.Series(best_config).to_json(config_file)
        
        return model
    
    def integrate_data(self, model: scvi.model.SCVI, adata: AnnData) -> AnnData:
        """Generate integrated data with embeddings."""
        
        # Get latent representation
        latent = model.get_latent_representation()
        adata.obsm["X_scVI"] = latent
        
        # Get normalized expression
        adata.layers["denoised"] = model.get_normalized_expression()
        
        # Run UMAP on latent space
        sc.pp.neighbors(adata, use_rep="X_scVI")
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=1.0)
        
        # Add integration metadata
        adata.uns['integration'] = {
            'method': 'scVI',
            'n_latent': latent.shape[1],
            'n_samples': len(adata.obs['sample'].unique()),
            'processed_at': pd.Timestamp.now().isoformat()
        }
        
        return adata
    
    def create_integration_plots(self, adata: AnnData) -> None:
        """Create integration visualization plots."""
        
        plots_dir = self.output_dir / "integration_plots"
        plots_dir.mkdir(exist_ok=True)
        
        # UMAP colored by batch
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        
        # UMAP by sample
        sc.pl.umap(adata, color='sample', ax=axes[0], show=False, title='UMAP by Sample')
        
        # UMAP by leiden clusters
        sc.pl.umap(adata, color='leiden', ax=axes[1], show=False, title='UMAP by Leiden Clusters')
        
        plt.tight_layout()
        plt.savefig(plots_dir / "umap_comparison.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Create batch mixing plots
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Calculate batch mixing metrics
        unique_batches = adata.obs['sample'].unique()
        if len(unique_batches) > 1:
            # Create batch mixing visualization
            batch_mixing = adata.obs.groupby(['leiden', 'sample']).size().unstack(fill_value=0)
            batch_mixing_pct = batch_mixing.div(batch_mixing.sum(axis=1), axis=0)
            
            sns.heatmap(batch_mixing_pct, annot=True, fmt='.2f', ax=ax)
            ax.set_title('Batch Distribution per Cluster')
            plt.tight_layout()
            plt.savefig(plots_dir / "batch_mixing_heatmap.png", dpi=300, bbox_inches='tight')
            plt.close()
        
        logging.info("Created integration plots")
    
    def process_all_samples(self, sample_list: list = None, 
                          do_hyperparameter_tuning: bool = True,
                          force: bool = False) -> Dict[str, Path]:
        """Process all samples through integration pipeline."""
        
        if sample_list is None:
            sample_list = [f.stem.replace('_filtered', '') 
                          for f in self.input_dir.glob('*_filtered.h5ad')]
        
        if DEBUG_MODE and len(sample_list) > 5:
            sample_list = sample_list[:5]
            logging.info(f"Debug mode: processing only {len(sample_list)} samples")
        
        if not sample_list:
            raise ValueError("No samples found")
        
        # Load and concatenate data
        adata = self.load_and_concatenate_data(sample_list)
        
        # Prepare data
        adata = self.prepare_data_for_scvi(adata)
        
        # Hyperparameter tuning or use defaults
        if do_hyperparameter_tuning:
            analysis = self.hyperparameter_search(adata)
            model = self.train_best_model(adata, analysis)
        else:
            model = self.train_scvi_model(adata, {
                'n_layers': 2,
                'n_latent': self.latent_dim,
                'n_hidden': 128,
                'max_epochs': self.max_epochs,
                'batch_size': self.batch_size,
                'lr': self.learning_rate
            })
        
        # Integrate data
        integrated_adata = self.integrate_data(model, adata)
        
        # Save integrated data
        output_file = self.output_dir / "integrated.h5ad"
        integrated_adata.write_h5ad(output_file)
        
        # Create plots
        self.create_integration_plots(integrated_adata)
        
        # Save summary
        summary = {
            'n_samples': len(sample_list),
            'n_cells': integrated_adata.n_obs,
            'n_genes': integrated_adata.n_vars,
            'n_clusters': len(integrated_adata.obs['leiden'].unique()),
            'latent_dim': integrated_adata.obsm["X_scVI"].shape[1]
        }
        
        summary_file = self.output_dir / "integration_summary.json"
        pd.Series(summary).to_json(summary_file)
        
        logging.info(f"Integration completed: {summary}")
        
        return {
            'integrated': output_file,
            'model': self.models_dir / "best_scvi_model",
            'summary': summary_file
        }


def main():
    """Main function for integration processing."""
    
    parser = argparse.ArgumentParser(description="scVI integration with hyperparameter tuning")
    parser.add_argument("--input", default="data/processed", help="Input directory")
    parser.add_argument("--output", default="results/integrated", help="Output directory")
    parser.add_argument("--samples", nargs="+", help="Specific samples to process")
    parser.add_argument("--no-tuning", action="store_true", help="Skip hyperparameter tuning")
    parser.add_argument("--latent-dim", type=int, default=SCVI_LATENT_DIM)
    parser.add_argument("--max-epochs", type=int, default=SCVI_MAX_EPOCHS)
    parser.add_argument("--batch-size", type=int, default=SCVI_BATCH_SIZE)
    parser.add_argument("--learning-rate", type=float, default=SCVI_LEARNING_RATE)
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging("integration.log")
    
    processor = IntegrationProcessor(
        input_dir=args.input,
        output_dir=args.output,
        latent_dim=args.latent_dim,
        max_epochs=args.max_epochs,
        batch_size=args.batch_size,
        learning_rate=args.learning_rate
    )
    
    results = processor.process_all_samples(
        sample_list=args.samples,
        do_hyperparameter_tuning=not args.no_tuning
    )
    
    logger.info("Integration completed successfully")


if __name__ == "__main__":
    main()