# Single-Cell Analysis Pipeline for GSE235063

A production-grade, reproducible single-cell RNA-seq analysis pipeline for the GSE235063 dataset (75 AML samples) with ambient RNA removal, QC, integration, and annotation capabilities.

## 🎯 Project Overview

This pipeline processes 75 AML bone marrow samples from GSE235063 using modern single-cell analysis tools:
- **Ambient RNA removal**: CellBender (CPU mode)
- **Quality control**: MAD outlier detection + doublet removal
- **Integration**: scVI-tools with hyperparameter tuning
- **Annotation**: CellTypist + scVI label transfer

## 🚀 Quick Start

### Debug Mode (3 samples, 2-3 hours)
```bash
# Activate environment
module load python/3.10.2
source .venv/bin/activate

# Run debug pipeline
export DEBUG_MODE=true
uv run nextflow run main.nf --debug true --samples 3
```

### Production Mode (75 samples, 1-2 days)
```bash
# Full pipeline
export DEBUG_MODE=false
uv run nextflow run main.nf -profile slurm --resume
```

## 📁 Project Structure

```
single-cell-analysis/
├── nextflow/               # Nextflow workflow files
│   ├── main.nf            # Core pipeline
│   ├── nextflow.config    # SLURM configuration
│   ├── modules/           # Pipeline modules
│   └── bin/              # Custom scripts
├── src/                   # Python analysis modules
│   ├── preprocessing.py   # CellBender ambient RNA removal
│   ├── qc_filter.py      # QC and filtering
│   ├── integration.py    # scVI integration
│   ├── annotation.py     # Cell annotation
│   ├── config.py         # Configuration settings
│   └── utils.py          # Utility functions
├── data/                  # Data directories
│   ├── raw/              # Raw GEO data
│   ├── processed/        # Processed AnnData files
│   └── references/       # Reference datasets
├── results/              # Output directories
│   ├── qc_reports/       # QC visualizations
│   ├── integrated/       # Integrated datasets
│   └── annotations/      # Annotated datasets
└── logs/                 # Execution logs
```

## 🛠️ Installation

### HPC Setup (Linux)
```bash
# Clone repository
git clone <repository-url>
cd single-cell-analysis

# Load Python module
module load python/3.10.2

# Install UV (if not available)
pip install uv

# Create environment
uv sync  # Uses pyproject.toml

# Activate environment
source .venv/bin/activate
```

### Dependencies
All dependencies are managed via UV and locked in `uv.lock`:
- **Core**: scanpy, scvi-tools, celltypist, cellbender
- **Analysis**: doubletdetection, seaborn, matplotlib
- **HPC**: ray[tune], nextflow
- **Environment**: Python 3.10.2 (locked)

## 📊 Data Pipeline

### 1. Data Acquisition
```bash
# Download GSE235063
python src/download_data.py --study GSE235063 --output data/raw/

# Validate downloads
python src/validate_data.py --input data/raw/
```

### 2. Ambient RNA Removal
- **Tool**: CellBender (CPU mode)
- **Expected cells**: 5,000 per sample
- **Total droplets**: 50,000 (debug: 10,000)
- **Output**: Denoised AnnData files

### 3. Quality Control
- **Metrics**: n_genes, n_counts, mt_percent
- **Outliers**: MAD-based detection (threshold=3)
- **Doublets**: DoubletDetection
- **Filters**: min_genes=200, max_genes=6000, max_mt_percent=20

### 4. Integration
- **Tool**: scVI-tools
- **Hyperparameters**: 
  - Latent dimensions: 30
  - Hidden layers: 2
  - Hidden units: 128
  - Max epochs: 400 (debug: 50)
- **Tuning**: Ray Tune with 50 trials

### 5. Annotation
- **Tools**: CellTypist + scVI label transfer
- **Models**: 
  - Custom AML models
  - Immune_All_Low
- **Validation**: Manual review of markers

## 🔧 Configuration

### Debug vs Production
| Setting | Debug | Production |
|---------|--------|------------|
| Samples | 3 | 75 |
| CellBender memory | 8GB | 32GB |
| CellBender time | 2h | 8h |
| scVI epochs | 50 | 400 |
| QC memory | 4GB | 16GB |

### Environment Variables
```bash
export DEBUG_MODE=true/false     # Debug mode toggle
export SLURM_QUEUE=your-queue    # SLURM queue name
export MAX_CPUS=32               # Max CPU cores
export MAX_MEMORY=128G           # Max memory
```

## ⚙️ Nextflow Execution

### Local Testing
```bash
# Test with small dataset
nextflow run main.nf -with-report debug_report.html

# Resume from checkpoint
nextflow run main.nf -resume
```

### HPC SLURM
```bash
# Submit to SLURM
sbatch nextflow_slurm.sh

# Monitor jobs
squeue -u $USER
```

## 📈 Output Files

### Per-Sample Outputs
- `{sample}_denoised.h5ad`: CellBender output
- `{sample}_qc_report.html`: QC visualizations
- `{sample}_filtered.h5ad`: QC-filtered data

### Integrated Outputs
- `integrated.h5ad`: scVI-integrated dataset
- `annotated.h5ad`: Cell-type annotated dataset
- `qc_summary.html`: Pipeline summary report

### Reports
- `timeline.html`: Execution timeline
- `report.html`: Detailed execution report
- `dag.png`: Workflow DAG

## 🧪 Development

### Adding New Modules
1. Create module in `src/`
2. Add Nextflow process in `nextflow/modules/`
3. Update configuration in `src/config.py`
4. Add tests in `tests/`

### Debug Testing
```bash
# Quick syntax validation
nextflow run nextflow/main.nf --debug true -stub

# Full debug test with 3 samples
./scripts/test_debug.sh

# Monitor debug run
nextflow run nextflow/main.nf --debug true --with-report debug_report.html
```

## 🔍 Troubleshooting

### Common Issues
- **Memory errors**: Reduce batch size in debug mode
- **SLURM timeouts**: Increase time limits in config
- **Missing dependencies**: Run `uv sync --reinstall`
- **Nextflow resume**: Use `-resume` flag for interrupted runs

### Debug Commands
```bash
# Check sample manifest
python src/utils.py --manifest data/raw/

# Validate environment
python -c "import scanpy; print(scanpy.__version__)"

# Test CellBender
uv run cellbender remove-background --help
```

## 📚 References

- **GSE235063**: [GEO Study](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235063)
- **CellBender**: [Documentation](https://cellbender.readthedocs.io/)
- **scVI-tools**: [Documentation](https://docs.scvi-tools.org/)
- **CellTypist**: [Documentation](https://celltypist.readthedocs.io/)

## 🤝 Contributing

1. Create feature branch: `git checkout -b feature/new-module`
2. Add tests for new functionality
3. Run quality checks: `black src/ && isort src/`
4. Submit pull request with clear description

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.