# Single-Cell Analysis Pipeline Usage Guide

## Quick Start

### 1. Initial Setup
```bash
# Validate setup
./scripts/validate_setup.sh

# Install dependencies (if needed)
uv sync --no-cache
```

### 2. Testing Pipeline (Recommended First Step)
```bash
# Run with 3 samples for testing (2-3 hours)
./scripts/test_pipeline.sh
```

### 3. Production Run
```bash
# Full pipeline with all 75 samples (1-2 days)
sbatch scripts/submit_slurm.sh false
```

### 4. Debug Mode
```bash
# Quick debug run with 3 samples
sbatch scripts/submit_slurm.sh true
```

## Pipeline Overview

The pipeline processes GSE235063 (75 AML bone marrow samples) through 4 stages:

1. **CellBender** - Ambient RNA removal
2. **QC & Filtering** - Quality control with MAD outlier detection
3. **Integration** - scVI with hyperparameter tuning
4. **Annotation** - CellTypist + scVI label transfer

## Resource Requirements

| Mode | Samples | RAM | Time | CPUs |
|------|---------|-----|------|------|
| Debug | 3 | 4-8GB | 2-3 hours | 4 |
| Production | 75 | 32-128GB | 1-2 days | 32 |

## Directory Structure

```
single_cell_analysis/
├── data/
│   ├── raw/           # Raw 10x data
│   └── processed/     # Processed outputs
├── results/
│   ├── cellbender/    # CellBender outputs
│   ├── qc/           # QC reports
│   ├── integrated/   # scVI integration
│   └── annotations/  # Final annotations
├── nextflow/
│   ├── main.nf       # Main workflow
│   └── nextflow.config # SLURM config
├── src/
│   ├── preprocessing.py # CellBender module
│   ├── qc_filter.py     # QC module
│   ├── integration.py   # scVI module
│   └── annotation.py    # Annotation module
└── scripts/
    ├── test_pipeline.sh   # Test runner
    ├── submit_slurm.sh    # SLURM submission
    └── validate_setup.sh  # Setup validation
```

## Commands Reference

### Individual Module Testing
```bash
# Test CellBender
uv run python src/preprocessing.py --sample GSM123456 --debug

# Test QC
uv run python src/qc_filter.py --sample GSM123456 --debug

# Test Integration (reduced)
uv run python src/integration.py --no-tuning --max-epochs 50

# Test Annotation
uv run python src/annotation.py --model Immune_All_Low --debug
```

### Nextflow Commands
```bash
# Run specific samples
nextflow run nextflow/main.nf --samples GSM123456,GSM123457,GSM123458

# Resume interrupted run
nextflow run nextflow/main.nf -resume

# Generate reports
nextflow run nextflow/main.nf -with-report -with-timeline -with-dag
```

## SLURM Integration

### Submit Jobs
```bash
# Debug mode (3 samples)
sbatch scripts/submit_slurm.sh true

# Production mode (all samples)
sbatch scripts/submit_slurm.sh false

# Custom samples
sbatch scripts/submit_slurm.sh false GSM123456,GSM123457
```

### Monitor Jobs
```bash
squeue -u $USER
sacct -j JOBID --format=JobID,JobName%20,State,Elapsed,MaxRSS
```

## Output Files

### Key Outputs
- `results/annotations/annotated.h5ad` - Final annotated dataset
- `results/integrated/best_scvi_model/` - Trained scVI model
- `results/reports/` - HTML reports and visualizations

### Reports
- `pipeline_report.html` - Comprehensive pipeline report
- `timeline.html` - Execution timeline
- `pipeline_dag.html` - Workflow DAG
- `trace.txt` - Resource usage trace

## Debug Mode Features

- **Reduced samples**: 3 samples instead of 75
- **Faster processing**: 2-3 hours vs 1-2 days
- **Lower resources**: 4-8GB RAM vs 32-128GB
- **Fewer epochs**: 10-50 vs 150-400
- **Limited tuning**: 5 trials vs 50 trials

## Troubleshooting

### Common Issues

1. **Memory errors**: Use debug mode or increase SLURM memory
2. **CellBender not found**: Install with `pip install cellbender`
3. **Sample not found**: Check data/raw/ directory structure
4. **Permission errors**: Ensure scripts are executable

### Debug Commands
```bash
# Check job status
squeue -u $USER

# Check logs
tail -f logs/slurm_JOBID.out
tail -f logs/slurm_JOBID.err

# Test individual components
uv run python -c "import src.config; print('Config loaded')"
```

## Performance Tips

1. **Use debug mode** for initial testing
2. **Resume capability** - add `-resume` to continue interrupted runs
3. **Resource scaling** - adjust SLURM parameters in nextflow.config
4. **Monitor memory** - check trace.txt for resource usage

## Example Workflow

```bash
# 1. Setup
./scripts/validate_setup.sh

# 2. Test with 3 samples
./scripts/test_pipeline.sh

# 3. If successful, run production
sbatch scripts/submit_slurm.sh false

# 4. Monitor progress
squeue -u $USER
ls -la results/annotations/
```

## Contact and Support

For issues or questions:
- Check logs in `logs/` directory
- Review reports in `results/reports/`
- Validate setup with `./scripts/validate_setup.sh`