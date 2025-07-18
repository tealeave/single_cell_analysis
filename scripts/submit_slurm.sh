#!/bin/bash
#SBATCH --job-name=scRNA_pipeline
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=free

# Load required modules
module load python/3.10.2
module load nextflow/23.04.1

# Set up environment
export PYTHONPATH="$PWD/src:$PYTHONPATH"
export RAY_DISABLE_IMPORT_WARNING=1

# Create directories
mkdir -p logs reports results

# Check if debug mode requested
DEBUG_MODE=${1:-false}
SAMPLES=${2:-""}

# Set debug mode parameters
if [ "$DEBUG_MODE" = "true" ]; then
    echo "Running in DEBUG mode with 3 samples"
    DEBUG_FLAG="--debug true"
    TIME_LIMIT="4:00:00"
else
    echo "Running in PRODUCTION mode with all samples"
    DEBUG_FLAG="--debug false"
    TIME_LIMIT="24:00:00"
fi

# Update SLURM time limit if needed
if [ "$DEBUG_MODE" = "true" ]; then
    sbatch --time=$TIME_LIMIT --job-name=scRNA_debug --wrap="nextflow run nextflow/main.nf $DEBUG_FLAG --samples $SAMPLES"
    exit 0
fi

# Run the pipeline
echo "Starting Nextflow pipeline at $(date)"
echo "Debug mode: $DEBUG_MODE"
echo "Samples: $SAMPLES"

nextflow run nextflow/main.nf \
    -config nextflow/nextflow.config \
    $DEBUG_FLAG \
    --samples "$SAMPLES" \
    -with-report reports/pipeline_report.html \
    -with-timeline reports/timeline.html \
    -with-dag reports/pipeline_dag.html \
    -resume

echo "Pipeline completed at $(date)"