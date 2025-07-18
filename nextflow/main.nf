#!/usr/bin/env nextflow

/*
 * Single-Cell Analysis Pipeline for GSE235063
 * Nextflow workflow with debug mode support
 */

nextflow.enable.dsl=2

// Workflow parameters
params.input_dir = "data/raw"
params.output_dir = "results"
params.debug = false
params.samples = null
params.max_cpus = 32
params.max_memory = "128GB"
params.slurm_queue = "your-queue"

// Debug mode configuration
if (params.debug) {
    params.max_cpus = 4
    params.max_memory = "16GB"
    params.cellbender_cpus = 2
    params.cellbender_memory = "8GB"
    params.cellbender_time = "2h"
    params.qc_cpus = 2
    params.qc_memory = "4GB"
    params.qc_time = "1h"
    params.integration_cpus = 4
    params.integration_memory = "8GB"
    params.integration_time = "2h"
    params.annotation_cpus = 2
    params.annotation_memory = "4GB"
    params.annotation_time = "1h"
} else {
    params.cellbender_cpus = 8
    params.cellbender_memory = "32GB"
    params.cellbender_time = "8h"
    params.qc_cpus = 4
    params.qc_memory = "16GB"
    params.qc_time = "4h"
    params.integration_cpus = 8
    params.integration_memory = "16GB"
    params.integration_time = "4h"
    params.annotation_cpus = 4
    params.annotation_memory = "8GB"
    params.annotation_time = "2h"
}

// Sample discovery
samples_ch = Channel.fromPath("${params.input_dir}/*")
    .map { file -> file.name }
    .ifEmpty { error "No samples found in ${params.input_dir}" }

// Limit samples for debug mode
if (params.debug) {
    samples_ch = samples_ch.take(3)
    log.info "Debug mode: processing only 3 samples"
}

// If specific samples provided
if (params.samples) {
    samples_ch = Channel.from(params.samples.split(','))
}

/*
 * Process 1: CellBender Ambient RNA Removal
 */
process cellbender {
    tag "${sample_id}"
    
    cpus params.cellbender_cpus
    memory params.cellbender_memory
    time params.cellbender_time
    
    input:
    val sample_id
    
    output:
    path "${sample_id}_denoised.h5ad", emit: denoised
    path "${sample_id}_cellbender_metrics.csv", emit: metrics
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "Processing ${sample_id} with CellBender"
    
    # Run CellBender through Python wrapper
    uv run python ${baseDir}/../src/preprocessing.py \
        --input ${params.input_dir} \
        --output ${params.output_dir}/cellbender \
        --sample ${sample_id} \
        --expected-cells 5000 \
        --total-droplets 50000 \
        --fpr 0.01
    
    # Move outputs to expected locations
    mv ${params.output_dir}/cellbender/${sample_id}_denoised.h5ad .
    mv ${params.output_dir}/cellbender/${sample_id}_cellbender_metrics.csv .
    """
    
    stub:
    """
    touch ${sample_id}_denoised.h5ad
    touch ${sample_id}_cellbender_metrics.csv
    """
}

/*
 * Process 2: Quality Control and Filtering
 */
process qc_filter {
    tag "${sample_id}"
    
    cpus params.qc_cpus
    memory params.qc_memory
    time params.qc_time
    
    input:
    path denoised_file
    val sample_id
    
    output:
    path "${sample_id}_filtered.h5ad", emit: filtered
    path "${sample_id}_qc_report.html", emit: qc_report
    path "${sample_id}_filtering_report.csv", emit: filtering_report
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "QC analysis for ${sample_id}"
    
    uv run python ${baseDir}/../src/qc_filter.py \
        --input ${params.output_dir}/cellbender \
        --output ${params.output_dir}/qc \
        --sample ${sample_id} \
        --min-genes 200 \
        --max-genes 6000 \
        --max-mt-percent 20 \
        --mad-threshold 3
    
    # Move outputs to expected locations
    mv ${params.output_dir}/qc/${sample_id}_filtered.h5ad .
    mv ${params.output_dir}/qc/qc_reports/${sample_id}_qc_report.html .
    mv ${params.output_dir}/qc/qc_reports/${sample_id}_filtering_report.csv .
    """
    
    stub:
    """
    touch ${sample_id}_filtered.h5ad
    touch ${sample_id}_qc_report.html
    touch ${sample_id}_filtering_report.csv
    """
}

/*
 * Process 3: Integration with scVI
 */
process integration {
    
    cpus params.integration_cpus
    memory params.integration_memory
    time params.integration_time
    
    input:
    path filtered_files
    
    output:
    path "integrated.h5ad", emit: integrated
    path "best_scvi_model", emit: model
    path "integration_summary.json", emit: summary
    path "integration_plots", emit: plots
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "Running scVI integration"
    
    # Collect all filtered files
    mkdir -p temp_input
    cp ${filtered_files} temp_input/
    
    uv run python ${baseDir}/../src/integration.py \
        --input temp_input \
        --output ${params.output_dir}/integrated \
        ${params.debug ? '--no-tuning' : ''}
    
    # Move outputs to expected locations
    mv ${params.output_dir}/integrated/integrated.h5ad .
    mv ${params.output_dir}/integrated/models/best_scvi_model .
    mv ${params.output_dir}/integrated/integration_summary.json .
    mv ${params.output_dir}/integrated/integration_plots .
    """
    
    stub:
    """
    touch integrated.h5ad
    mkdir -p best_scvi_model
    echo '{}' > integration_summary.json
    mkdir -p integration_plots
    """
}

/*
 * Process 4: Annotation
 */
process annotation {
    
    cpus params.annotation_cpus
    memory params.annotation_memory
    time params.annotation_time
    
    input:
    path integrated_file
    
    output:
    path "annotated.h5ad", emit: annotated
    path "annotation_report.html", emit: report
    path "validation_report.csv", emit: validation
    path "annotation_plots", emit: plots
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "Running annotation pipeline"
    
    uv run python ${baseDir}/../src/annotation.py \
        --input ${params.output_dir}/integrated \
        --output ${params.output_dir}/annotations
    
    # Move outputs to expected locations
    mv ${params.output_dir}/annotations/annotated.h5ad .
    mv ${params.output_dir}/annotations/annotation_report.html .
    mv ${params.output_dir}/annotations/validation_report.csv .
    mv ${params.output_dir}/annotations/annotation_plots .
    """
    
    stub:
    """
    touch annotated.h5ad
    touch annotation_report.html
    touch validation_report.csv
    mkdir -p annotation_plots
    """
}

/*
 * Workflow definition
 */
workflow {
    
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║ Single-Cell Analysis Pipeline for GSE235063              ║
    ║ Debug Mode: ${params.debug}                              ║
    ║ Samples: ${params.samples ?: 'All available'}            ║
    ╚══════════════════════════════════════════════════════════╝
    """
    
    // Main workflow
    cellbender_out = cellbender(samples_ch)
    
    qc_out = qc_filter(cellbender_out.denoised, cellbender_out.denoised.map { it.name.replace('_denoised.h5ad', '') }.flatten())
    
    // Collect all filtered files for integration
    filtered_files = qc_out.filtered.collect()
    
    integration_out = integration(filtered_files)
    
    annotation_out = annotation(integration_out.integrated)
    
    // Create final summary
    create_summary_report(
        cellbender_out.metrics.collect(),
        qc_out.filtering_report.collect(),
        integration_out.summary,
        annotation_out.validation
    )
}

/*
 * Process to create final summary report
 */
process create_summary_report {
    
    publishDir "${params.output_dir}/reports", mode: 'copy'
    
    input:
    path cellbender_metrics
    path qc_reports
    path integration_summary
    path validation_report
    
    output:
    path "pipeline_summary.html"
    path "summary_stats.json"
    
    script:
    """
    #!/bin/bash
    set -e
    
    echo "Creating final summary report"
    
    # Create summary HTML report
    cat > pipeline_summary.html << 'EOF'
    <!DOCTYPE html>
    <html>
    <head>
        <title>Single-Cell Analysis Pipeline Summary</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 40px; }
            .header { background: #f0f0f0; padding: 20px; border-radius: 5px; }
            .section { margin: 20px 0; }
            .metrics { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; }
            .metric-card { border: 1px solid #ddd; padding: 15px; border-radius: 5px; }
            .success { color: #28a745; }
            .warning { color: #ffc107; }
            .error { color: #dc3545; }
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Single-Cell Analysis Pipeline Summary</h1>
            <p>Analysis completed successfully</p>
            <p>Debug mode: DEBUG_MODE_PLACEHOLDER</p>
        </div>
        
        <div class="section">
            <h2>Pipeline Status</h2>
            <div class="success">✅ All processes completed successfully</div>
        </div>
        
        <div class="section">
            <h2>Output Files</h2>
            <ul>
                <li>Annotated dataset: annotated.h5ad</li>
                <li>Integration model: best_scvi_model/</li>
                <li>QC reports: Available in results/reports/</li>
                <li>Visualization plots: Available in results/*/plots/</li>
            </ul>
        </div>
        
        <div class="section">
            <h2>Next Steps</h2>
            <ul>
                <li>Review annotation results in annotated.h5ad</li>
                <li>Examine QC plots for data quality</li>
                <li>Run downstream analysis as needed</li>
            </ul>
        </div>
    </body>
    </html>
    EOF
    
    # Create summary JSON with current timestamp
    python3 -c "
import json
import datetime
with open('summary_stats.json', 'w') as f:
    json.dump({
        'status': 'completed',
        'timestamp': datetime.datetime.now().isoformat(),
        'debug_mode': '${params.debug}'
    }, f, indent=2)
"
    """
}

/*
 * Helper functions
 */
def getRuntimeConfig() {
    return [
        debug: params.debug,
        cellbender_resources: [cpus: params.cellbender_cpus, memory: params.cellbender_memory, time: params.cellbender_time],
        qc_resources: [cpus: params.qc_cpus, memory: params.qc_memory, time: params.qc_time],
        integration_resources: [cpus: params.integration_cpus, memory: params.integration_memory, time: params.integration_time],
        annotation_resources: [cpus: params.annotation_cpus, memory: params.annotation_memory, time: params.annotation_time]
    ]
}