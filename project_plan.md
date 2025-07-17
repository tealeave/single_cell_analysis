### Comprehensive Project Plan for Single-Cell Analysis Pipeline

As a senior principal data scientist with over 10 years of experience in bioinformatics and single-cell genomics, I'll outline a robust, reproducible, and scalable plan for your single-cell analysis project. This plan draws from the provided tutorial transcripts (focusing on scRNA-seq preprocessing, QC, ambient RNA removal, integration, and annotation) while incorporating best practices for modern workflows. Key principles guiding this plan:

- **Reproducibility**: Use UV for dependency management (locked to Python 3.10.2), Git for version control, and Nextflow for orchestration to ensure resumability, provenance tracking, and environment-agnostic execution.
- **Modularity**: Break the pipeline into discrete, testable stages (e.g., preprocessing per sample, integration across samples) to allow easy debugging and reuse.
- **Scalability and HPC Optimization**: Leverage SLURM for parallelization on your HPC cluster. Identify embarrassingly parallel steps (e.g., processing multiple files/samples) and distribute them as array jobs or Nextflow processes. No GPU assumptions, as per your setup—use CPU-friendly tools like Scanpy, scVI-tools (CPU mode), and CellBender (CPU fallback).
- **Resumability**: Nextflow's built-in caching and resume features will handle interruptions, similar to your RNA-seq example. Use `--resume` flags and configure SLURM executors.
- **Error Handling and Quality**: Include validation checks, logging (via Nextflow reports and Python logging), and automated QC reports. Adhere to FAIR principles (Findable, Accessible, Interoperable, Reusable) for data and code.
- **Efficiency**: Minimize manual intervention; automate where possible. Handle large datasets (e.g., multiple samples) by processing in batches.
- **Extension to Multi-Modal**: The plan focuses on scRNA-seq (as per transcripts), but includes hooks for scATAC-seq integration (e.g., via Signac or ArchR) in future phases.
- **Testing and Iteration**: Start with a small subset of data for local testing before scaling to HPC.

The pipeline will cover:
- Data download and preparation.
- Ambient RNA removal (CellBender, CPU mode).
- QC and filtering (iterative, per-sample).
- Integration (scVI-tools with hyperparameter tuning).
- Annotation (CellTypist with custom models, scVI label transfer).
- Downstream: Differential expression/abundance (placeholder for Part 3 if available).

Estimated timeline: 2-4 weeks for setup and initial run (assuming 1-2 FTEs), depending on data size and HPC queue times.

---

#### 1. Project Setup and Environment Management
**Goal**: Create a reproducible, version-controlled project structure.

- **Directory Structure**:
  ```
  single-cell-project/
  ├── nextflow/               # Nextflow workflow files
  │   ├── main.nf             # Core pipeline
  │   ├── nextflow.config     # Config for SLURM executor, resources
  │   ├── modules/            # Modular processes (e.g., qc.nf, integrate.nf)
  │   └── bin/                # Custom scripts (e.g., Python wrappers)
  ├── src/                    # Python scripts/notebooks
  │   ├── preprocessing.py    # Core preprocessing script
  │   ├── qc_filter.py        # QC and filtering
  │   ├── integration.py      # scVI integration
  │   ├── annotation.py       # CellTypist/scVI annotation
  │   └── utils.py            # Helper functions (e.g., MAD outlier detection)
  ├── data/                   # Raw and intermediate data (symlink to HPC storage if needed)
  │   ├── raw/                # Downloaded GEO data
  │   ├── processed/          # Filtered matrices
  │   └── references/         # Reference datasets/models
  ├── results/                # Outputs (QC plots, AnnData files, reports)
  │   ├── qc_reports/         # Per-sample QC visualizations
  │   ├── integrated/         # Integrated AnnData
  │   └── annotations/        # Labeled datasets
  ├── logs/                   # Nextflow/SLURM logs
  ├── uv.lock                 # UV lockfile for dependencies
  ├── README.md               # Documentation, usage instructions
  └── .gitignore              # Ignore large data files, results
  ```

- **Version Control**: Initialize with Git (`git init`). Track code/scripts; ignore large data/results. Use branches (e.g., `dev` for testing, `main` for stable). Commit often with semantic messages (e.g., "feat: add QC module").

- **Dependency Management with UV**:
  - Install UV: `pip install uv` (or via Homebrew/PyPI).
  - Create project: `uv init --python 3.10.2`.
  - Add dependencies: `uv add scanpy scvi-tools celltypist doubletdetection seaborn matplotlib pandas numpy scipy`.
  - For HPC: Generate a `requirements.txt` with `uv export --format requirements-txt > requirements.txt` for SLURM jobs.
  - Lockfile: Run `uv sync` to generate `uv.lock` for exact reproducibility.
  - Activation: Use `uv run` for scripts (e.g., `uv run python src/preprocessing.py`).

- **Documentation**: README.md with setup instructions, pipeline overview, SLURM/Nextflow usage. Include Jupyter notebooks in `src/` for exploratory analysis (e.g., QC visualization).

---

#### 2. Data Acquisition and Preparation
**Goal**: Download, organize, and prepare data for pipeline input.

- **Steps**:
  1. Download GEO data (GSE235063) via script: Use `wget` or `aria2` for parallel downloads in Nextflow.
  2. Extract and rename files (as in tutorial: rename genes to features, add "Gene Expression" column).
  3. Convert to AnnData: Use Scanpy's `read_10x_mtx` in a Python script; save as `.h5ad`.
  4. Reference Data: Download AML/bone marrow references (as in transcript). Process into custom CellTypist models.

- **Nextflow Integration**: A `download.nf` process to fetch/extract data, resumable if interrupted.
- **HPC**: Run as a single SLURM job (low resource: 1 core, 4GB RAM).

- **Best Practices**: Validate downloads with MD5 checksums. Store raw data immutably; use symlinks for large files on HPC shared storage.

---

#### 3. Pipeline Stages
The pipeline is modular, with each stage as a Nextflow process. Inputs/outputs are AnnData (.h5ad) files for resumability.

- **Stage 1: Ambient RNA Removal (CellBender)**:
  - Input: Raw .h5ad per sample.
  - Script: Wrapper around CellBender (CPU mode: `--total-droplets-included 50000`).
  - Parallelization: Per-sample (75 samples → SLURM array jobs).
  - Output: Denoised .h5ad + metrics.csv.
  - QC: Generate histograms of fraction_counts_removed (as in transcript).

- **Stage 2: Quality Control and Filtering**:
  - Input: Denoised .h5ad per sample.
  - Script: `qc_filter.py` – Compute metrics (mt%, n_genes, etc.), MAD outlier detection, DoubletDetection.
  - Iterative: Annotate outliers/doublets but filter conservatively; visualize distributions (Ridge plots).
  - Parallelization: Per-sample.
  - Output: Filtered .h5ad per sample + QC plots/reports.

- **Stage 3: Integration (scVI-tools)**:
  - Input: Filtered .h5ad (concatenated across samples).
  - Script: `integration.py` – Hyperparameter tuning (Ray Tune: search space for layers, latent dims, etc.), train scVI model.
  - Tuning: Use CPU (no GPU); limit trials to 50-100 for efficiency.
  - Output: Integrated AnnData with embeddings (X_scVI).
  - Best Practice: Early stopping; save model for reproducibility.

- **Stage 4: Annotation**:
  - Input: Integrated AnnData.
  - Script: `annotation.py` – CellTypist (custom AML models + Immune_All_Low), scVI label transfer.
  - Manual Verification: Overcluster (Leiden res=2.0), score gene sets (e.g., AML blasts), visualize UMAPs/markers.
  - Output: Annotated AnnData (.h5ad).

- **Stage 5: Downstream Analysis (Placeholder)**:
  - Differential expression/abundance (e.g., Milo or DESeq2-like via Scanpy).
  - Multi-modal hooks: Add processes for scATAC (e.g., Signac integration).

- **Reporting**: Nextflow HTML reports + custom Jupyter notebook for visualizations (e.g., UMAPs, Ridge plots).

---

#### 4. HPC and SLURM Integration
- **Parallelizable Steps**: Preprocessing/QC/ambient removal per sample (array jobs: `#SBATCH --array=0-74`).
- **Resources per Job**: 4-8 cores, 16-32GB RAM, 4-24 hours (tune via Nextflow).
- **SLURM Config**: In `nextflow.config`:
  ```
  process {
    executor = 'slurm'
    queue = 'your-queue'
    cpus = 4
    memory = '16 GB'
    time = '4h'
  }
  ```
- **Job Submission**: Nextflow handles SLURM submissions (`nextflow run main.nf -profile slurm`).
- **Data Management**: Use shared HPC storage; Nextflow stages files automatically.

---

#### 5. Nextflow Orchestration
**Goal**: Resumable, portable pipeline like your RNA-seq example.

- **main.nf** Outline:
  ```
  params.input_dir = 'data/raw'
  params.output_dir = 'results'

  process download_data { ... }  // wget + extract
  process ambient_removal { input: path(file); output: path(denoised.h5ad); script: "cellbender ..." }
  process qc_filter { input: path(denoised); output: path(filtered.h5ad); script: "python qc_filter.py ..." }
  process integrate { input: path('filtered/*.h5ad'); output: path(integrated.h5ad); script: "python integration.py ..." }
  process annotate { input: path(integrated); output: path(annotated.h5ad); script: "python annotation.py ..." }

  workflow {
    raw_files = Channel.fromPath("${params.input_dir}/*.h5ad")
    denoised = ambient_removal(raw_files)
    filtered = qc_filter(denoised)
    integrated = integrate(filtered.collect())  // Collect all for integration
    annotated = annotate(integrated)
  }
  ```
- **Config**: SLURM executor, resource limits, resume enabled.
- **Execution**: `nextflow run main.nf --resume` on HPC. Reports: Timeline, DAG, execution logs.
- **Best Practices**: Containerize with Singularity (Nextflow supports it for HPC). Use params for flexibility (e.g., `--samples 10` for testing).

---

#### 6. Tools and Dependencies
- **Core**: Scanpy, scVI-tools, CellTypist, DoubletDetection.
- **Viz/QC**: Seaborn, Matplotlib.
- **Utils**: Pandas, NumPy, SciPy.
- **Tuning/HPC**: Ray[Tune], Nextflow.
- **UV Lock**: Pin to Python 3.10.2; include versions for all packages.

---

#### 7. Implementation Steps and Timeline
1. **Week 1: Setup (1-2 days)**: Git init, UV environment, directory structure, basic Nextflow skeleton.
2. **Week 1-2: Core Scripts (3-5 days)**: Develop Python scripts based on transcripts (preprocessing, QC, etc.). Test locally on subset data.
3. **Week 2: Nextflow Integration (2-3 days)**: Modularize into processes; test resume/SLURM on small jobs.
4. **Week 3: Full Run and Validation (3-5 days)**: Run on HPC; generate QC reports; iterate on annotation.
5. **Week 4: Extensions/Refinement (2-4 days)**: Add scATAC hooks; document; handle errors (e.g., retries in Nextflow).

**Risks/Mitigations**: High memory usage—process in batches. SLURM queues—use `--maxForks` in Nextflow. Annotation accuracy—manual review step.

This plan ensures a professional, scalable pipeline. If you provide more details (e.g., data size, specific analyses), I can refine it further!