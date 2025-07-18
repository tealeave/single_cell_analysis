#!/bin/bash

# Validation script for pipeline setup

echo "🔍 Validating pipeline setup..."

# Check Python and dependencies
echo "📋 Checking Python environment..."
python --version 2>/dev/null || echo "❌ Python not found"
python -c "import sys; print(f'Python {sys.version}')" 2>/dev/null

# Check UV installation
echo "📋 Checking UV..."
uv --version 2>/dev/null || echo "⚠️  UV not found, using system Python"

# Check required packages
echo "📋 Checking required packages..."
python -c "import scanpy; print(f'✅ scanpy {scanpy.__version__}')" 2>/dev/null || echo "❌ scanpy not found"
python -c "import scvi; print(f'✅ scvi {scvi.__version__}')" 2>/dev/null || echo "❌ scvi not found"
python -c "import celltypist; print(f'✅ celltypist {celltypist.__version__}')" 2>/dev/null || echo "❌ celltypist not found"

# Check Nextflow
echo "📋 Checking Nextflow..."
nextflow -version 2>/dev/null || echo "❌ Nextflow not found"

# Check directory structure
echo "📋 Checking directory structure..."
[ -d "src" ] && echo "✅ src/ directory exists" || echo "❌ src/ missing"
[ -d "data" ] && echo "✅ data/ directory exists" || echo "❌ data/ missing"
[ -d "nextflow" ] && echo "✅ nextflow/ directory exists" || echo "❌ nextflow/ missing"
[ -d "scripts" ] && echo "✅ scripts/ directory exists" || echo "❌ scripts/ missing"

# Check data availability
echo "📋 Checking data availability..."
if [ -d "data/raw" ]; then
    SAMPLE_COUNT=$(ls data/raw/ | wc -l)
    echo "✅ Found $SAMPLE_COUNT samples in data/raw/"
    
    if [ $SAMPLE_COUNT -ge 3 ]; then
        echo "✅ Sufficient samples for testing"
    else
        echo "⚠️  Fewer than 3 samples found"
    fi
else
    echo "❌ data/raw/ directory not found"
fi

# Check configuration files
echo "📋 Checking configuration files..."
[ -f "nextflow/main.nf" ] && echo "✅ Nextflow main workflow exists" || echo "❌ nextflow/main.nf missing"
[ -f "nextflow/nextflow.config" ] && echo "✅ Nextflow config exists" || echo "❌ nextflow.config missing"
[ -f "src/config.py" ] && echo "✅ Python config exists" || echo "❌ src/config.py missing"

# Check script permissions
echo "📋 Checking script permissions..."
[ -x "scripts/test_pipeline.sh" ] && echo "✅ test_pipeline.sh is executable" || echo "⚠️  test_pipeline.sh not executable"
[ -x "scripts/submit_slurm.sh" ] && echo "✅ submit_slurm.sh is executable" || echo "⚠️  submit_slurm.sh not executable"

# Check CellBender availability
echo "📋 Checking CellBender..."
cellbender --help >/dev/null 2>&1 && echo "✅ CellBender available" || echo "⚠️  CellBender not found - install with: pip install cellbender"

echo ""
echo "🎯 Setup validation complete!"
echo ""
echo "🚀 Next steps:"
echo "   1. Run test: ./scripts/test_pipeline.sh"
echo "   2. Submit to SLURM: sbatch scripts/submit_slurm.sh true"
echo "   3. Run full pipeline: sbatch scripts/submit_slurm.sh false"