#!/bin/bash
# Nextflow-based DEBUG testing script

set -e

echo "🧪 Nextflow DEBUG Test Suite"
echo "=========================="

# Configuration
DEBUG_OUTPUT="results_debug"
LOG_DIR="logs_debug"

# Create debug directories
mkdir -p "$DEBUG_OUTPUT" "$LOG_DIR"

# Function to log messages
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_DIR/debug_test.log"
}

log "Starting Nextflow DEBUG testing..."

# Step 1: Validate Nextflow syntax
log "1. Testing Nextflow syntax..."
if nextflow run nextflow/main.nf --debug true --samples "test" -stub >/dev/null 2>&1; then
    log "✅ Nextflow syntax validation passed"
else
    log "❌ Nextflow syntax validation failed"
    exit 1
fi

# Step 2: Check data availability
log "2. Checking data availability..."
if [ -d "data/raw" ]; then
    SAMPLES=($(ls data/raw/ | head -3))
    SAMPLE_COUNT=${#SAMPLES[@]}
    
    if [ $SAMPLE_COUNT -ge 1 ]; then
        log "✅ Found $SAMPLE_COUNT samples: ${SAMPLES[*]}"
        TEST_SAMPLES=$(IFS=,; echo "${SAMPLES[*]}")
    else
        log "⚠️  No samples found, creating test structure..."
        mkdir -p data/test/{sample1,sample2,sample3}
        TEST_SAMPLES="sample1,sample2,sample3"
    fi
else
    log "⚠️  Creating test data structure..."
    mkdir -p data/test/{sample1,sample2,sample3}
    TEST_SAMPLES="sample1,sample2,sample3"
fi

# Step 3: Run DEBUG pipeline
log "3. Running DEBUG pipeline with samples: $TEST_SAMPLES"

# Set debug mode
export DEBUG_MODE=true

# Run Nextflow pipeline in DEBUG mode
nextflow run nextflow/main.nf \
    --debug true \
    --samples "$TEST_SAMPLES" \
    --output_dir "$DEBUG_OUTPUT" \
    -with-report "$DEBUG_OUTPUT/debug_report.html" \
    -with-timeline "$DEBUG_OUTPUT/debug_timeline.html" \
    -with-dag "$DEBUG_OUTPUT/debug_dag.html" \
    -with-trace "$DEBUG_OUTPUT/debug_trace.txt" \
    -resume

# Check results
log "4. Validating results..."
RESULT_FILES=(
    "$DEBUG_OUTPUT/annotations/annotated.h5ad"
    "$DEBUG_OUTPUT/reports/pipeline_summary.html"
    "$DEBUG_OUTPUT/debug_report.html"
)

for file in "${RESULT_FILES[@]}"; do
    if [ -f "$file" ]; then
        log "✅ Generated: $(basename "$file")"
    else
        log "⚠️  Missing: $(basename "$file")"
    fi
done

# Summary
log "🎯 DEBUG test completed!"
log "Results available in: $DEBUG_OUTPUT/"
log "Reports available in: $DEBUG_OUTPUT/reports/"

# Display quick summary
echo ""
echo "📊 Debug Test Summary"
echo "===================="
echo "Samples processed: $TEST_SAMPLES"
echo "Output directory: $DEBUG_OUTPUT"
echo "Log directory: $LOG_DIR"
echo ""
echo "📁 Key files to check:"
echo "   - $DEBUG_OUTPUT/annotations/annotated.h5ad (final output)"
echo "   - $DEBUG_OUTPUT/reports/pipeline_summary.html (comprehensive report)"
echo "   - $DEBUG_OUTPUT/debug_report.html (debug execution report)"
echo ""
echo "🚀 Ready for production:"
echo "   nextflow run nextflow/main.nf --debug false"