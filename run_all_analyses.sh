#!/bin/bash

# ================================================================
# Brain Microprotein Discovery - Summary Analysis Pipeline
# ================================================================
# This script runs ONLY the final summary analysis steps to generate
# key results from the processed master data file.
#
# ⚠️ DATA NOTICE:
# - Raw input data are NOT included (privacy + size constraints)
# - Analyses run on a pre-processed master file
#
# PURPOSE:
# Demonstrate methodology and reproduce key findings without
# requiring raw sequencing or proteomics data.
#
# MODES:
# 1. docs (default) – show analysis steps and methodology
# 2. run            – execute analyses and generate result files
#
# Usage:
#   bash run_all_analyses.sh --mode=docs   # Documentation mode
#   bash run_all_analyses.sh --mode=run    # Generate summary results
# ================================================================

set -e  # Exit on any error

# Parse command line arguments
MODE="docs"  # Default to documentation mode
while [[ $# -gt 0 ]]; do
  case $1 in
    --mode=*)
      MODE="${1#*=}"
      shift
      ;;
    *)
      echo "Unknown option $1"
      echo "Usage: $0 [--mode=docs|run]"
      exit 1
      ;;
  esac
done

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print colored messages
print_step() {
    echo -e "${BLUE}[STEP]${NC} $1"
}

print_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_docs() {
    echo -e "${BLUE}[DOCS]${NC} $1"
}

# Function to run command or show documentation
run_or_show() {
    local cmd="$1"
    local description="$2"
    local file_to_check="$3"
    
    if [ "$MODE" = "docs" ]; then
        print_docs "$description"
        echo "    Command: $cmd"
        if [ -n "$file_to_check" ] && [ ! -f "$file_to_check" ]; then
            echo "    Status: Script file not found (expected for documentation mode)"
        else
            echo "    Status: Script file available"
        fi
        echo ""
    else
        if [ -n "$file_to_check" ] && [ -f "$file_to_check" ]; then
            print_step "$description"
            eval "$cmd"
            print_success "$description completed"
        else
            print_warning "$description - file not found: $file_to_check"
        fi
    fi
}

# Get the absolute path to the repository root
REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
CODE_DIR="${REPO_ROOT}/Code"
RESULTS_DIR="${REPO_ROOT}/Results"

echo "================================================================"
echo "Brain Microprotein Discovery Project - Results Summary Pipeline"
echo "================================================================"
echo "Repository: ${REPO_ROOT}"
echo "Code Directory: ${CODE_DIR}"
echo "Results Directory: ${RESULTS_DIR}"
echo "Mode: ${MODE}"
echo "================================================================"

if [ "$MODE" = "docs" ]; then
    echo ""
    echo "📖 DOCUMENTATION MODE"
    echo "   This mode shows the analysis methodology and summary steps."
    echo "   Processing scripts are documented but not executed (require raw data)."
    echo "   To generate summary results: bash run_all_analyses.sh --mode=run"
    echo ""
elif [ "$MODE" = "run" ]; then
    echo ""
    echo "🚀 SUMMARY RESULTS MODE"
    echo "   This mode runs summary analysis scripts to generate key findings."
    echo "   Uses our processed master data file (microprotein_master.csv)."
    echo ""
else
    print_error "Invalid mode: $MODE. Use 'docs' or 'run'"
    exit 1
fi

# Check if required directories exist
if [ ! -d "$CODE_DIR" ]; then
    print_error "Code directory not found: $CODE_DIR"
    exit 1
fi

if [ ! -d "$RESULTS_DIR" ]; then
    if [ "$MODE" = "run" ]; then
        print_warning "Results directory not found, creating: $RESULTS_DIR"
        mkdir -p "$RESULTS_DIR"
    else
        print_docs "Results directory would be created at: $RESULTS_DIR"
    fi
fi

# Change to repository root
cd "$REPO_ROOT"

# ================================================================
# DATA FILE PREPARATION
# ================================================================
# Check for master data file and unzip if necessary

print_step "Checking for master data file..."

MASTER_CSV="${CODE_DIR}/data/microprotein_master.csv"
MASTER_ZIP="${CODE_DIR}/data/micrprotein_master.zip"

if [ -f "$MASTER_CSV" ]; then
    # CSV file already exists (previously unzipped)
    if [ "$MODE" = "docs" ]; then
        print_docs "Master data file available: $MASTER_CSV"
    else
        print_success "Master data file found: $MASTER_CSV"
    fi
elif [ -f "$MASTER_ZIP" ]; then
    # ZIP file exists, need to unzip it
    if [ "$MODE" = "docs" ]; then
        print_docs "Master data file would be extracted from: $MASTER_ZIP"
        echo "    Command: unzip -o \"$MASTER_ZIP\" -d \"${CODE_DIR}/data/\""
        echo "    Status: ZIP file available for extraction"
        echo ""
    else
        print_step "Extracting master data file from ZIP archive..."
        unzip -o "$MASTER_ZIP" -d "${CODE_DIR}/data/"
        if [ -f "$MASTER_CSV" ]; then
            print_success "Master data file extracted successfully: $MASTER_CSV"
        else
            print_error "Failed to extract master data file from ZIP"
            exit 1
        fi
    fi
else
    # Neither CSV nor ZIP exists
    print_error "Master data file not found: $MASTER_CSV"
    print_error "ZIP archive not found: $MASTER_ZIP"
    echo "Please ensure the master data ZIP file (micrprotein_master.zip) is available in Code/data/"
    exit 1
fi

# ================================================================
# ANALYSIS RESULTS SUMMARY
# ================================================================
# This script runs only the summary analysis scripts that generate
# final results from the master data file. These scripts demonstrate
# the key findings and methodology of our brain microprotein discovery study.

print_step "Phase 1: Annotation and Discovery Results"

# 1.1 smORF Annotation Analysis Results
print_step "1.1 Running smORF annotation analysis summaries..."
cd "${CODE_DIR}/Microprotein_annotation_summary"

if [ "$MODE" = "docs" ]; then
    print_docs "smORF Annotator pipeline (data processing - not run in summary mode)"
    echo "    Command: bash run_Annotator.sh"
    echo "    Status: Processing script available (requires raw data)"
    echo ""
fi

run_or_show "python Brain_Microproteins_Discovery_summary.py" "Brain Microproteins Discovery summary" "Brain_Microproteins_Discovery_summary.py"
run_or_show "python ShortStop_Microproteins_summary.py" "ShortStop Microproteins summary" "ShortStop_Microproteins_summary.py"
run_or_show "python Create_BED_GTF_FASTA_files.py" "Create BED, GTF, and FASTA coordinates for genomic visualization" "Create_BED_GTF_FASTA_files.py"

# 1.2 RP3 Analysis Results
print_step "1.2 Running RP3 analysis summary..."
cd "${CODE_DIR}/RP3_analysis"

run_or_show "python RP3_Results_summary.py" "RP3 Results summary" "RP3_Results_summary.py"

# ================================================================
# PHASE 2: Transcriptomics Results Summary
# ================================================================

print_step "Phase 2: Transcriptomics Results Summary"

# 2.1 Long-read RNA Analysis Results
print_step "2.1 Running long-read RNA analysis summary..."
cd "${CODE_DIR}/Longread_RNA_analysis"

# Document ESPRESSO data processing scripts without running them
if [ "$MODE" = "docs" ]; then
    print_docs "ESPRESSO data processing (processing scripts - not run in summary mode)"
    echo "    Scripts: convert_ESPERSSO_to_CPM.sh, convert_ESPRESSO_to_CPM_and_filter.py, deseq_brain_espresso.r"
    echo "    Status: Processing scripts available (require raw sequencing data)"
    echo ""
fi

run_or_show "python Long-Read_Transcriptomics_Results_summary.py" "Long-Read Transcriptomics Results summary" "Long-Read_Transcriptomics_Results_summary.py"

# 2.2 Short-read RNA Analysis Results
print_step "2.2 Running short-read RNA analysis summary..."
cd "${CODE_DIR}/Shortread_RNA_analysis"

# Document shortread DESeq processing scripts without running them
if [ "$MODE" = "docs" ]; then
    print_docs "Short-read DESeq processing (processing scripts - not run in summary mode)"
    echo "    Scripts: RNA_differential_expression.R"
    echo "    Status: Processing scripts available (require raw sequencing data)"
    echo ""
fi

run_or_show "python Short-Read_Transcriptomics_Results_summary.py" "Short-Read Transcriptomics Results summary" "Short-Read_Transcriptomics_Results_summary.py"

# ================================================================
# PHASE 3: Proteomics Results Summary
# ================================================================

print_step "Phase 3: Proteomics Results Summary"

cd "${CODE_DIR}/Peptide_TMT_analysis"

# 3.1 TMT Data Processing Documentation
print_step "3.1 TMT data processing overview..."

# Document fragpipe results processing scripts without running them
if [ "$MODE" = "docs" ]; then
    print_docs "FragPipe results processing (processing scripts - not run in summary mode)"
    echo "    Processing scripts: process_proteinID_from_TMT_round1.py, process_proteinID_frmo_TMT_round2.py, find_unique_tryptic_peptides.py"
    echo "    Matrix generation: generate_raw_TMT_intensity_matrix_across_batches.R, TMT_regressed_corrrected_matrix.R"
    echo "    TAMPOR analysis: TAMPOR_round1.R, TAMPOR_round2.R, TAMPOR_combined_rounds.R, TAMPOR_microproteins.R"
    echo "    Statistical tests: TMT_ANOVA.R"
    echo "    Status: Processing scripts available (require raw proteomics data)"
    echo ""
fi

# 3.2 Proteomics Summary
run_or_show "python Proteomics_Results_summary.py" "Proteomics Results summary" "Proteomics_Results_summary.py"

# ================================================================
# PHASE 4: Single-cell RNA-seq Results Summary
# ================================================================

print_step "Phase 4: Single-cell RNA-seq Results Summary"

cd "${CODE_DIR}/scRNAseq_summary_merging_analysis"

if [ "$MODE" = "docs" ]; then
    print_docs "scRNA-seq summary analysis (processing script - not run in summary mode)"
    echo "    Command: Rscript scRNAseq_summary.R"
    echo "    Status: Processing script available (requires scRNA-seq data)"
    echo ""
fi

# Run the scRNA-seq summary script
run_or_show "Rscript scRNAseq_summary.R" "scRNA-seq Summary analysis" "scRNAseq_summary.R"

# ================================================================
# PHASE 5: Master Supplemental Tables
# ================================================================

print_step "Phase 5: Master Supplemental Tables"

cd "${RESULTS_DIR}"

if [ "$MODE" = "docs" ]; then
    print_docs "Master supplemental tables workbook generation"
    echo "    Command: python generate_supplemental_tables.py --include-scrna"
    echo "    Status: Script file available"
    echo ""
fi

# Generate supplemental workbook and always include scRNA-seq table
run_or_show "python generate_supplemental_tables.py --include-scrna" "Generate Supplemental Tables workbook (including scRNA-seq table)" "generate_supplemental_tables.py"

# ================================================================
# FINAL STEPS
# ================================================================

print_step "Final Steps: Summary Results Generated"

# Return to repository root
cd "$REPO_ROOT"

# Check if key output files were generated
echo ""
echo "================================================================"
echo "Summary Analysis Pipeline Completed!"
echo "================================================================"
echo ""

if [ "$MODE" = "docs" ]; then
    print_docs "Expected summary output files:"
    echo "    Results/Annotations/Brain_Microproteins_Discovery_summary.csv"
    echo "    Results/Annotations/ShortStop_Microproteins_summary.csv"
    echo "    GTF_and_BED_files/Unreviewed_Brain_Microproteins_CDS_absent_from_UniProt.bed"
    echo "    Results/RP3/RP3_Results_summary.csv"
    echo "    Results/Proteomics/Proteomics_Results_summary.csv"
    echo "    Results/Transcriptomics/Long-Read_Transcriptomics_Results_summary.csv"
    echo "    Results/Transcriptomics/Short-Read_Transcriptomics_Results_summary.csv"
    echo "    Results/scRNA_Enrichment/scRNA_Enrichment_summary.csv"
    echo "    ../supplementary/Supplemental_Tables.xlsx"
    echo ""
    echo "================================================================"
    echo "Documentation Mode Summary"
    echo "================================================================"
    echo ""
    print_success "Methodology documentation complete!"
    echo ""
    echo "This pipeline demonstrates our brain microprotein discovery approach:"
    echo "  1. smORF annotation and classification (from assembled transcriptomes)"
    echo "  2. Ribosome profiling validation (RP3 pipeline)"
    echo "  3. Mass spectrometry evidence (TMT proteomics)"
    echo "  4. Transcriptomic context (RNA-seq differential expression)"
    echo "  5. Single-cell expression patterns (scRNA-seq enrichment)"
    echo ""
    echo "To generate summary results:"
    echo "  bash run_all_analyses.sh --mode=run"
    echo ""
else
    # Check which output files exist
    TOTAL_OUTPUTS=9
    FOUND_OUTPUTS=0
    
    echo "Checking for generated output files:"
    
    OUTPUT_FILES=(
        "Results/Annotations/Brain_Microproteins_Discovery_summary.csv"
        "Results/Annotations/ShortStop_Microproteins_summary.csv"
        "GTF_and_BED_files/Unreviewed_Brain_Microproteins_CDS_absent_from_UniProt.bed"
        "Results/RP3/RP3_Results_summary.csv"
        "Results/Proteomics/Proteomics_Results_summary.csv"
        "Results/Transcriptomics/Long-Read_Transcriptomics_Results_summary.csv"
        "Results/Transcriptomics/Short-Read_Transcriptomics_Results_summary.csv"
        "Results/scRNA_Enrichment/scRNA_Enrichment_summary.csv"
        "../supplementary/Supplemental_Tables.xlsx"
    )
    
    for file in "${OUTPUT_FILES[@]}"; do
        if [ -f "$file" ]; then
            print_success "Found: $file"
            ((FOUND_OUTPUTS++))
        else
            print_warning "Missing: $file"
        fi
    done
    
    echo ""
    echo "================================================================"
    echo "Summary: $FOUND_OUTPUTS/$TOTAL_OUTPUTS key output files found"
    echo "================================================================"
    
    if [ $FOUND_OUTPUTS -eq $TOTAL_OUTPUTS ]; then
        print_success "All summary analyses completed successfully!"
        echo ""
        echo "You can now explore the summary results in the Results/ directory:"
        echo "  - Results/Annotations/     : smORF annotations and classifications"
        echo "  - Results/Proteomics/      : TMT-MS quantification results"
        echo "  - Results/RP3/             : RiboCode ribosome profiling results"
        echo "  - Results/Transcriptomics/ : RNA-seq differential expression summaries"
        echo ""
        echo "These files demonstrate the key findings from our brain microprotein discovery study."
        echo ""
    else
        print_warning "Some summary analyses may have failed. Check the output above for details."
        echo ""
        echo "Common issues:"
        echo "  - Missing master data file (check Code/data/microprotein_master.csv)"
        echo "  - Missing Python dependencies (pandas, numpy)"
        echo "  - File path issues (should use relative paths)"
        echo ""
    fi
fi

echo "For questions or issues, contact: brmiller@salk.edu"
echo "================================================================"
