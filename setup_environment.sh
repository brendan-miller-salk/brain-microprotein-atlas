#!/bin/bash

# ================================================================
# Brain Microprotein Discovery Project - Environment Setup
# ================================================================
# This script sets up the required environment and dependencies
# for reproducing the brain microprotein discovery analyses.
#
# Usage: bash setup_environment.sh
# ================================================================

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

print_step() {
    echo -e "${BLUE}[SETUP]${NC} $1"
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

echo "================================================================"
echo "Brain Microprotein Discovery Project - Environment Setup"
echo "================================================================"

# Check if conda is available
if command -v conda &> /dev/null; then
    print_step "Conda found. Setting up conda environment..."
    
    # Create environment from yml file
    if [ -f "environment.yml" ]; then
        print_step "Creating conda environment from environment.yml..."
        conda env create -f environment.yml
        print_success "Conda environment 'brain_microproteins' created!"
        
        echo ""
        echo "To activate the environment, run:"
        echo "  conda activate brain_microproteins"
        echo ""
    else
        print_warning "environment.yml not found. Creating basic environment..."
        conda create -n brain_microproteins python=3.10 -y
        conda activate brain_microproteins
        
        # Install Python packages
        if [ -f "requirements.txt" ]; then
            print_step "Installing Python packages from requirements.txt..."
            pip install -r requirements.txt
            print_success "Python packages installed!"
        fi
        
        # Install R packages manually
        print_step "Installing R packages..."
        conda install -c conda-forge -c bioconda -c r r-base r-essentials bioconductor-deseq2 r-dplyr r-ggplot2 bedtools -y
        print_success "Basic R packages installed!"
    fi
    
elif command -v pip &> /dev/null; then
    print_step "Conda not found. Using pip for Python packages..."
    
    if [ -f "requirements.txt" ]; then
        print_step "Installing Python packages from requirements.txt..."
        pip install -r requirements.txt
        print_success "Python packages installed!"
    else
        print_error "requirements.txt not found!"
        exit 1
    fi
    
    print_warning "Please install R and required R packages manually:"
    echo "  R packages needed: DESeq2, dplyr, ggplot2, biomaRt, etc."
    echo "  Also install bedtools for annotation pipeline"
    
else
    print_error "Neither conda nor pip found! Please install Python package manager."
    exit 1
fi

# Check for bedtools
if command -v bedtools &> /dev/null; then
    print_success "bedtools found!"
else
    print_warning "bedtools not found. Please install bedtools for annotation pipeline."
    echo "  Ubuntu/Debian: sudo apt-get install bedtools"
    echo "  macOS: brew install bedtools"
    echo "  Conda: conda install -c bioconda bedtools"
fi

# Check for R
if command -v R &> /dev/null; then
    print_success "R found!"
else
    print_warning "R not found. Please install R for statistical analyses."
    echo "  https://www.r-project.org/"
fi

echo ""
echo "================================================================"
echo "Setup Summary"
echo "================================================================"
echo ""
print_step "Next steps to reproduce the analysis:"
echo ""
echo "1. Activate the environment (if using conda):"
echo "   conda activate brain_microproteins"
echo ""
echo "2. Ensure you have the required input data in Code/data/"
echo "   (microprotein_master.csv ships in the repo and is sufficient"
echo "    to run the summary pipeline and dashboard)."
echo ""
echo "3. Run the complete summary pipeline:"
echo "   bash run_all_analyses.sh --mode=run"
echo ""
echo "4. Or run individual analysis modules:"
echo "   cd Code/[analysis_module]/"
echo "   # See README.md in each module for specific instructions"
echo ""
echo "5. Launch the interactive Streamlit dashboard:"
echo "   cd Results"
echo "   pip install -r requirements.txt   # streamlit + plotly only"
echo "   bash launch_dashboard.sh          # http://localhost:8505"
echo ""
echo "Notes:"
echo "  - PROSIT spectral validation"
echo "    (Code/Peptide_TMT_analysis/prosit/prosit_pipeline.py) needs"
echo "    pyteomics, koinapy, tqdm and network access to"
echo "    koina.wilhelmlab.org. Install koinapy manually if needed:"
echo "       pip install koinapy"
echo "  - Microscopy actin pipeline needs nd2reader / scikit-image."
echo ""
echo "================================================================"
echo "For questions or issues, contact: brmiller@salk.edu"
echo "================================================================"
