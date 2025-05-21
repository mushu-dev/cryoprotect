#!/bin/bash
# Run the enhanced property verification script and output results

echo "Running enhanced property verification..."
python3 enhanced_property_verification.py --output property_verification_report.json

if [ $? -eq 0 ]; then
    echo "Verification completed successfully."
    echo "Detailed report saved to property_verification_report.json"
    
    # Display a summary of the verification results
    echo ""
    echo "Summary of verification results:"
    echo "================================"
    cat property_verification_report.json | grep -A 10 "\"overall_accuracy\":"
    
    # Check if we should also try running with RDKit
    if command -v rdkit > /dev/null 2>&1; then
        echo ""
        echo "RDKit is available. Comprehensive verification was performed."
    else
        echo ""
        echo "RDKit is not available in this environment."
        echo "Only reference value comparison was performed."
        echo "For more comprehensive verification, run this in the RDKit container."
    fi
else
    echo "Verification failed. Check logs for details."
fi