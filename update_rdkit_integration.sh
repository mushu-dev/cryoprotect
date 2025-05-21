#!/bin/bash
# Script to update existing RDKit usage to the new wrapper approach

set -e

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}Updating RDKit integration in CryoProtect project${NC}"
echo "========================================================="

# Make sure wrapper and mock implementation are executable
chmod +x rdkit_wrapper.py mock_rdkit_formula.py 2>/dev/null || true

# Create a summary file to log changes
SUMMARY_FILE="rdkit_integration_update_summary.md"

cat > $SUMMARY_FILE << 'EOF'
# RDKit Integration Update Summary

This document summarizes the changes made to integrate the robust RDKit wrapper throughout the CryoProtect project.

## Files Updated

| Filename | Changes Made | Status |
|----------|--------------|--------|
EOF

# Function to find files with direct RDKit imports
find_rdkit_imports() {
    echo -e "${BLUE}Finding files with direct RDKit imports...${NC}"
    
    # Find Python files with direct RDKit imports
    FILES=$(grep -l "import rdkit\|from rdkit" --include="*.py" -r . | grep -v "rdkit_wrapper.py\|mock_rdkit\|test_rdkit")
    
    if [ -z "$FILES" ]; then
        echo "No files with direct RDKit imports found"
        return
    fi
    
    echo -e "Found $(echo "$FILES" | wc -l) files with direct RDKit imports"
    echo "$FILES"
}

# Function to make a backup of a file before modifying
backup_file() {
    local file=$1
    local backup="${file}.bak"
    
    if [ ! -f "$backup" ]; then
        cp "$file" "$backup"
        echo "Created backup: $backup"
    fi
}

# Function to update import statements in a file
update_imports() {
    local file=$1
    local status="❌ Failed"
    
    echo -e "${BLUE}Updating imports in $file...${NC}"
    backup_file "$file"
    
    # Replace direct RDKit imports with wrapper import
    if grep -q "from rdkit import" "$file" || grep -q "import rdkit\." "$file"; then
        # Add import of the wrapper at the top of imports section
        if ! grep -q "import rdkit_wrapper" "$file" && ! grep -q "from rdkit_wrapper import" "$file"; then
            sed -i '/import rdkit\|from rdkit/i # Import RDKit wrapper for unified access\nimport rdkit_wrapper' "$file"
        fi
        
        # Add comment to mark old imports
        sed -i 's/from rdkit import/# LEGACY: from rdkit import/g' "$file"
        sed -i 's/import rdkit/# LEGACY: import rdkit/g' "$file"
        
        # Update the summary file
        status="✅ Updated"
    else
        status="⚠️ No RDKit imports found"
    fi
    
    # Add to summary
    echo "| $file | Updated imports | $status |" >> $SUMMARY_FILE
}

# Function to suggest function replacements
suggest_replacements() {
    local file=$1
    
    echo -e "${BLUE}Suggesting function replacements in $file...${NC}"
    
    # Common RDKit patterns to replace
    patterns=(
        "Chem\.MolFromSmiles" "rdkit_wrapper.create_molecule_from_smiles"
        "Descriptors\.MolWt" "props = rdkit_wrapper.calculate_properties(mol); mw = props['molecular_weight']"
        "Lipinski\.NumHDonors" "props = rdkit_wrapper.calculate_properties(mol); h_donors = props['h_donors']"
        "Chem\.MolToSmiles" "rdkit_wrapper.mol_to_smiles"
        "AllChem\.GetMorganFingerprintAsBitVect" "rdkit_wrapper.generate_fingerprint"
        "DataStructs\.TanimotoSimilarity" "rdkit_wrapper.calculate_similarity"
        "Chem\.MolToInchiKey" "rdkit_wrapper.smiles_to_inchi_key"
    )
    
    # Look for patterns in the file
    SUGGESTIONS=""
    
    for ((i=0; i<${#patterns[@]}; i+=2)); do
        pattern="${patterns[$i]}"
        replacement="${patterns[$i+1]}"
        
        if grep -q "$pattern" "$file"; then
            SUGGESTIONS+="$pattern → $replacement\n"
        fi
    done
    
    if [ -n "$SUGGESTIONS" ]; then
        echo -e "${YELLOW}Suggested replacements for $file:${NC}"
        echo -e "$SUGGESTIONS"
        
        # Add suggestions to the summary
        echo -e "Suggested replacements:\n$SUGGESTIONS" | sed 's/$/  /' >> $SUMMARY_FILE
    fi
}

# Function to check file for complicated RDKit usage
check_complexity() {
    local file=$1
    
    echo -e "${BLUE}Checking complexity of RDKit usage in $file...${NC}"
    
    # Check for advanced usage that might need manual attention
    COMPLEX_USAGE=""
    
    # Check for custom descriptors or advanced functionality
    if grep -q "rdMolDescriptors" "$file" || grep -q "MolSurf" "$file" || grep -q "FragmentCatalog" "$file"; then
        COMPLEX_USAGE+="- Contains advanced RDKit descriptors\n"
    fi
    
    # Check for 3D operations
    if grep -q "AllChem\.UFFOptimizeMolecule" "$file" || grep -q "AllChem\.EmbedMolecule" "$file"; then
        COMPLEX_USAGE+="- Contains 3D molecule operations\n"
    fi
    
    # Check for fingerprint operations
    if grep -q "GetMorganFingerprintAsBitVect" "$file" || grep -q "DataStructs\.TanimotoSimilarity" "$file"; then
        COMPLEX_USAGE+="- Contains fingerprint or similarity operations\n"
    fi
    
    # Check for drawing operations
    if grep -q "Draw\." "$file" || grep -q "rdMolDraw2D" "$file"; then
        COMPLEX_USAGE+="- Contains molecule drawing/visualization\n"
    fi
    
    if [ -n "$COMPLEX_USAGE" ]; then
        echo -e "${YELLOW}Complex RDKit usage in $file needs manual attention:${NC}"
        echo -e "$COMPLEX_USAGE"
        
        # Add to summary
        echo -e "Complex usage that needs manual attention:\n$COMPLEX_USAGE" | sed 's/$/  /' >> $SUMMARY_FILE
        
        return 1
    fi
    
    return 0
}

# Function to perform full update on a file
update_file() {
    local file=$1
    
    echo -e "${BLUE}Processing $file...${NC}"
    
    # Update imports
    update_imports "$file"
    
    # Check complexity
    if check_complexity "$file"; then
        # Suggest function replacements
        suggest_replacements "$file"
    else
        echo -e "${YELLOW}This file requires manual review due to complex RDKit usage${NC}"
    fi
    
    echo ""
}

# Main execution
FOUND_FILES=$(find_rdkit_imports)
if [ -n "$FOUND_FILES" ]; then
    echo "$FOUND_FILES" | while read -r file; do
        if [ -n "$file" ]; then
            update_file "$file"
        fi
    done
else
    echo -e "${YELLOW}No files with direct RDKit imports found to update${NC}"
    echo "| No files | No updates needed | ✅ |" >> $SUMMARY_FILE
fi

# Finish summary
cat >> $SUMMARY_FILE << 'EOF'

## Next Steps

1. Install the updated files:
   - Copy `rdkit_wrapper.py` and `mock_rdkit_formula.py` to all deployment environments
   - Make sure both files are in the Python path

2. Review complex files:
   - Files marked with complex usage need manual review and updates
   - Test each file thoroughly after updating

3. Update tests:
   - Run `rdkit_integration_test.py` in each environment
   - Fix any issues found during testing

4. Update documentation:
   - Add instructions for using the RDKit wrapper
   - Update developer guidelines

## Long-term Improvements

1. Gradually refactor direct RDKit calls to use the wrapper
2. Extend the wrapper with additional functionality as needed
3. Set up CI/CD pipeline to run tests in the RDKit container
4. Consider implementing a microservice architecture for RDKit operations

For complete details on the RDKit integration approach, see the `RDKIT_INTEGRATION_GUIDE.md` file.
EOF

echo -e "${GREEN}Update summary written to $SUMMARY_FILE${NC}"
echo "Review this file for details on changes made and next steps"

echo -e "${BLUE}Setup verification tests...${NC}"
chmod +x rdkit_integration_test.py build_rdkit_container.sh run_with_rdkit_container.sh
echo -e "${GREEN}RDKit integration update completed!${NC}"
echo
echo -e "${BLUE}Next steps:${NC}"
echo "1. Review the $SUMMARY_FILE file for detailed information"
echo "2. Run tests in both host and container environments:"
echo "   - Host: ./rdkit_integration_test.py"
echo "   - Container: ./run_with_rdkit_container.sh rdkit_integration_test.py"
echo "3. Build dedicated RDKit container for production:"
echo "   - ./build_rdkit_container.sh"
echo "4. Update CI/CD pipeline to use the RDKit container"
echo "5. Test all functionality before deploying to production"