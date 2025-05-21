#!/bin/bash
# Run the RDKit property calculation script with different options

# Colors for output
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

print_header() {
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}========================================${NC}"
}

# Ensure we have the right environment
if [ ! -f "rdkit_wrapper.py" ]; then
    echo -e "${RED}Error: rdkit_wrapper.py not found!${NC}"
    echo "Please make sure you're in the correct directory."
    exit 1
fi

if [ ! -f "populate_rdkit_molecular_properties.py" ]; then
    echo -e "${RED}Error: populate_rdkit_molecular_properties.py not found!${NC}"
    echo "Please make sure you're in the correct directory."
    exit 1
fi

# Make the script executable
chmod +x populate_rdkit_molecular_properties.py

# Function to run with options
run_calculation() {
    local where_clause="$1"
    local description="$2"
    local batch_size="${3:-100}"
    
    print_header "Running RDKit Property Calculation: $description"
    
    if [ -n "$where_clause" ]; then
        echo -e "${YELLOW}WHERE clause: $where_clause${NC}"
        python3 populate_rdkit_molecular_properties.py --where "$where_clause" --batch-size "$batch_size"
    else
        python3 populate_rdkit_molecular_properties.py --batch-size "$batch_size"
    fi
    
    echo
}

# Check command line arguments
if [ "$1" == "--help" ] || [ "$1" == "-h" ]; then
    echo "Usage: $0 [option]"
    echo
    echo "Options:"
    echo "  --known        Process only known cryoprotectants"
    echo "  --all          Process all molecules in database (may take a long time)"
    echo "  --recent       Process recently added molecules (last 7 days)"
    echo "  --container    Run inside the RDKit container"
    echo "  --sample N     Process a random sample of N molecules"
    echo "  --help, -h     Show this help message"
    exit 0
fi

# Default action if no arguments
if [ $# -eq 0 ]; then
    print_header "CryoProtect RDKit Property Calculator"
    echo "This script calculates molecular properties using RDKit for molecules in the database."
    echo
    echo "Please specify an option:"
    echo "  ./run_rdkit_property_calculation.sh --known    (Process only known cryoprotectants)"
    echo "  ./run_rdkit_property_calculation.sh --all      (Process all molecules - may take time)"
    echo "  ./run_rdkit_property_calculation.sh --recent   (Process recently added molecules)"
    echo "  ./run_rdkit_property_calculation.sh --container (Run inside RDKit container)"
    echo "  ./run_rdkit_property_calculation.sh --sample 100 (Process 100 random molecules)"
    echo
    echo "For more options: ./run_rdkit_property_calculation.sh --help"
    exit 0
fi

# Process options
case "$1" in
    --known)
        # Search for names matching known cryoprotectants
        known_names="name ILIKE '%glycerol%' OR name ILIKE '%dmso%' OR name ILIKE '%dimethyl sulfoxide%' OR "\
"name ILIKE '%ethylene glycol%' OR name ILIKE '%propylene glycol%' OR "\
"name ILIKE '%trehalose%' OR name ILIKE '%sucrose%' OR name ILIKE '%mannitol%' OR "\
"name ILIKE '%dextran%' OR name ILIKE '%ficoll%' OR name ILIKE '%polyvinylpyrrolidone%' OR "\
"name ILIKE '%hydroxyethyl starch%' OR name ILIKE '%formamide%' OR name ILIKE '%methanol%' OR "\
"name ILIKE '%glucose%' OR name ILIKE '%sorbitol%'"
        
        run_calculation "$known_names" "Known Cryoprotectants" 20
        ;;
    --all)
        echo -e "${YELLOW}Warning: This will process all molecules in the database and may take a long time.${NC}"
        read -p "Proceed? (y/n) " answer
        if [[ $answer == "y" ]]; then
            run_calculation "" "All Molecules" 200
        else
            echo "Operation cancelled."
        fi
        ;;
    --recent)
        # Process molecules added in the last 7 days
        run_calculation "created_at > NOW() - INTERVAL '7 days'" "Recently Added Molecules" 100
        ;;
    --container)
        print_header "Running in RDKit Container"
        if command -v podman &> /dev/null; then
            echo "Using podman to run the container..."
            # Make sure our files are accessible to the container
            chmod +x populate_rdkit_molecular_properties.py
            chmod +x rdkit_wrapper.py
            chmod +x mock_rdkit_formula.py
            
            # Create a temporary directory for mounting
            TEMP_DIR=$(mktemp -d)
            cp populate_rdkit_molecular_properties.py rdkit_wrapper.py mock_rdkit_formula.py "$TEMP_DIR"
            
            # Run in the container
            podman run -it --rm -v "$TEMP_DIR:/app" localhost/cryoprotect-rdkit python3 /app/populate_rdkit_molecular_properties.py --where "name ILIKE '%glycerol%'" --batch-size 10
            
            # Clean up
            rm -rf "$TEMP_DIR"
        else
            echo -e "${RED}Error: podman not found. Please install podman or use a native Python environment.${NC}"
            exit 1
        fi
        ;;
    --sample)
        # Process a random sample of molecules
        if [ -z "$2" ]; then
            echo -e "${RED}Error: Missing sample size${NC}"
            echo "Usage: $0 --sample N"
            exit 1
        fi
        sample_size="$2"
        run_calculation "random() < 0.1 LIMIT $sample_size" "Random Sample ($sample_size molecules)" 50
        ;;
    *)
        echo -e "${RED}Error: Unknown option $1${NC}"
        echo "Use --help to see available options."
        exit 1
        ;;
esac

print_header "Property Calculation Complete"
echo -e "${GREEN}Results are stored in the molecular_properties table.${NC}"
echo "For scientific interpretation of these properties, see CRYOPROTECTANT_PROPERTY_GUIDE.md"