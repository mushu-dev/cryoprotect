# RDKit Cryoprotectant Property Calculator Guide

This guide explains how to use the RDKit-based tools to calculate and store molecular properties relevant to cryoprotectant effectiveness.

## Overview

The CryoProtect project includes an RDKit-based tool for calculating and storing cryoprotectant-specific molecular properties:

**`rdkit_property_calculator.py`** - Uses the Supabase client library with service role authentication

This tool uses the `rdkit_wrapper.py` module, which provides a unified interface to RDKit with fallback to a mock implementation when RDKit is not available.

## Cryoprotectant Properties

These tools calculate the following cryoprotectant-specific properties based on standard molecular descriptors:

| Property | Description | Scientific Significance |
|----------|-------------|-------------------------|
| H-Bond Donor/Acceptor Ratio | Ratio of hydrogen bond donors to acceptors | Ideal cryoprotectants have a balanced ratio (0.5-2.0) for optimal ice-binding and water replacement |
| Total H-Bonding Capacity | Sum of H-bond donors and acceptors | Higher values generally correlate with better cryoprotection due to increased water interaction |
| Polarity Index | Ratio of TPSA to molecular surface area | Indicates the molecule's ability to interact with water and displace it in cellular structures |
| Membrane Interaction Score | Estimated ability to interact with cell membranes | Critical for penetrating cells and providing intracellular protection |
| Ice Interaction Potential | Estimated ability to disrupt ice crystal formation | Direct measure of a molecule's ability to prevent damaging ice crystallization |
| Vitrification Potential | Estimated ability to promote vitrification | Indicates how well a compound can promote the formation of an amorphous glassy state instead of ice |
| Estimated Toxicity | Simplified toxicity estimation | Lower values are better; based on molecular properties associated with toxicity |
| Cryoprotectant Score | Overall effectiveness score (0-10) | Weighted combination of all properties, with higher values indicating better predicted performance |

The cryoprotectant score calculation is a weighted combination of these properties, designed to identify molecules with the best overall profile for cryoprotection.

## Installation Requirements

### Core Dependencies
- Python 3.6 or higher
- RDKit (recommended for accurate calculations)
- Supabase Python client (`supabase-py`)
- python-dotenv (for loading environment variables)
- tqdm (for progress bars)

### Installation with Conda (Recommended)
```bash
conda create -n cryoprotect python=3.9
conda activate cryoprotect
conda install -c conda-forge rdkit
pip install supabase python-dotenv tqdm
```

### Without RDKit
If RDKit cannot be installed, the tools will use a mock implementation with reduced accuracy, but all functionality will still work.

## Authentication Setup

The `rdkit_property_calculator.py` script requires Supabase authentication credentials:

1. Create a `.env` file in the project root (based on `.env-example`):
   ```
   SUPABASE_URL=your-project-id.supabase.co
   SUPABASE_SERVICE_ROLE_KEY=your-service-role-key
   ```

2. If no `.env` file is found, you will be prompted for credentials.

3. You can also pass credentials via command-line arguments:
   ```bash
   ./run_rdkit_calculator.sh --url your-project-id.supabase.co --key your-service-role-key
   ```

## Usage

### Running the Property Calculator

#### Using the Wrapper Script

The easiest way to run the tool is using the wrapper script:

```bash
./run_rdkit_calculator.sh [options]
```

#### Command-line Options

The calculator supports these options:

- `--known`: Process only known cryoprotectants
- `--sample N`: Process a random sample of N molecules
- `--limit N`: Limit processing to N molecules
- `--url URL`: Supabase URL
- `--key KEY`: Supabase service role key
- `--help`: Show help message

### Example Uses

```bash
# Calculate properties for known cryoprotectants
./run_rdkit_calculator.sh --known

# Process a random sample of 100 molecules
./run_rdkit_calculator.sh --sample 100

# Process first 50 molecules in the database
./run_rdkit_calculator.sh --limit 50

# Specify Supabase credentials on the command line
./run_rdkit_calculator.sh --known --url your-project-id.supabase.co --key your-service-role-key
```

## Testing Property Calculations

To verify that the property calculations work correctly without connecting to the database:

```bash
./test_rdkit_cryoprotectant_properties.py
```

This will calculate properties for known cryoprotectants and save the results to `test_results/cryoprotectant_property_calculations.json`.

## Troubleshooting

### Authentication Issues

If you encounter authentication errors:

1. **Check your credentials**: Ensure your Supabase URL and service role key are correct
2. **RLS Policies**: The service role key is used to bypass Row-Level Security (RLS)
3. **API Limits**: Be aware of any API rate limits on your Supabase project

### RDKit Issues

If RDKit is not available:

1. The tools will use a mock implementation with reduced accuracy
2. Install RDKit: `conda install -c conda-forge rdkit`
3. Verify RDKit is installed: `python -c "import rdkit; print(rdkit.__version__)"`

### Database Issues

If database updates fail:

1. **Check table structure**: Ensure the required tables exist (property_types, calculation_methods, etc.)
2. **Check permissions**: The service role must have write access to the tables
3. **Review error messages**: Look for specific error messages in the logs

## Customizing Property Calculations

To customize the property calculations:

1. Edit the `CRYOPROTECTANT_PROPERTIES` list in either script
2. Modify the `calculate_cryoprotectant_score` function to change the weighting
3. Add new property types to the database if needed

## Scientific Background

### Key Factors for Effective Cryoprotectants

1. **Water interaction**: Ability to form hydrogen bonds with water molecules
2. **Ice disruption**: Ability to prevent or disrupt ice crystal formation
3. **Membrane penetration**: Ability to cross cell membranes for intracellular protection
4. **Vitrification promotion**: Ability to promote formation of a glassy state
5. **Low toxicity**: Minimal toxicity to cells at effective concentrations

### Property Calculation Formulas

The property calculations are based on these core molecular descriptors:
- Molecular weight
- LogP (octanol-water partition coefficient)
- TPSA (topological polar surface area)
- H-bond donors and acceptors
- Ring count
- Rotatable bonds

These are combined in scientifically informed formulas to predict cryoprotectant effectiveness.

## References

1. Fuller, B. J. (2004). Cryoprotectants: the essential antifreezes to protect life in the frozen state. CryoLetters, 25(6), 375-388.
2. Best, B. P. (2015). Cryoprotectant toxicity: facts, issues, and questions. Rejuvenation research, 18(5), 422-436.
3. Hubálek, Z. (2003). Protectants used in the cryopreservation of microorganisms. Cryobiology, 46(3), 205-229.
4. Elliott, G. D., Wang, S., & Fuller, B. J. (2017). Cryoprotectants: A review of the actions and applications of cryoprotective solutes that modulate cell recovery from ultra-low temperatures. Cryobiology, 76, 74-91.