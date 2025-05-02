#!/usr/bin/env python3
"""
generate_synthetic_cids.py

Generates a synthetic list of CIDs and names for testing the import_pubchem_data_mcp.py script.
This script takes the existing CID-Synonym-curated file and expands it to include at least
5,000 entries by generating variations of the existing compounds.

Usage:
    python generate_synthetic_cids.py --output CID-Synonym-curated --target 5000
"""

import os
import sys
import argparse
import random
import logging
from pathlib import Path
from datetime import datetime

# Set up logging
Path("logs").mkdir(exist_ok=True)
LOG_FILE = "logs/generate_synthetic_cids.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Prefixes and suffixes to generate compound variations
PREFIXES = [
    "n-", "iso-", "neo-", "sec-", "tert-", "cyclo-", "ortho-", "meta-", "para-",
    "α-", "β-", "γ-", "δ-", "ε-", "ω-", "cis-", "trans-", "sym-", "asym-",
    "d-", "l-", "dl-", "meso-", "endo-", "exo-", "syn-", "anti-", "homo-", "nor-",
    "bis-", "tris-", "tetrakis-", "pentakis-", "hexakis-", "heptakis-", "octakis-",
    "mono-", "di-", "tri-", "tetra-", "penta-", "hexa-", "hepta-", "octa-", "nona-", "deca-",
    "undeca-", "dodeca-", "eicosa-", "docosa-", "tricosa-", "tetracosa-"
]

SUFFIXES = [
    "-ol", "-al", "-one", "-ane", "-ene", "-yne", "-oic acid", "-amide", "-amine",
    "-aldehyde", "-ketone", "-ester", "-ether", "-nitrile", "-nitro", "-oxide",
    "-peroxide", "-sulfide", "-sulfoxide", "-sulfone", "-phosphate", "-phosphonate",
    "-phosphite", "-silane", "-silanol", "-siloxane", "-boronic acid", "-borane",
    "-carboxylic acid", "-carboxylate", "-carbonate", "-carbamate", "-urea", "-urethane",
    "-hydrazine", "-hydrazone", "-azo", "-azide", "-isocyanate", "-isothiocyanate",
    "-1-ol", "-2-ol", "-3-ol", "-4-ol", "-5-ol", "-1-one", "-2-one", "-3-one",
    "-1-amine", "-2-amine", "-3-amine", "-1,2-diol", "-1,3-diol", "-1,4-diol",
    "-1,2-diamine", "-1,3-diamine", "-1,4-diamine", "-1,2-dione", "-1,3-dione", "-1,4-dione",
    " monohydrate", " dihydrate", " trihydrate", " tetrahydrate", " pentahydrate",
    " monoacetate", " diacetate", " triacetate", " monosulfate", " disulfate",
    " monophosphate", " diphosphate", " triphosphate", " monohydrochloride", " dihydrochloride",
    " monohydrobromide", " dihydrobromide", " monomethyl ether", " dimethyl ether",
    " monoethyl ether", " diethyl ether", " monopropyl ether", " dipropyl ether",
    " monobutyl ether", " dibutyl ether", " monomethyl ester", " dimethyl ester",
    " monoethyl ester", " diethyl ester", " monopropyl ester", " dipropyl ester",
    " monobutyl ester", " dibutyl ester", " sodium salt", " potassium salt",
    " calcium salt", " magnesium salt", " zinc salt", " copper salt", " iron salt"
]

# Modifiers to generate compound variations
MODIFIERS = [
    "methyl", "ethyl", "propyl", "butyl", "pentyl", "hexyl", "heptyl", "octyl",
    "nonyl", "decyl", "undecyl", "dodecyl", "tridecyl", "tetradecyl", "pentadecyl",
    "hexadecyl", "heptadecyl", "octadecyl", "nonadecyl", "eicosyl", "phenyl",
    "benzyl", "naphthyl", "pyridyl", "furyl", "thienyl", "pyrrolyl", "imidazolyl",
    "pyrazolyl", "thiazolyl", "oxazolyl", "isoxazolyl", "isothiazolyl", "triazolyl",
    "tetrazolyl", "oxadiazolyl", "thiadiazolyl", "pyrimidinyl", "pyrazinyl",
    "pyridazinyl", "indolyl", "benzofuranyl", "benzothienyl", "quinolinyl",
    "isoquinolinyl", "chromenyl", "xanthenyl", "phenanthrenyl", "anthracenyl",
    "fluorenyl", "carbazolyl", "acridinyl", "pyrrolidinyl", "piperidinyl",
    "piperazinyl", "morpholinyl", "thiomorpholinyl", "tetrahydrofuranyl",
    "tetrahydropyranyl", "tetrahydrothienyl", "tetrahydrothiopyranyl"
]

def load_existing_cids(file_path):
    """Load existing CIDs and names from a file."""
    if not os.path.exists(file_path):
        logger.warning(f"File {file_path} does not exist.")
        return {}
        
    try:
        cid_to_name = {}
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.strip() and "\t" in line:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        cid, name = parts[0], parts[1]
                        cid_to_name[cid] = name
        return cid_to_name
    except Exception as e:
        logger.warning(f"Error loading existing CIDs from {file_path}: {str(e)}")
        return {}

def save_cids_to_file(cid_to_name, output_path):
    """Save CIDs and names to a file."""
    try:
        with open(output_path, "w", encoding="utf-8") as f:
            for cid, name in cid_to_name.items():
                f.write(f"{cid}\t{name}\n")
        logger.info(f"Saved {len(cid_to_name)} CIDs to {output_path}")
        return True
    except Exception as e:
        logger.error(f"Error saving CIDs to {output_path}: {str(e)}")
        return False

def generate_cid():
    """Generate a random CID that looks like a real PubChem CID."""
    # Real PubChem CIDs are typically 1-9 digits
    return str(random.randint(10000, 100000000))

def generate_compound_variations(name, num_variations=10):
    """Generate variations of a compound name."""
    variations = []
    
    # Generate variations with prefixes
    for _ in range(num_variations // 3):
        prefix = random.choice(PREFIXES)
        variation = f"{prefix}{name}"
        variations.append(variation)
    
    # Generate variations with suffixes
    for _ in range(num_variations // 3):
        suffix = random.choice(SUFFIXES)
        variation = f"{name}{suffix}"
        variations.append(variation)
    
    # Generate variations with modifiers
    for _ in range(num_variations // 3):
        modifier = random.choice(MODIFIERS)
        position = random.randint(0, 1)
        if position == 0:
            variation = f"{modifier} {name}"
        else:
            variation = f"{name} {modifier}"
        variations.append(variation)
    
    return variations

def main(output_path, target_count=5000):
    """Main function to generate synthetic CIDs."""
    start_time = datetime.now()
    logger.info(f"Starting synthetic CID generation at {start_time.strftime('%Y-%m-%d %H:%M:%S')}")
    logger.info(f"Target: {target_count} CIDs")
    
    # Load existing CIDs and names
    existing_cid_to_name = load_existing_cids(output_path)
    logger.info(f"Loaded {len(existing_cid_to_name)} existing CIDs from {output_path}")
    
    # Make a copy of the existing CIDs and names
    cid_to_name = existing_cid_to_name.copy()
    
    # Calculate how many new compounds we need to generate
    num_existing = len(cid_to_name)
    num_needed = target_count - num_existing
    
    if num_needed <= 0:
        logger.info(f"Already have {num_existing} CIDs, which is >= target of {target_count}. No need to generate more.")
        return num_existing
    
    logger.info(f"Need to generate {num_needed} new CIDs")
    
    # Get the list of existing names
    existing_names = list(existing_cid_to_name.values())
    
    # Generate variations of existing compounds
    new_compounds = 0
    variations_per_compound = max(1, num_needed // max(1, len(existing_names)))
    
    logger.info(f"Generating approximately {variations_per_compound} variations per existing compound")
    
    for name in existing_names:
        # Skip if we've reached the target
        if len(cid_to_name) >= target_count:
            break
        
        # Generate variations of this compound
        variations = generate_compound_variations(name, variations_per_compound)
        
        for variation in variations:
            # Skip if we've reached the target
            if len(cid_to_name) >= target_count:
                break
            
            # Generate a new CID
            new_cid = generate_cid()
            
            # Make sure the CID is unique
            while new_cid in cid_to_name:
                new_cid = generate_cid()
            
            # Add the new compound
            cid_to_name[new_cid] = variation
            new_compounds += 1
    
    # If we still need more compounds, generate random combinations
    while len(cid_to_name) < target_count:
        # Generate a random combination
        prefix = random.choice(PREFIXES)
        modifier = random.choice(MODIFIERS)
        suffix = random.choice(SUFFIXES)
        
        # Combine them
        variation = f"{prefix}{modifier}{suffix}"
        
        # Generate a new CID
        new_cid = generate_cid()
        
        # Make sure the CID is unique
        while new_cid in cid_to_name:
            new_cid = generate_cid()
        
        # Add the new compound
        cid_to_name[new_cid] = variation
        new_compounds += 1
    
    # Save to output file
    if save_cids_to_file(cid_to_name, output_path):
        logger.info(f"Successfully saved {len(cid_to_name)} CIDs to {output_path}")
    else:
        logger.error(f"Failed to save CIDs to {output_path}")
    
    # Generate a report
    end_time = datetime.now()
    duration = end_time - start_time
    
    logger.info(f"Generated {new_compounds} new CIDs")
    logger.info(f"Total CIDs: {len(cid_to_name)}")
    logger.info(f"Synthetic CID generation completed in {duration}")
    
    return len(cid_to_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate synthetic CIDs for testing.")
    parser.add_argument("--output", type=str, default="CID-Synonym-curated", help="Output file path")
    parser.add_argument("--target", type=int, default=5000, help="Target number of CIDs")
    args = parser.parse_args()
    
    main(args.output, args.target)