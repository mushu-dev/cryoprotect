#!/usr/bin/env python3
"""
Mock RDKit implementation for molecular formula calculations and property predictions.
To be used when RDKit is not available in the environment.
"""

import re
import logging
import math

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Define atomic weights for common elements
ATOMIC_WEIGHTS = {
    'H': 1.00794,
    'C': 12.0107,
    'N': 14.0067,
    'O': 15.9994,
    'F': 18.9984,
    'P': 30.9738,
    'S': 32.065,
    'Cl': 35.453,
    'Br': 79.904,
    'I': 126.9045
}

def calculate_molecular_formula(smiles):
    """
    Calculate molecular formula from SMILES string.
    This is a simplified implementation that:
    1. Counts atom types from the SMILES string
    2. Makes adjustments for implicit hydrogens
    3. Returns a formula string in Hill notation (C first, then H, then others alphabetically)
    
    Note: This is a rough approximation and may not be perfect for complex molecules.
    """
    if not smiles or not isinstance(smiles, str):
        return None
    
    # Clean up the SMILES string (remove charges, stereochemistry, etc.)
    # This is a simplified approach
    clean_smiles = re.sub(r'[+\-\[\]()]', '', smiles)
    clean_smiles = re.sub(r'@|/', '', clean_smiles)
    
    # Extract atoms (single letters followed by numbers, or single letters)
    atom_dict = {}
    
    # Define pattern for atom extraction
    # This handles both single-letter elements (C, H, O, N) and 
    # two-letter elements (Cl, Br, etc.)
    atom_pattern = r'([A-Z][a-z]?)(\d*)'
    
    for match in re.finditer(atom_pattern, clean_smiles):
        atom, count = match.groups()
        count = int(count) if count else 1
        atom_dict[atom] = atom_dict.get(atom, 0) + count
    
    # Estimate implicit hydrogens - this is a very rough approximation
    # For a proper implementation, we would need to analyze the molecular structure
    c_count = atom_dict.get('C', 0)
    o_count = atom_dict.get('O', 0)
    n_count = atom_dict.get('N', 0)
    
    # Rough estimate: each C can have up to 4 bonds, each N up to 3, each O up to 2
    # Implicit H is often what's needed to fill these bonds
    # This is very approximate and will be wrong for many molecules, especially those with rings
    implicit_h = (4 * c_count + 3 * n_count + 2 * o_count) - (2 * (c_count + o_count + n_count - 1))
    implicit_h = max(0, implicit_h)  # Don't allow negative H count
    
    # Add explicit H if present in SMILES
    h_count = atom_dict.get('H', 0) + implicit_h
    if h_count > 0:
        atom_dict['H'] = h_count
    
    # Build formula string in Hill notation (C, H, then others alphabetically)
    formula_parts = []
    
    # Carbon first
    if 'C' in atom_dict:
        count = atom_dict.pop('C')
        formula_parts.append(f"C{count if count > 1 else ''}")
    
    # Hydrogen second
    if 'H' in atom_dict:
        count = atom_dict.pop('H')
        formula_parts.append(f"H{count if count > 1 else ''}")
    
    # Then other elements alphabetically
    for atom in sorted(atom_dict.keys()):
        count = atom_dict[atom]
        formula_parts.append(f"{atom}{count if count > 1 else ''}")
    
    formula = "".join(formula_parts)
    logger.debug(f"Calculated formula: {formula} from SMILES: {smiles}")
    
    return formula

def calculate_molecular_weight(smiles):
    """
    Estimate molecular weight from SMILES string using atomic weights.
    """
    if not smiles or not isinstance(smiles, str):
        return None
    
    # Use the formula calculator to get an estimated molecular formula
    formula = calculate_molecular_formula(smiles)
    if not formula:
        return None
    
    # Parse the formula to calculate molecular weight
    total_weight = 0
    
    # This regex will match element symbols followed by optional numbers
    # e.g., C2H6O will match C2, H6, O1
    formula_pattern = r'([A-Z][a-z]?)(\d*)'
    
    for match in re.finditer(formula_pattern, formula):
        element, count = match.groups()
        count = int(count) if count else 1
        
        if element in ATOMIC_WEIGHTS:
            total_weight += ATOMIC_WEIGHTS[element] * count
        else:
            # Default weight for unknown elements
            total_weight += 12.0 * count
            logger.warning(f"Unknown element {element} in formula {formula}, using default weight")
    
    # Round to 2 decimal places
    return round(total_weight, 2)

def calculate_logp(smiles):
    """
    Estimate LogP (octanol-water partition coefficient) from SMILES.
    This is a simple heuristic based on fragment contributions.
    """
    if not smiles or not isinstance(smiles, str):
        return None
    
    # Count hydrophobic vs. hydrophilic features
    # Carbon and halogens increase logP (lipophilicity)
    lipophilic = (
        smiles.count('C') + smiles.count('c') + 
        0.5 * (smiles.count('Cl') + smiles.count('Br') + smiles.count('I')) +
        0.2 * smiles.count('F')
    )
    
    # Oxygen and nitrogen decrease logP (increase hydrophilicity)
    hydrophilic = (
        1.5 * (smiles.count('O') + smiles.count('o')) + 
        1.0 * (smiles.count('N') + smiles.count('n')) +
        0.5 * smiles.count('S')
    )
    
    # Basic estimation formula
    logp = 0.5 * lipophilic - 1.0 * hydrophilic
    
    # Additional adjustments for common groups
    if 'OH' in smiles or 'oh' in smiles:
        logp -= 0.5  # Hydroxyl groups reduce logP
    
    if 'COOH' in smiles or 'C(=O)O' in smiles:
        logp -= 0.8  # Carboxylic acid groups significantly reduce logP
    
    if 'NH2' in smiles:
        logp -= 0.5  # Amine groups reduce logP
    
    # Round to 2 decimal places
    return round(logp, 2)

def calculate_tpsa(smiles):
    """
    Estimate TPSA (Topological Polar Surface Area) from SMILES.
    TPSA is sum of surface contributions of polar fragments (N, O, etc.)
    """
    if not smiles or not isinstance(smiles, str):
        return None
    
    # Simplified contribution values for polar atoms
    tpsa = (
        30.5 * (smiles.count('N') + smiles.count('n')) +  # Nitrogen contribution
        28.5 * (smiles.count('O') + smiles.count('o')) +  # Oxygen contribution
        34.0 * smiles.count('S') +  # Sulfur contribution
        30.0 * smiles.count('P')    # Phosphorus contribution
    )
    
    # Additional adjustments for common functional groups
    if 'OH' in smiles:
        tpsa += 10.0  # Additional for hydroxyl
    
    if 'NH2' in smiles:
        tpsa += 15.0  # Additional for primary amine
    
    if 'C(=O)O' in smiles:
        tpsa += 15.0  # Additional for carboxylic acid
    
    # Round to 2 decimal places
    return round(tpsa, 2)

def calculate_h_donors(smiles):
    """
    Estimate the number of hydrogen bond donors from SMILES.
    H-bond donors are typically O-H or N-H groups.
    """
    if not smiles or not isinstance(smiles, str):
        return None
    
    # Count OH and NH occurrences
    oh_count = smiles.count('OH') + smiles.count('oh')
    nh_count = smiles.count('NH') + smiles.count('nh')
    
    # Additional checks for NH2, NH3, etc.
    if 'NH2' in smiles:
        nh_count += 1
    
    if 'NH3' in smiles:
        nh_count += 1
    
    # Default minimum value
    h_donors = max(1, oh_count + nh_count)
    
    # For very large molecules, add a proportional factor
    c_count = smiles.count('C') + smiles.count('c')
    if c_count > 20:
        # For large molecules, estimate additional donors based on size
        h_donors += max(0, int(c_count / 20))
    
    return h_donors

def calculate_h_acceptors(smiles):
    """
    Estimate the number of hydrogen bond acceptors from SMILES.
    H-bond acceptors are typically N and O atoms.
    """
    if not smiles or not isinstance(smiles, str):
        return None
    
    # Count N and O atoms as potential acceptors
    n_count = smiles.count('N') + smiles.count('n')
    o_count = smiles.count('O') + smiles.count('o')
    
    # Basic acceptor count
    h_acceptors = n_count + o_count
    
    # Add minimum value
    h_acceptors = max(1, h_acceptors)
    
    return h_acceptors

def calculate_rotatable_bonds(smiles):
    """
    Estimate the number of rotatable bonds from SMILES.
    This is a simplified approximation based on SMILES pattern analysis.
    """
    if not smiles or not isinstance(smiles, str):
        return None
    
    # A simplified approach: estimate based on molecule size and bond types
    c_count = smiles.count('C') + smiles.count('c')
    
    # Rotatable bonds often scale with the number of carbon atoms
    # but not linearly (rings and multiple bonds reduce rotatability)
    if c_count <= 2:
        rot_bonds = 0
    elif c_count <= 5:
        rot_bonds = c_count - 2
    else:
        # For larger molecules, not all C-C bonds are rotatable
        rot_bonds = int(c_count * 0.5)
    
    # Adjust for rings (reduce rotatable bonds)
    ring_count = smiles.count('1') + smiles.count('2') + smiles.count('3') 
    rot_bonds = max(0, rot_bonds - ring_count * 2)
    
    # Final value should be at least 1 for most organic molecules
    rot_bonds = max(1, rot_bonds)
    
    return rot_bonds

def calculate_ring_count(smiles):
    """
    Estimate the number of rings from SMILES.
    Based on the number of ring closure digits in the SMILES.
    """
    if not smiles or not isinstance(smiles, str):
        return None
    
    # Count ring closures (numbers in SMILES often indicate ring closures)
    # SMILES uses numbers and % notation for ring closures
    digits = set('123456789')
    ring_closures = sum(1 for c in smiles if c in digits)
    
    # Each ring needs 2 closures in SMILES, so divide by 2
    ring_count = max(0, ring_closures // 2)
    
    # If we have % ring closures, count those too (for >9 rings)
    ring_count += smiles.count('%')
    
    return ring_count

def calculate_aromatic_ring_count(smiles):
    """
    Estimate the number of aromatic rings from SMILES.
    Based on lowercase c, n, o, etc. which indicate aromatic atoms.
    """
    if not smiles or not isinstance(smiles, str):
        return None
    
    # Count lowercase letters which typically represent aromatic atoms
    aromatic_atoms = sum(1 for c in smiles if c in 'cnops')
    
    # Estimate aromatic rings based on aromatic atom count
    # An aromatic ring typically has 5-6 aromatic atoms
    ring_estimate = max(0, aromatic_atoms // 5)
    
    return ring_estimate

def calculate_heavy_atom_count(smiles):
    """
    Count the number of heavy atoms (non-hydrogen) from SMILES.
    """
    if not smiles or not isinstance(smiles, str):
        return None
    
    # Count uppercase letters which typically start element symbols
    heavy_atoms = sum(1 for c in smiles if c.isupper())
    
    # Count lowercase letters that aren't preceded by uppercase (aromatic atoms)
    for i in range(1, len(smiles)):
        if smiles[i].islower() and not smiles[i-1].isupper():
            heavy_atoms += 1
    
    # Add first character if lowercase (typically aromatic c)
    if len(smiles) > 0 and smiles[0].islower():
        heavy_atoms += 1
    
    return max(1, heavy_atoms)

def test_formula_calculator():
    """Test the formula calculator with some known molecules."""
    test_cases = [
        ("C", "CH4"),  # Methane
        ("CC", "C2H6"),  # Ethane
        ("CCO", "C2H6O"),  # Ethanol
        ("c1ccccc1", "C6H6"),  # Benzene
        ("C1CCCCC1", "C6H12"),  # Cyclohexane
        ("CC(=O)O", "C2H4O2"),  # Acetic acid
        ("c1ccccc1C(=O)O", "C7H6O2"),  # Benzoic acid
        ("CCN", "C2H7N"),  # Ethylamine
        ("c1cc(O)ccc1", "C6H6O"),  # Phenol
        ("ClCCCl", "C2H4Cl2"),  # 1,2-Dichloroethane
    ]
    
    results = []
    for smiles, expected in test_cases:
        calculated = calculate_molecular_formula(smiles)
        match = calculated == expected
        results.append((smiles, expected, calculated, match))
        print(f"SMILES: {smiles}, Expected: {expected}, Calculated: {calculated}, Match: {match}")
    
    # Count matches
    successful = sum(1 for _, _, _, match in results if match)
    print(f"Success rate: {successful}/{len(test_cases)} ({successful*100/len(test_cases):.1f}%)")
    
    return results

def test_property_calculator():
    """Test the property calculator with some known molecules."""
    test_cases = [
        # SMILES, Name for readability
        ("C", "Methane"),
        ("CC", "Ethane"),
        ("CCO", "Ethanol"),
        ("c1ccccc1", "Benzene"),
        ("C1CCCCC1", "Cyclohexane"),
        ("CC(=O)O", "Acetic acid"),
        ("c1ccccc1C(=O)O", "Benzoic acid"),
        ("CCN", "Ethylamine"),
        ("c1cc(O)ccc1", "Phenol"),
        ("ClCCCl", "1,2-Dichloroethane"),
    ]
    
    print("\nTesting property calculations:")
    print("------------------------------")
    for smiles, name in test_cases:
        mw = calculate_molecular_weight(smiles)
        logp = calculate_logp(smiles)
        tpsa = calculate_tpsa(smiles)
        h_donors = calculate_h_donors(smiles)
        h_acceptors = calculate_h_acceptors(smiles)
        rot_bonds = calculate_rotatable_bonds(smiles)
        ring_count = calculate_ring_count(smiles)
        aromatic_rings = calculate_aromatic_ring_count(smiles)
        heavy_atoms = calculate_heavy_atom_count(smiles)
        
        print(f"\n{name} ({smiles}):")
        print(f"  Molecular Weight: {mw}")
        print(f"  LogP: {logp}")
        print(f"  TPSA: {tpsa}")
        print(f"  H-Donors: {h_donors}")
        print(f"  H-Acceptors: {h_acceptors}")
        print(f"  Rotatable Bonds: {rot_bonds}")
        print(f"  Ring Count: {ring_count}")
        print(f"  Aromatic Rings: {aromatic_rings}")
        print(f"  Heavy Atoms: {heavy_atoms}")

if __name__ == "__main__":
    # Run tests when executed directly
    print("Testing mock RDKit molecular formula calculator")
    test_formula_calculator()
    test_property_calculator()