# RDKit-Based Cryoprotectant Property Analysis

This document explains how we use RDKit to calculate molecular properties that predict a compound's potential as a cryoprotectant.

## Introduction

Cryoprotectants are substances that prevent ice crystal formation during freezing, protecting cells and tissues from damage. Their effectiveness depends on several molecular properties that can be calculated using RDKit's cheminformatics capabilities.

Our implementation uses a unified approach through the `rdkit_wrapper` module, which provides both authentic RDKit calculations and a fallback mock implementation when RDKit isn't available.

## Key Molecular Properties for Cryoprotection

### Basic Properties (Calculated by RDKit)

| Property | Description | Cryoprotective Significance |
|----------|-------------|---------------------------|
| `molecular_weight` | Mass of the molecule | Affects cell penetration and vitrification capability |
| `logp` | Octanol-water partition coefficient | Measures hydrophobicity/hydrophilicity balance |
| `tpsa` | Topological Polar Surface Area | Indicates potential for hydrogen bonding with water |
| `h_donors` | Hydrogen bond donors | Critical for disrupting ice lattice formation |
| `h_acceptors` | Hydrogen bond acceptors | Important for water interaction and vitrification |
| `rotatable_bonds` | Number of rotatable bonds | Affects molecular flexibility and solution entropy |
| `ring_count` | Number of rings | Contributes to overall molecular structure |
| `aromatic_ring_count` | Number of aromatic rings | Affects solubility and membrane interactions |
| `heavy_atom_count` | Non-hydrogen atoms | Related to molecular size and complexity |

### Cryoprotectant-Specific Properties (Derived)

1. **H-Bond Donor/Acceptor Ratio** (`h_bond_donor_acceptor_ratio`)
   - Calculation: `h_donors / h_acceptors`
   - Significance: Balanced H-bonding (ratios 0.5-2.0) is optimal for many cryoprotectants
   - Example: Glycerol has 3 donors and 3 acceptors (ratio=1.0), making it very effective

2. **Total H-Bonding Capacity** (`total_h_bonding_capacity`)
   - Calculation: `h_donors + h_acceptors`
   - Significance: Higher values correlate with better water interaction
   - Example: DMSO has 0 donors and 1 acceptor (total=1), but still effective due to strong acceptor

3. **Polarity Index** (`polarity_index`)
   - Calculation: `tpsa / (molecular_weight^(2/3) * 10)`
   - Significance: Indicates distribution of polar surface relative to size
   - Example: Balanced polarity is important for membrane interaction

4. **Membrane Interaction Score** (`membrane_interaction_score`)
   - Calculation: Complex formula weighing LogP, size, and H-bonding
   - Significance: Predicts interaction with cell membranes
   - Example: Penetrating cryoprotectants like DMSO score high

5. **Ice Interaction Potential** (`ice_interaction_potential`)
   - Calculation: Based on hydroxyl groups, TPSA, and molecular size
   - Significance: Predicts ability to disrupt ice crystal formation
   - Example: Polyols like glycerol score high due to multiple hydroxyl groups

6. **Vitrification Potential** (`vitrification_potential`)
   - Calculation: Complex formula based on H-bonding, LogP, and ring structures
   - Significance: Predicts ability to form glass-like state rather than ice crystals
   - Example: DMSO has excellent vitrification properties

7. **Estimated Toxicity** (`estimated_toxicity`)
   - Calculation: Based on LogP, molecular weight, and H-bonding
   - Significance: Predicts potential cellular toxicity
   - Example: Very hydrophobic compounds (high LogP) often show higher toxicity

8. **Overall Cryoprotectant Score** (`cryoprotectant_score`)
   - Calculation: Weighted combination of other properties (scale 0-10)
   - Significance: Overall predicted effectiveness as a cryoprotectant
   - Examples:
     - Glycerol scores ~7.8
     - DMSO scores ~8.2
     - Ethylene glycol scores ~7.5
     - Trehalose scores ~6.9

## Types of Cryoprotectants

### Penetrating Cryoprotectants

These molecules enter cells and prevent intracellular ice formation:

- **Small polyols**: Glycerol, ethylene glycol, propylene glycol
- **Amides**: Formamide, acetamide
- **Sulfoxides**: DMSO
- **Alcohols**: Methanol, ethanol (limited use due to toxicity)

Common properties:
- MW < 100-150 Da
- Moderate LogP (-2 to 0)
- Multiple hydrogen bonding sites
- High ice interaction and vitrification potential

### Non-Penetrating Cryoprotectants

These molecules remain outside cells and control osmotic pressure:

- **Sugars**: Sucrose, trehalose, glucose
- **Polymers**: Polyvinylpyrrolidone (PVP), hydroxyethyl starch
- **Proteins**: Albumin

Common properties:
- Larger molecules (MW > 180 Da, often much higher)
- Very hydrophilic (lower LogP)
- Multiple hydrogen bonding sites
- High vitrification potential but don't enter cells

## Interpreting Results

The `cryoprotectant_score` provides an overall assessment:

- **8-10**: Excellent potential, comparable to known effective cryoprotectants
- **6-8**: Good potential, worth experimental investigation
- **4-6**: Moderate potential, may work in combination with others
- **2-4**: Low potential, likely ineffective alone
- **0-2**: Very poor potential

## Using the RDKit Property Calculator

Our `rdkit_property_calculator.py` script:

1. Connects to the Supabase database
2. Calculates standard RDKit properties for molecules
3. Derives cryoprotectant-specific properties
4. Stores all calculated properties in the database
5. Identifies promising cryoprotectant candidates

Run it with:
```bash
python rdkit_property_calculator.py --known  # Calculate for known cryoprotectants
python rdkit_property_calculator.py --sample 100  # Random sample of 100 molecules
```

## References

1. Elliott, G. D., Wang, S., & Fuller, B. J. (2017). Cryoprotectants: A review of the actions and applications of cryoprotective solutes that modulate cell recovery from ultra-low temperatures. Cryobiology, 76, 74-91.

2. Hubalek, Z. (2003). Protectants used in the cryopreservation of microorganisms. Cryobiology, 46(3), 205-229.

3. Polge, C., Smith, A. U., & Parkes, A. S. (1949). Revival of spermatozoa after vitrification and dehydration at low temperatures. Nature, 164(4172), 666-666.