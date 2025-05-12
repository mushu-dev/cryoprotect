# Molecule Name Standardization Conventions

This document establishes naming conventions for molecules in the CryoProtect database to ensure consistency, clarity, and proper differentiation between similar compounds.

## General Principles

1. **Canonical Naming**: Use the most widely accepted name for common compounds
2. **Capitalization**: Capitalize only the first letter of the name
3. **IUPAC for Specificity**: Use IUPAC names when specificity is required
4. **PubChem Alignment**: Prefer names that match PubChem's primary name
5. **Stereochemistry**: Include stereochemistry indicators for specific isomers
6. **Consistency**: Use consistent formatting for similar types of molecules

## Name Structure Standards

### 1. Common Cryoprotectants

Use the widely recognized common name for well-known cryoprotectants:

- Glycerol (not propane-1,2,3-triol)
- Ethylene glycol (not ethane-1,2-diol)
- Propylene glycol (not propane-1,2-diol)
- DMSO or Dimethyl sulfoxide (not methylsulfinylmethane)
- Trehalose (not α,α-trehalose)

### 2. Amino Acids

Use the standard three-letter code with stereochemistry specified:

- Glycine (achiral)
- L-Alanine (not alanine or (S)-alanine)
- D-Alanine (not (R)-alanine)

For amino acid derivatives, use the format: `[prefix]-[amino acid]`:
- N-Acetyl-L-cysteine
- O-Phospho-L-serine

### 3. Stereochemistry and Isomerism

Include stereochemistry in the name when relevant:

- L-Glutamic acid (not just glutamic acid)
- D-Glucose (not just glucose)
- cis-9-Octadecenoic acid (not just octadecenoic acid)

### 4. Salts and Ionization States

Name format: `[compound name] [cation/anion]`

- Sodium chloride (not NaCl)
- Potassium glutamate (not glutamic acid potassium salt)
- Calcium acetate (not Ca(CH3COO)2)

For acid/base forms:
- Glutamic acid (neutral)
- Glutamate (negatively charged form)
- Protonated glycine (positively charged form)

### 5. Hydration States

Include hydration state for compounds where it affects cryoprotective properties:

- Trehalose dihydrate (not just trehalose with water molecules)
- Anhydrous glycerol (for explicitly water-free)
- Glucose monohydrate

### 6. Mixtures and Blends

For defined mixtures, use the format: `[compound] / [compound] ([ratio])`:

- Glycerol/DMSO (50:50)
- Propylene glycol/ethanol (7:3)

### 7. Compounds with Functional Groups

For compounds with important functional groups, ensure the functional group is clearly indicated:

- 2-Mercaptoethanol (not just mercaptoethanol)
- 3-Hydroxypropanoic acid (not just hydroxypropanoic acid)

### 8. Polymers

For polymers, include average molecular weight or degree of polymerization:

- Polyethylene glycol 4000 (not just PEG)
- Polyvinylpyrrolidone K30 (not just PVP)

### 9. Test Compounds

All test compounds should follow this format:
- TEST_[compound name]

### 10. Duplicate Resolution

When molecules appear to be duplicates but have different PubChem CIDs or structures:

1. Check if they are actually different isomers, salts, or hydration states
2. If different, add appropriate qualifiers: `[name] ([qualifier])`
   - Glycerol (pharmaceutical grade)
   - Glycerol (technical grade)
   - Trehalose (alpha form)
   - Trehalose (beta form)

3. If exact duplicates, keep the one with the most complete data and mark the others for consolidation

## Implementation Guidelines

When implementing these standards:

1. Compare the current name with the standardized form
2. Update names that do not match the standard
3. Record the original name in the properties field
4. For every rename, update the modification_history

Automatic name standardization should be applied with caution:
- Manually review any changes to names of common cryoprotectants
- Use automated processes for systematic cases (capitalization, formatting)
- Maintain a mapping of old names to new names for reference

## Special Cases

### Compounds with "None" as Name

For compounds with "None" as name:
1. If PubChem CID exists:
   - Retrieve the preferred name from PubChem
   - Use IUPAC name if no common name exists
2. If no PubChem CID:
   - Generate name from SMILES if possible
   - Use molecular formula as fallback with prefix "Compound-"
   - Example: "Compound-C6H12O6" for a glucose with no proper name

### Compounds with Multiple Common Names

When a compound has multiple common names, prioritize:
1. Name used in cryobiology literature
2. PubChem's preferred name
3. Most commonly used name in scientific literature
4. IUPAC systematic name