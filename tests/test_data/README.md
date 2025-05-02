# CryoProtect v2 Test Datasets

This directory contains test datasets for the CryoProtect v2 application. These datasets are designed for rapid integration testing and development, providing a variety of cryoprotectants and mixtures with diverse properties.

## Datasets Overview

### 1. Core Cryoprotectants (`core_cryoprotectants.json`)

A collection of 14 well-known cryoprotectants with diverse properties:

- **Small penetrating cryoprotectants**: DMSO, glycerol, ethylene glycol, propylene glycol, methanol, ethanol, 1,3-propanediol, formamide
- **Large non-penetrating cryoprotectants**: trehalose, sucrose, polyvinyl alcohol, hydroxyethyl starch
- **Natural osmolytes**: proline, betaine

Each molecule includes:
- Basic identifiers (name, SMILES, InChI, InChIKey, formula)
- Physical properties (molecular weight, melting/boiling points, glass transition temperature)
- Chemical properties (LogP, TPSA, hydrogen bonding)
- Description of its use as a cryoprotectant

### 2. Mixtures (`mixtures.json`)

A set of 7 predefined mixtures using compounds from the core dataset:

- **DMSO/EG Mixture**: Standard 1:1 DMSO/Ethylene glycol mixture for vitrification
- **Glycerol/Trehalose Mixture**: Combination of penetrating and non-penetrating cryoprotectants
- **PG/EG/DMSO Mixture**: Triple cryoprotectant mixture for improved vitrification
- **Sucrose/PVA Solution**: Non-penetrating cryoprotectant mixture
- **Methanol/Glycerol Mixture**: Low molecular weight cryoprotectant mixture
- **Proline/Betaine/Trehalose Mixture**: Natural cryoprotectant mixture inspired by extremophile organisms
- **1,3-Propanediol/HES Mixture**: Mixture combining small penetrating and large non-penetrating cryoprotectants

Each mixture includes:
- Name and description
- Component molecules with concentrations
- Predicted properties (glass transition temperature, cryoprotection scores)
- Typical application

### 3. Edge Cases (`edge_cases.json`)

A collection of 10 compounds that test boundary conditions:

- **Extreme sizes**: Water (smallest), Dextran (very large)
- **Extreme LogP values**: Hexane (highly hydrophobic), Perfluorooctane (extremely hydrophobic)
- **Extreme hydrogen bonding**: Phosphoric acid (high density of H-bonds)
- **Extreme reactivity**: Sulfuric acid (strong acid)
- **Ionic compounds**: Sodium hydroxide (strong base), Ammonium chloride (salt)
- **No hydrogen bonding**: Benzene (aromatic), Cyclohexane (cyclic alkane)

Each molecule includes the same properties as the core dataset, plus an "edge_case_type" field indicating what boundary condition it tests.

## Using the Test Datasets

### Loading Data

Use the provided `load_test_data.py` script to load the datasets into the CryoProtect database:

```bash
# Load all datasets
python load_test_data.py

# Load only core cryoprotectants
python load_test_data.py --dataset core

# Load only mixtures
python load_test_data.py --dataset mixtures

# Load only edge cases
python load_test_data.py --dataset edge

# Clear existing test data before loading
python load_test_data.py --clear
```

### Integration Testing

These datasets are designed for testing various aspects of the CryoProtect application:

1. **Core functionality testing**: Use the core dataset to test basic molecule handling, property calculation, and visualization.

2. **Mixture analysis testing**: Use the mixtures dataset to test:
   - Mixture property prediction
   - Compatibility analysis
   - Synergy detection
   - Optimization algorithms

3. **Edge case handling**: Use the edge cases dataset to test how the application handles extreme values and unusual compounds.

4. **Performance testing**: The datasets are small enough for rapid testing but diverse enough to exercise all application features.

### Example Test Cases

1. **Property calculation validation**:
   ```python
   # Test that LogP values are calculated correctly
   def test_logp_calculation():
       molecule = Molecule.get_by_name("DMSO")
       logp = MolecularProperty.get_property(molecule["id"], "LogP")
       assert abs(logp["numeric_value"] - (-1.35)) < 0.1
   ```

2. **Mixture prediction testing**:
   ```python
   # Test that mixture properties are predicted correctly
   def test_mixture_property_prediction():
       mixture = Mixture.get_by_name("DMSO/EG Mixture")
       components = MixtureComponent.get_by_mixture_id(mixture["id"])
       properties = MixtureProperty.predict_mixture_properties(components)
       assert "Cryoprotection Score" in properties
       assert properties["Cryoprotection Score"] > 8.0
   ```

3. **Edge case handling**:
   ```python
   # Test that the application handles extreme molecular weights
   def test_extreme_molecular_weight():
       dextran = Molecule.get_by_name("Dextran")
       water = Molecule.get_by_name("Water")
       
       # Test that both molecules can be processed
       assert calculate_all_properties(dextran["smiles"]) is not None
       assert calculate_all_properties(water["smiles"]) is not None
   ```

## Data Sources and Validation

The test datasets were created based on:

1. **Scientific literature**: Properties and uses of common cryoprotectants were derived from peer-reviewed publications.

2. **Chemical databases**: Molecular identifiers (SMILES, InChI) were validated against PubChem and ChemSpider.

3. **RDKit validation**: All molecules were validated to ensure they can be parsed by RDKit for property calculation.

4. **Application requirements**: The datasets were designed to test all features of the CryoProtect application.

## Limitations and Considerations

1. **Simplified data**: These datasets are simplified for testing purposes and may not represent the full complexity of real-world cryoprotectant research.

2. **Property approximations**: Some properties (especially glass transition temperatures) are approximations based on literature values.

3. **Test scope**: The datasets are designed for integration testing, not for scientific validation of cryoprotectant efficacy.

4. **Database impact**: The `--clear` option should be used with caution as it will remove existing test data from the database.