# Molecular Similarity and Filtering Mechanisms

## Structure-Based Similarity Algorithms

CryoProtect implements advanced algorithms for molecular similarity searching, enabling researchers to identify structurally related compounds with potential cryoprotectant properties.

### Fingerprint-Based Similarity Methods

#### Morgan/Circular Fingerprints (ECFP)
```python
def calculate_morgan_similarity(molecule1, molecule2, radius=2, bits=2048):
    """
    Calculate Morgan fingerprint (ECFP) similarity between two molecules.
    
    Parameters:
    -----------
    molecule1, molecule2 : rdkit.Chem.rdchem.Mol
        RDKit molecule objects
    radius : int
        Fingerprint radius (2=ECFP4, 3=ECFP6)
    bits : int
        Bit vector length
        
    Returns:
    --------
    float
        Tanimoto similarity score (0-1)
    """
    from rdkit.Chem import AllChem
    from rdkit import DataStructs
    
    # Generate fingerprints
    fp1 = AllChem.GetMorganFingerprintAsBitVect(molecule1, radius, nBits=bits)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(molecule2, radius, nBits=bits)
    
    # Calculate similarity
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity
```

This method is particularly effective for cryoprotectant screening as it captures local structural features that often determine hydrogen bonding patterns and water interaction properties.

#### Topological Fingerprints
These fingerprints encode path-based substructural features and provide complementary information to circular fingerprints:

```python
def calculate_topological_similarity(molecule1, molecule2, bits=2048):
    """Calculate topological fingerprint similarity between molecules"""
    from rdkit.Chem import RDKFingerprint
    from rdkit import DataStructs
    
    # Generate fingerprints
    fp1 = RDKFingerprint(molecule1, fpSize=bits)
    fp2 = RDKFingerprint(molecule2, fpSize=bits)
    
    # Calculate similarity
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity
```

#### MACCS Keys
Industry-standard 166-bit structural keys focusing on functional groups highly relevant for cryoprotectant behavior:

```python
def calculate_maccs_similarity(molecule1, molecule2):
    """Calculate MACCS keys similarity between molecules"""
    from rdkit.Chem import MACCSkeys
    from rdkit import DataStructs
    
    # Generate fingerprints
    fp1 = MACCSkeys.GenMACCSKeys(molecule1)
    fp2 = MACCSkeys.GenMACCSKeys(molecule2)
    
    # Calculate similarity
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return similarity
```

### Pharmacophore-Based Methods

For cryoprotectants, specific 3D arrangement of functional groups is critical. Our custom pharmacophore models capture key features:

```python
def generate_cryoprotectant_pharmacophore(molecule):
    """
    Generate a pharmacophore fingerprint optimized for cryoprotectant features.
    
    Captures:
    - Hydrogen bond donors (essential for water interaction)
    - Hydrogen bond acceptors (essential for water interaction)
    - Hydrophobic regions (cell membrane interaction)
    - Conformational flexibility points
    """
    from rdkit.Chem import ChemicalFeatures
    from rdkit import RDConfig
    import os
    
    # Load feature factory
    featFactory = ChemicalFeatures.BuildFeatureFactory(
        os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
    
    # Generate 3D coordinates if not present
    from rdkit.Chem import AllChem
    if not molecule.GetNumConformers():
        AllChem.EmbedMolecule(molecule)
        AllChem.UFFOptimizeMolecule(molecule)
    
    # Find features
    features = featFactory.GetFeaturesForMol(molecule)
    
    # Generate pharmacophore fingerprint
    # Custom implementation for cryoprotectant relevance
    # ...
    
    return pharmacophore_fp
```

### Maximum Common Substructure Analysis

Identifies common structural elements between known effective cryoprotectants:

```python
def find_common_cryoprotectant_scaffolds(molecules):
    """
    Identify common substructures among known cryoprotectants.
    
    Parameters:
    -----------
    molecules : list of rdkit.Chem.rdchem.Mol
        RDKit molecule objects of known cryoprotectants
        
    Returns:
    --------
    list of rdkit.Chem.rdchem.Mol
        Common substructures found
    """
    from rdkit.Chem import rdFMCS
    
    # Find MCS among all molecules
    mcs_result = rdFMCS.FindMCS(
        molecules,
        threshold=0.8,  # Allow some molecules to not match
        timeout=60,     # Maximum computation time
        atomCompare=rdFMCS.AtomCompare.CompareElements,
        bondCompare=rdFMCS.BondCompare.CompareOrder,
        matchValences=True,
        ringMatchesRingOnly=True,
        completeRingsOnly=True
    )
    
    # Convert SMARTS to molecule
    from rdkit import Chem
    mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
    
    return mcs_mol
```

## Molecular Data Pipeline

### End-to-End Data Flow

The CryoProtect molecular data pipeline integrates multiple sources and filtering mechanisms:

```
+----------------+     +----------------+     +----------------+
|                |     |                |     |                |
|  PubChem API   +---->+  ChEMBL API    +---->+  Custom Input  |
|                |     |                |     |                |
+-------+--------+     +-------+--------+     +-------+--------+
        |                      |                      |
        v                      v                      v
+-------+------------------------+------------------------+
|                                                         |
|           Data Normalization & Standardization          |
|                                                         |
+-------+------------------------+------------------------+
        |                        |                        |
        v                        v                        v
+-------+------+    +------------+-------------+    +-----+-------+
|              |    |                          |    |             |
| De-duplicate |    | Structure Standardization|    | Property    |
| Compounds    |    | - Canonical SMILES       |    | Calculation |
|              |    | - Salt stripping         |    |             |
+-------+------+    +------------+-------------+    +-----+-------+
        |                        |                        |
        v                        v                        v
+-------+------------------------+------------------------+
|                                                         |
|                Initial Filtering Stage                  |
|                                                         |
| - Basic physicochemical filters                         |
| - Structural filters (required/excluded fragments)      |
| - Known toxicity and reactivity filters                 |
|                                                         |
+-------+------------------------+------------------------+
        |
        v
+-------+-------------------------------+
|                                       |
|       Cryoprotectant Relevance Score  |
|                                       |
| - Similarity to known cryoprotectants |
| - Property-based scoring              |
| - Machine learning models             |
|                                       |
+-------+-------------------------------+
        |
        v
+-------+-------------------------------+
|                                       |
|     Scientific Database Storage       |
|                                       |
| - Molecules table                     |
| - Molecular_properties table          |
| - Classification and tagging          |
|                                       |
+-------+---------------+---------------+
        |               |
        v               v
+-------+------+    +---+------------+
|              |    |                |
| Search Index |    | Analytics Data |
|              |    |                |
+--------------+    +----------------+
```

### Input Data Processing

1. **Source Data Extraction**:
   - PubChem API extraction with targeted compound classes
   - ChEMBL API extraction focusing on cryoprotectant-relevant activities
   - Custom data input from literature and experimental sources

2. **Data Standardization**:
   ```python
   def standardize_molecules(molecules):
       """
       Standardize molecule structures for consistency.
       
       Operations:
       - Remove salts and solvents
       - Normalize functional groups
       - Generate canonical tautomers
       - Standardize stereochemistry
       - Assign unique identifiers
       """
       from rdkit.Chem.MolStandardize import rdMolStandardize
       
       standardized_molecules = []
       for mol in molecules:
           # Remove salts and get largest fragment
           clean_mol = rdMolStandardize.FragmentParent(mol)
           
           # Normalize functional groups
           normalized_mol = rdMolStandardize.Normalize(clean_mol)
           
           # Generate canonical tautomer
           tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
           canonical_taut = tautomer_enumerator.Canonicalize(normalized_mol)
           
           # Standardize stereochemistry (if needed)
           # ...
           
           standardized_molecules.append(canonical_taut)
       
       return standardized_molecules
   ```

3. **Property Calculation**:
   - Physical properties (MW, logP, TPSA, etc.)
   - Hydrogen bonding features
   - Custom cryoprotectant-specific properties

### Advanced Property Calculations

CryoProtect implements specialized property calculations relevant for cryoprotectant activity:

```python
def calculate_cryoprotectant_properties(molecule):
    """
    Calculate properties specifically relevant for cryoprotectant activity.
    
    Returns:
    --------
    dict
        Dictionary of calculated properties
    """
    properties = {}
    
    # Basic physicochemical properties
    from rdkit.Chem import Descriptors, Lipinski
    properties['molecular_weight'] = Descriptors.MolWt(molecule)
    properties['logP'] = Descriptors.MolLogP(molecule)
    properties['hydrogen_bond_donors'] = Lipinski.NumHDonors(molecule)
    properties['hydrogen_bond_acceptors'] = Lipinski.NumHAcceptors(molecule)
    properties['rotatable_bonds'] = Descriptors.NumRotatableBonds(molecule)
    properties['polar_surface_area'] = Descriptors.TPSA(molecule)
    
    # Cryoprotectant-specific properties
    # 1. Hydroxyl density (number of OH groups per heavy atom)
    from rdkit.Chem import AllChem
    oh_pattern = AllChem.MolFromSmarts('[OH]')
    num_OH = len(molecule.GetSubstructMatches(oh_pattern))
    num_heavy_atoms = molecule.GetNumHeavyAtoms()
    properties['hydroxyl_density'] = num_OH / num_heavy_atoms if num_heavy_atoms > 0 else 0
    
    # 2. Hydrogen bonding capacity
    properties['hydrogen_bonding_capacity'] = (
        properties['hydrogen_bond_donors'] + properties['hydrogen_bond_acceptors']
    )
    
    # 3. Conformational flexibility index
    properties['flexibility_index'] = (
        properties['rotatable_bonds'] / num_heavy_atoms if num_heavy_atoms > 0 else 0
    )
    
    # 4. Water interaction potential
    # Custom calculation based on hydrogen bonding and polar surface area
    properties['water_interaction_potential'] = (
        0.5 * properties['hydrogen_bonding_capacity'] + 
        0.01 * properties['polar_surface_area']
    )
    
    return properties
```

## Scientific Filtering Criteria for Cryoprotectants

### Primary Filtering Parameters

CryoProtect applies scientifically validated filters to identify potential cryoprotectants:

| Property | Range | Rationale |
|----------|-------|-----------|
| Molecular Weight | 32-500 Da | Balance between size and permeability; smaller molecules typically have better penetration |
| LogP | -3.0 to 1.5 | Water solubility with moderate membrane permeability |
| H-Bond Donors | 1-6 | Sufficient H-bond donors to interact with water and biomolecules |
| H-Bond Acceptors | 2-10 | Sufficient H-bond acceptors to disrupt water crystallization |
| Rotatable Bonds | 0-8 | Conformational flexibility for optimal interactions |
| Polar Surface Area | 20-140 Å² | Balance between hydrophilicity and cell permeability |

### Secondary Scientific Criteria

Beyond basic filters, CryoProtect applies more nuanced scientific criteria:

1. **Glass Transition Temperature Prediction**:
   - Predictive models for glass transition temperature (Tg)
   - Higher Tg values preferred (typically > -45°C)
   - Critical for vitrification-based cryopreservation

2. **Cell Membrane Interaction Models**:
   - Predictions of membrane permeability
   - Models for lipid bilayer interaction
   - Balance between penetration and membrane disruption

3. **Ice Crystal Nucleation Inhibition**:
   - Structural features associated with ice nucleation inhibition
   - Hydrogen bonding patterns that disrupt water crystallization
   - Surface topography interference with crystal lattice formation

4. **Toxicity and Biocompatibility Filters**:
   - Exclusion of known cytotoxic structural elements
   - Reactivity prediction to exclude compounds that may form adducts with biomolecules
   - Cell viability impact prediction

## Machine Learning Integration for Molecule Filtering

CryoProtect integrates specialized machine learning models for molecule filtering and scoring:

### 1. Supervised Classification Models

```python
class CryoprotectantClassifier:
    """
    Machine learning model for classifying potential cryoprotectants.
    
    Uses established cryoprotectants as positive examples and
    compounds known to be ineffective as negative examples.
    """
    def __init__(self, model_path=None):
        """Initialize classifier, optionally loading a saved model"""
        from sklearn.ensemble import RandomForestClassifier
        
        if model_path:
            import joblib
            self.model = joblib.load(model_path)
        else:
            self.model = RandomForestClassifier(
                n_estimators=100,
                max_depth=10,
                random_state=42
            )
            
        # Feature names for descriptor calculation
        self.features = [
            'MolWt', 'MolLogP', 'NumHDonors', 'NumHAcceptors', 
            'NumRotatableBonds', 'TPSA', 'NumAromaticRings',
            'FractionCSP3', 'HeavyAtomCount', 'NumAliphaticRings',
            'hydroxyl_density', 'water_interaction_potential', 
            'flexibility_index'
        ]
    
    def calculate_features(self, molecule):
        """Calculate molecular descriptors for the classifier"""
        from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors
        
        features = []
        for feature_name in self.features:
            if hasattr(Descriptors, feature_name):
                feature_value = getattr(Descriptors, feature_name)(molecule)
            elif hasattr(Lipinski, feature_name):
                feature_value = getattr(Lipinski, feature_name)(molecule)
            elif hasattr(rdMolDescriptors, feature_name):
                feature_value = getattr(rdMolDescriptors, feature_name)(molecule)
            elif feature_name == 'hydroxyl_density':
                # Custom feature calculation (example)
                from rdkit.Chem import AllChem
                oh_pattern = AllChem.MolFromSmarts('[OH]')
                num_OH = len(molecule.GetSubstructMatches(oh_pattern))
                num_heavy_atoms = molecule.GetNumHeavyAtoms()
                feature_value = num_OH / num_heavy_atoms if num_heavy_atoms > 0 else 0
            else:
                # Default value for unknown features
                feature_value = 0
                
            features.append(feature_value)
            
        return features
    
    def predict_probability(self, molecule):
        """
        Predict probability of being an effective cryoprotectant.
        
        Parameters:
        -----------
        molecule : rdkit.Chem.rdchem.Mol
            RDKit molecule object
            
        Returns:
        --------
        float
            Probability (0-1) of being an effective cryoprotectant
        """
        features = self.calculate_features(molecule)
        prob = self.model.predict_proba([features])[0][1]  # Probability of positive class
        return prob
    
    def train(self, molecules, labels):
        """
        Train the classifier on provided molecules.
        
        Parameters:
        -----------
        molecules : list of rdkit.Chem.rdchem.Mol
            RDKit molecule objects
        labels : list of int
            Binary labels (1 = effective cryoprotectant, 0 = not effective)
        """
        import numpy as np
        
        # Calculate features for all molecules
        X = np.array([self.calculate_features(mol) for mol in molecules])
        y = np.array(labels)
        
        # Train the model
        self.model.fit(X, y)
    
    def save(self, model_path):
        """Save trained model to disk"""
        import joblib
        joblib.dump(self.model, model_path)
```

### 2. Graph Neural Networks for Property Prediction

CryoProtect uses graph neural networks for direct prediction of cryoprotectant properties from molecular structure:

```python
class CryoprotectantPropertyPredictor:
    """
    Deep learning model for predicting cryoprotectant-relevant properties.
    
    Uses graph neural networks to predict properties directly from
    molecular structure.
    """
    def __init__(self, model_path=None):
        """Initialize model, optionally loading weights from disk"""
        try:
            import torch
            from torch_geometric.nn import GCNConv, global_mean_pool
            import torch.nn.functional as F
            self.torch_available = True
        except ImportError:
            self.torch_available = False
            return
        
        class GNN(torch.nn.Module):
            def __init__(self):
                super(GNN, self).__init__()
                self.conv1 = GCNConv(9, 64)  # 9 = atom features
                self.conv2 = GCNConv(64, 64)
                self.conv3 = GCNConv(64, 64)
                self.lin1 = torch.nn.Linear(64, 32)
                self.lin2 = torch.nn.Linear(32, 5)  # 5 output properties
                
            def forward(self, data):
                x, edge_index, batch = data.x, data.edge_index, data.batch
                
                # Graph convolution layers
                x = F.relu(self.conv1(x, edge_index))
                x = F.relu(self.conv2(x, edge_index))
                x = F.relu(self.conv3(x, edge_index))
                
                # Global pooling
                x = global_mean_pool(x, batch)
                
                # Fully connected layers
                x = F.relu(self.lin1(x))
                x = self.lin2(x)
                
                return x
        
        self.model = GNN()
        
        if model_path and self.torch_available:
            import torch
            self.model.load_state_dict(torch.load(model_path))
            self.model.eval()
    
    def molecule_to_graph(self, molecule):
        """Convert RDKit molecule to PyTorch Geometric graph"""
        if not self.torch_available:
            return None
            
        import torch
        from torch_geometric.data import Data
        
        # Get atom features
        atomic_nums = []
        for atom in molecule.GetAtoms():
            atomic_nums.append(atom.GetAtomicNum())
            
        # One-hot encoding of atom features
        atom_features = []
        for atomic_num in atomic_nums:
            feature = [0] * 9  # 9 atom types: C, N, O, F, P, S, Cl, Br, I
            if atomic_num == 6:  # Carbon
                feature[0] = 1
            elif atomic_num == 7:  # Nitrogen
                feature[1] = 1
            elif atomic_num == 8:  # Oxygen
                feature[2] = 1
            elif atomic_num == 9:  # Fluorine
                feature[3] = 1
            elif atomic_num == 15:  # Phosphorus
                feature[4] = 1
            elif atomic_num == 16:  # Sulfur
                feature[5] = 1
            elif atomic_num == 17:  # Chlorine
                feature[6] = 1
            elif atomic_num == 35:  # Bromine
                feature[7] = 1
            elif atomic_num == 53:  # Iodine
                feature[8] = 1
            atom_features.append(feature)
            
        # Get bond list
        edges = []
        for bond in molecule.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            edges.append([i, j])
            edges.append([j, i])  # Add reverse edge for undirected graph
            
        # Create PyTorch tensors
        x = torch.tensor(atom_features, dtype=torch.float)
        edge_index = torch.tensor(edges, dtype=torch.long).t().contiguous()
        
        # Create graph data object
        data = Data(x=x, edge_index=edge_index)
        
        return data
    
    def predict_properties(self, molecule):
        """
        Predict cryoprotectant-relevant properties from molecule.
        
        Properties predicted:
        1. Glass transition temperature (Tg)
        2. Cell membrane permeability
        3. Ice nucleation inhibition factor
        4. Cytotoxicity estimate
        5. Viscosity impact factor
        
        Parameters:
        -----------
        molecule : rdkit.Chem.rdchem.Mol
            RDKit molecule object
            
        Returns:
        --------
        dict
            Dictionary of predicted properties
        """
        if not self.torch_available:
            return {
                "error": "PyTorch and PyTorch Geometric required for GNN predictions"
            }
            
        import torch
        
        # Convert molecule to graph
        graph = self.molecule_to_graph(molecule)
        if graph is None:
            return {}
            
        # Add batch dimension
        from torch_geometric.data import Batch
        batch = Batch.from_data_list([graph])
        
        # Predict properties
        with torch.no_grad():
            predictions = self.model(batch).numpy()[0]
            
        # Map predictions to property dictionary
        property_names = [
            "glass_transition_temperature",
            "membrane_permeability",
            "ice_inhibition_factor",
            "cytotoxicity_index",
            "viscosity_impact_factor"
        ]
        
        return dict(zip(property_names, predictions))
```

## API and User Interface for Molecule Filtering

### RESTful API Endpoints

CryoProtect provides dedicated API endpoints for molecule filtering and search:

1. **Basic Search Endpoint**:
   ```http
   GET /api/v1/molecules?name=glycerol&property_min=molecular_weight,60&property_max=molecular_weight,200
   ```

2. **Advanced Search Endpoint**:
   ```http
   POST /api/v1/molecules/search
   Content-Type: application/json
   
   {
     "filters": {
       "properties": [
         {"name": "molecular_weight", "min": 60, "max": 200},
         {"name": "logP", "max": 1.0}
       ],
       "structural": {
         "required": ["[OX2H]"],
         "excluded": ["[SX2]"]
       }
     },
     "similarity": {
       "reference": "C(CO)O", 
       "method": "morgan",
       "threshold": 0.6
     },
     "pagination": {"page": 1, "size": 20},
     "sort": {"field": "similarity", "order": "desc"}
   }
   ```

3. **Filtering Preset Endpoint**:
   ```http
   GET /api/v1/molecules/filter/cryoprotectants
   ```
   
   Applies a predefined set of cryoprotectant-relevant filters.

### User Interface Components

The web interface provides interactive components for molecule filtering:

1. **Property Range Sliders**:
   - Adjust molecular weight, logP, etc.
   - Real-time filter updates

2. **Structure Drawing Tool**:
   - Interactive sketcher for similarity searches
   - Template selection for common cryoprotectant scaffolds

3. **Filter Presets**:
   - One-click presets for different cryoprotection applications
   - Customizable and shareable filter configurations

4. **Results Visualization**:
   - Interactive charts showing property distributions
   - Structural clustering of results
   - Similarity network graphs