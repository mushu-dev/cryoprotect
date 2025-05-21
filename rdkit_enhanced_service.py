"""
Enhanced RDKit Service - Provides advanced molecule modeling capabilities.
"""
import os
import json
import logging
import tempfile
from datetime import datetime
from functools import wraps
import traceback

from flask import Flask, request, jsonify, Response
from flask_cors import CORS
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Draw, Descriptors, Descriptors3D, PandasTools
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.ML.Descriptors import MoleculeDescriptors
import redis
import pickle

# Configure RDKit logging
RDLogger.DisableLog('rdApp.*')

# Initialize Flask app
app = Flask(__name__)
CORS(app)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("rdkit_enhanced_service.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Initialize Redis for caching (if available)
try:
    cache = redis.Redis(host='localhost', port=6379, db=0)
    cache_available = True
except:
    cache_available = False
    logger.warning("Redis cache not available, continuing without caching")

# Cache decorator
def cached(ttl=86400):
    """Cache decorator for expensive calculations."""
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if not cache_available:
                return func(*args, **kwargs)
                
            # Create cache key from function name and arguments
            key = f"rdkit:{func.__name__}:{str(args)}:{str(sorted(kwargs.items()))}"
            result = None
            
            try:
                # Try to get from cache
                cached_result = cache.get(key)
                if cached_result:
                    result = pickle.loads(cached_result)
                    logger.debug(f"Cache hit for {key}")
                else:
                    # Calculate and cache
                    result = func(*args, **kwargs)
                    cache.setex(key, ttl, pickle.dumps(result))
                    logger.debug(f"Cache miss for {key}")
            except Exception as e:
                logger.error(f"Cache error: {str(e)}")
                result = func(*args, **kwargs)
                
            return result
        return wrapper
    return decorator

# ---- Utility Functions ----

def mol_from_input(mol_input):
    """Convert various input formats to RDKit Mol object."""
    mol = None
    
    if isinstance(mol_input, str):
        # Try SMILES first
        mol = Chem.MolFromSmiles(mol_input)
        
        # If not SMILES, try InChI
        if mol is None and mol_input.startswith('InChI='):
            mol = Chem.MolFromInchi(mol_input)
            
        # If still None, try SMARTS
        if mol is None:
            mol = Chem.MolFromSmarts(mol_input)
    
    elif isinstance(mol_input, dict) and 'smiles' in mol_input:
        mol = Chem.MolFromSmiles(mol_input['smiles'])
    
    elif isinstance(mol_input, dict) and 'molblock' in mol_input:
        mol = Chem.MolFromMolBlock(mol_input['molblock'])
        
    return mol

def get_descriptors_list():
    """Get list of all available descriptors."""
    # 2D descriptors
    descriptors_2d = {d[0]: d[1].__doc__ for d in Descriptors._descList}
    
    # 3D descriptors
    descriptors_3d = {}
    for name in dir(Descriptors3D):
        obj = getattr(Descriptors3D, name)
        if callable(obj) and not name.startswith('_'):
            descriptors_3d[name] = obj.__doc__
    
    return {
        '2d': descriptors_2d,
        '3d': descriptors_3d
    }

@cached(ttl=3600)
def calculate_all_descriptors(mol, include_3d=False):
    """Calculate all available descriptors for a molecule."""
    if mol is None:
        return None
        
    results = {}
    
    # Calculate 2D descriptors
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    descriptors = calculator.CalcDescriptors(mol)
    
    for i, name in enumerate([x[0] for x in Descriptors._descList]):
        value = descriptors[i]
        # Convert numpy values to Python native types for JSON serialization
        if isinstance(value, (np.int64, np.int32, np.int16, np.int8)):
            value = int(value)
        elif isinstance(value, (np.float64, np.float32, np.float16)):
            value = float(value)
        results[name] = value
    
    # Calculate 3D descriptors if requested
    if include_3d:
        # First need to ensure the molecule has 3D coordinates
        if not mol.GetNumConformers() or mol.GetConformer().Is3D() == False:
            mol = prepare_3d_molecule(mol)
            
        for name in dir(Descriptors3D):
            if name.startswith('_'):
                continue
                
            descriptor_function = getattr(Descriptors3D, name)
            if callable(descriptor_function):
                try:
                    value = descriptor_function(mol)
                    # Convert numpy values
                    if isinstance(value, (np.int64, np.int32, np.int16, np.int8)):
                        value = int(value)
                    elif isinstance(value, (np.float64, np.float32, np.float16)):
                        value = float(value)
                    results[f"3D_{name}"] = value
                except Exception as e:
                    logger.warning(f"Error calculating 3D descriptor {name}: {str(e)}")
                    results[f"3D_{name}"] = None
    
    return results

def prepare_3d_molecule(mol):
    """Prepare a molecule for 3D descriptor calculation."""
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol)
    return mol

@cached(ttl=86400)  # Cache for 24 hours
def generate_conformers(mol, n_conformers=10, max_iterations=1000, random_seed=42):
    """Generate multiple conformers for a molecule."""
    if mol is None:
        return None
        
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate conformers
    conformer_ids = AllChem.EmbedMultipleConfs(
        mol,
        numConfs=n_conformers,
        maxAttempts=max_iterations,
        randomSeed=random_seed,
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True,
        enforceChirality=True,
        numThreads=0  # Use all available CPUs
    )
    
    # Optimize conformers with MMFF
    mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
    if mmff_props is not None:
        for conf_id in conformer_ids:
            AllChem.MMFFOptimizeMolecule(mol, confId=conf_id, mmffVariant='MMFF94s')
    
    # Calculate energy for each conformer
    energies = []
    for conf_id in conformer_ids:
        if mmff_props is not None:
            energy = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=conf_id).CalcEnergy()
        else:
            energy = 0.0  # Fallback if MMFF fails
        energies.append(energy)
    
    # Sort conformers by energy
    sorted_energies = sorted(zip(conformer_ids, energies), key=lambda x: x[1])
    
    # Extract conformers as molblocks
    results = []
    for conf_id, energy in sorted_energies:
        conf_mol = Chem.Mol(mol)
        conf_mol.RemoveAllConformers()
        conf_mol.AddConformer(mol.GetConformer(conf_id))
        results.append({
            'id': int(conf_id),
            'energy': float(energy),
            'molblock': Chem.MolToMolBlock(conf_mol),
            'relative_energy': float(energy - sorted_energies[0][1])
        })
    
    return results

def calculate_pharmacophore_features(mol):
    """Calculate pharmacophore features for a molecule."""
    if mol is None:
        return None
        
    # Ensure molecule has 3D coordinates
    if not mol.GetNumConformers() or mol.GetConformer().Is3D() == False:
        mol = prepare_3d_molecule(mol)
    
    # Define feature factories
    feature_factory = AllChem.BuildFeatureFactory('BaseFeatures.fdef')
    features = feature_factory.GetFeaturesForMol(mol)
    
    result = []
    for feature in features:
        feature_type = feature.GetType()
        feature_positions = feature.GetPos()
        atoms = feature.GetAtomIds()
        
        result.append({
            'type': feature_type,
            'atoms': [int(atom) for atom in atoms],
            'position': [float(feature_positions.x), float(feature_positions.y), float(feature_positions.z)]
        })
    
    return result

def standardize_molecule(mol):
    """Standardize a molecule using RDKit's MolStandardize."""
    if mol is None:
        return None
        
    # Create standardization objects
    clean_mol = rdMolStandardize.Cleanup(mol)
    
    # Get parent molecule after removing salt, solvents
    parent_mol = rdMolStandardize.FragmentParent(clean_mol)
    
    # Normalize functional groups
    normalized_mol = rdMolStandardize.Normalize(parent_mol)
    
    # Apply charge corrections
    charge_corrected = rdMolStandardize.Reionize(normalized_mol)
    
    # Canonicalize tautomers
    tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
    canonical_tautomer = tautomer_enumerator.Canonicalize(charge_corrected)
    
    return canonical_tautomer

# ---- API Endpoints ----

@app.route('/api/v1/rdkit/descriptors', methods=['POST'])
def api_calculate_descriptors():
    """Calculate molecular descriptors."""
    try:
        data = request.get_json()
        
        if not data or 'molecule' not in data:
            return jsonify({'error': 'Missing molecule data'}), 400
            
        include_3d = data.get('include_3d', False)
        
        mol = mol_from_input(data['molecule'])
        if mol is None:
            return jsonify({'error': 'Invalid molecule input'}), 400
            
        descriptors = calculate_all_descriptors(mol, include_3d=include_3d)
        
        return jsonify({
            'molecule': data['molecule'],
            'descriptors': descriptors
        })
        
    except Exception as e:
        logger.error(f"Error calculating descriptors: {str(e)}\n{traceback.format_exc()}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/v1/rdkit/conformers', methods=['POST'])
def api_generate_conformers():
    """Generate 3D conformers for a molecule."""
    try:
        data = request.get_json()
        
        if not data or 'molecule' not in data:
            return jsonify({'error': 'Missing molecule data'}), 400
            
        n_conformers = int(data.get('n_conformers', 10))
        if n_conformers > 100:
            return jsonify({'error': 'Maximum 100 conformers allowed'}), 400
            
        mol = mol_from_input(data['molecule'])
        if mol is None:
            return jsonify({'error': 'Invalid molecule input'}), 400
            
        conformers = generate_conformers(mol, n_conformers=n_conformers)
        
        return jsonify({
            'molecule': data['molecule'],
            'conformers': conformers
        })
        
    except Exception as e:
        logger.error(f"Error generating conformers: {str(e)}\n{traceback.format_exc()}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/v1/rdkit/minimize', methods=['POST'])
def api_minimize_molecule():
    """Energy minimize a molecule."""
    try:
        data = request.get_json()
        
        if not data or 'molecule' not in data:
            return jsonify({'error': 'Missing molecule data'}), 400
            
        mol = mol_from_input(data['molecule'])
        if mol is None:
            return jsonify({'error': 'Invalid molecule input'}), 400
            
        # Add hydrogens
        mol = Chem.AddHs(mol)
        
        # Generate 3D coordinates if not present
        if not mol.GetNumConformers() or mol.GetConformer().Is3D() == False:
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            
        # Minimize with MMFF
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol)
        if mmff_props is None:
            return jsonify({'error': 'MMFF minimization failed'}), 400
            
        force_field = AllChem.MMFFGetMoleculeForceField(mol, mmff_props)
        energy_before = force_field.CalcEnergy()
        
        # Minimize energy
        force_field.Minimize(maxIts=2000)
        energy_after = force_field.CalcEnergy()
        
        return jsonify({
            'molecule': data['molecule'],
            'energy_before': float(energy_before),
            'energy_after': float(energy_after),
            'energy_delta': float(energy_before - energy_after),
            'molblock': Chem.MolToMolBlock(mol)
        })
        
    except Exception as e:
        logger.error(f"Error minimizing molecule: {str(e)}\n{traceback.format_exc()}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/v1/rdkit/similarity/3d', methods=['POST'])
def api_3d_similarity():
    """Calculate 3D similarity between molecules."""
    try:
        data = request.get_json()
        
        if not data or 'query' not in data or 'target' not in data:
            return jsonify({'error': 'Missing query or target molecule'}), 400
            
        query_mol = mol_from_input(data['query'])
        target_mol = mol_from_input(data['target'])
        
        if query_mol is None or target_mol is None:
            return jsonify({'error': 'Invalid molecule input'}), 400
            
        # Prepare 3D structures
        query_mol = prepare_3d_molecule(query_mol)
        target_mol = prepare_3d_molecule(target_mol)
        
        # Calculate 3D similarity using shape-based alignment
        pyO3A = AllChem.GetO3A(query_mol, target_mol)
        score = pyO3A.Score()
        pyO3A.Align()
        
        # Get aligned structure
        aligned_query = Chem.Mol(query_mol)
        
        return jsonify({
            'query': data['query'],
            'target': data['target'],
            'similarity_score': float(score),
            'aligned_query_molblock': Chem.MolToMolBlock(aligned_query)
        })
        
    except Exception as e:
        logger.error(f"Error calculating 3D similarity: {str(e)}\n{traceback.format_exc()}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/v1/rdkit/pharmacophore', methods=['POST'])
def api_pharmacophore():
    """Generate pharmacophore model for a molecule."""
    try:
        data = request.get_json()
        
        if not data or 'molecule' not in data:
            return jsonify({'error': 'Missing molecule data'}), 400
            
        mol = mol_from_input(data['molecule'])
        if mol is None:
            return jsonify({'error': 'Invalid molecule input'}), 400
            
        features = calculate_pharmacophore_features(mol)
        
        return jsonify({
            'molecule': data['molecule'],
            'pharmacophore_features': features
        })
        
    except Exception as e:
        logger.error(f"Error calculating pharmacophore: {str(e)}\n{traceback.format_exc()}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/v1/rdkit/standardize', methods=['POST'])
def api_standardize():
    """Standardize a molecule."""
    try:
        data = request.get_json()
        
        if not data or 'molecule' not in data:
            return jsonify({'error': 'Missing molecule data'}), 400
            
        mol = mol_from_input(data['molecule'])
        if mol is None:
            return jsonify({'error': 'Invalid molecule input'}), 400
            
        standardized = standardize_molecule(mol)
        
        return jsonify({
            'original': data['molecule'],
            'standardized_smiles': Chem.MolToSmiles(standardized),
            'standardized_molblock': Chem.MolToMolBlock(standardized),
            'inchi': Chem.MolToInchi(standardized),
            'inchikey': Chem.MolToInchiKey(standardized)
        })
        
    except Exception as e:
        logger.error(f"Error standardizing molecule: {str(e)}\n{traceback.format_exc()}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/v1/rdkit/batch', methods=['POST'])
def api_batch_process():
    """Process multiple molecules in batch."""
    try:
        data = request.get_json()
        
        if not data or 'molecules' not in data or not isinstance(data['molecules'], list):
            return jsonify({'error': 'Missing or invalid molecules list'}), 400
            
        operation = data.get('operation', 'descriptors')
        include_3d = data.get('include_3d', False)
        
        results = []
        
        for mol_input in data['molecules']:
            mol = mol_from_input(mol_input)
            if mol is None:
                results.append({'input': mol_input, 'error': 'Invalid molecule input'})
                continue
                
            if operation == 'descriptors':
                result = calculate_all_descriptors(mol, include_3d=include_3d)
                results.append({'input': mol_input, 'descriptors': result})
            
            elif operation == 'standardize':
                standardized = standardize_molecule(mol)
                results.append({
                    'input': mol_input,
                    'standardized_smiles': Chem.MolToSmiles(standardized),
                    'inchikey': Chem.MolToInchiKey(standardized)
                })
            
            elif operation == 'conformers':
                n_conformers = int(data.get('n_conformers', 3))
                conformers = generate_conformers(mol, n_conformers=n_conformers)
                results.append({'input': mol_input, 'conformers': conformers})
            
            else:
                return jsonify({'error': 'Invalid operation'}), 400
        
        return jsonify({
            'operation': operation,
            'results': results
        })
        
    except Exception as e:
        logger.error(f"Error in batch processing: {str(e)}\n{traceback.format_exc()}")
        return jsonify({'error': str(e)}), 500

@app.route('/api/v1/rdkit/info', methods=['GET'])
def api_info():
    """Get information about the RDKit service."""
    descriptors_list = get_descriptors_list()
    
    return jsonify({
        'version': '1.0.0',
        'rdkit_version': Chem.__version__,
        'descriptors_count': {
            '2d': len(descriptors_list['2d']),
            '3d': len(descriptors_list['3d']),
            'total': len(descriptors_list['2d']) + len(descriptors_list['3d'])
        },
        'cache_available': cache_available,
        'endpoints': [
            '/api/v1/rdkit/descriptors',
            '/api/v1/rdkit/conformers',
            '/api/v1/rdkit/minimize',
            '/api/v1/rdkit/similarity/3d',
            '/api/v1/rdkit/pharmacophore',
            '/api/v1/rdkit/standardize',
            '/api/v1/rdkit/batch',
            '/api/v1/rdkit/info'
        ]
    })

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint."""
    return jsonify({'status': 'healthy', 'timestamp': datetime.now().isoformat()})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)