#!/usr/bin/env python3
"""
Mock RDKit Service for CryoProtect
This service provides basic RDKit-like functionality without requiring the actual RDKit package.
"""

import os
from flask import Flask, jsonify, request

# Mock RDKit availability
RDKIT_AVAILABLE = False
RDKIT_VERSION = "mock-2022.09.5"

app = Flask(__name__)

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'rdkit_available': RDKIT_AVAILABLE,
        'rdkit_version': RDKIT_VERSION,
        'rdkit_type': 'mock',
        'environment': os.environ.get('FLASK_ENV', 'development')
    })

@app.route('/calculate/<smiles>', methods=['GET'])
def calculate_properties(smiles):
    """Calculate molecular properties using mock RDKit"""
    try:
        # Return mock properties
        molecular_weight = 0.0
        logp = 0.0
        
        # Simple property estimation based on SMILES length and composition
        if smiles:
            # Basic property estimation rules (very simplified)
            molecular_weight = len(smiles) * 10  # Rough approximation based on length
            
            # LogP estimation: more Os and Ns make it more hydrophilic (negative logP)
            o_count = smiles.count('O')
            n_count = smiles.count('N')
            c_count = smiles.count('C')
            logp = (c_count * 0.4) - (o_count * 0.3) - (n_count * 0.2)
            
            # H-bond donors (simplified: count O and N)
            h_donors = o_count + n_count
            
            # H-bond acceptors (simplified: count O and N)
            h_acceptors = o_count + n_count
            
            properties = {
                'molecular_weight': round(molecular_weight, 2),
                'logp': round(logp, 2),
                'tpsa': round((o_count + n_count) * 20, 2),  # Simple TPSA approximation
                'h_donors': h_donors,
                'h_acceptors': h_acceptors,
                'rotatable_bonds': max(1, len(smiles) // 5),  # Simple approximation
                'calculation_method': 'mock'
            }
        else:
            properties = {
                'molecular_weight': 0.0,
                'logp': 0.0,
                'tpsa': 0.0,
                'h_donors': 0,
                'h_acceptors': 0,
                'rotatable_bonds': 0,
                'calculation_method': 'mock'
            }
        
        return jsonify({
            'smiles': smiles,
            'properties': properties
        })
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)