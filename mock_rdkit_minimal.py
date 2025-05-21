#!/usr/bin/env python3
"""
Mock RDKit Service for CryoProtect - Minimal Version
"""

from flask import Flask, jsonify, request

app = Flask(__name__)

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        'status': 'healthy',
        'rdkit_available': False,
        'rdkit_version': 'mock-2022.09.5',
        'rdkit_type': 'mock',
        'environment': 'minimal'
    })

@app.route('/calculate/<smiles>', methods=['GET'])
def calculate_properties(smiles):
    """Calculate molecular properties using mock RDKit"""
    try:
        # Simple property estimation based on SMILES length and composition
        molecular_weight = len(smiles) * 10  # Rough approximation based on length
        o_count = smiles.count('O')
        n_count = smiles.count('N')
        c_count = smiles.count('C')
        logp = (c_count * 0.4) - (o_count * 0.3) - (n_count * 0.2)
        h_donors = o_count + n_count
        h_acceptors = o_count + n_count
        
        properties = {
            'molecular_weight': round(molecular_weight, 2),
            'logp': round(logp, 2),
            'tpsa': round((o_count + n_count) * 20, 2),
            'h_donors': h_donors,
            'h_acceptors': h_acceptors,
            'rotatable_bonds': max(1, len(smiles) // 5),
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
