#!/bin/bash
#
# Minimal Test Environment Setup for CryoProtect
# This script sets up a minimal test environment without volume mounts

set -e

echo "Setting up minimal CryoProtect test environment..."

# Clean up any existing containers
echo "Cleaning up existing containers..."
podman stop cryoprotect-rdkit-minimal 2>/dev/null || true
podman rm cryoprotect-rdkit-minimal 2>/dev/null || true

# Create the network if it doesn't exist
if ! podman network inspect cryoprotect-net &>/dev/null; then
    echo "Creating cryoprotect-net network..."
    podman network create cryoprotect-net
fi

# Create a simple mock RDKit service script
echo "Creating mock RDKit service script..."
cat > mock_rdkit_minimal.py << 'EOF'
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
EOF

# Start the RDKit service container
echo "Starting minimal RDKit service container..."
podman run -d --name=cryoprotect-rdkit-minimal --replace \
    --network=cryoprotect-net \
    -p 5002:5000 \
    python:3.10-slim \
    sh -c "pip install flask && cat > mock_rdkit.py << 'EOT'
$(cat mock_rdkit_minimal.py)
EOT
python mock_rdkit.py"

echo "Waiting for RDKit service to start..."
sleep 3

# Test the RDKit service
echo "Testing RDKit service..."
curl -s http://localhost:5002/health | python -m json.tool

echo "Minimal test environment setup complete."
echo "You can now test by calling the RDKit service at: http://localhost:5002"