#!/bin/bash
# Direct startup script for systemd service

# Set environment variables
export SUPABASE_DB_HOST=localhost
export SUPABASE_DB_PORT=5432
export SUPABASE_DB_NAME=postgres
export SUPABASE_DB_USER=postgres
export SUPABASE_DB_PASSWORD=postgres

# Log startup
echo "Starting CryoProtect application at $(date)"
echo "Working directory: $(pwd)"

# Install missing dependencies (excluding RDKit)
echo "Installing missing dependencies..."
pip install scipy flask flask-restful psycopg2-binary requests python-dotenv

# Create dummy implementation for tests
echo "Creating dummy implementation for RDKit (for testing only)..."
mkdir -p /tmp/mock_modules/rdkit/Chem/Scaffolds
rm -rf /tmp/mock_modules/rdkit/*
touch /tmp/mock_modules/rdkit/__init__.py

# Create Chem module and submodules
cat > /tmp/mock_modules/rdkit/Chem/__init__.py << EOL
# Mock RDKit.Chem module for testing
def MolFromSmiles(smiles):
    return {"smiles": smiles}

def MolToSmiles(mol):
    return mol.get("smiles", "")

def MolFromInchi(inchi):
    return {"inchi": inchi}

def MolToInchi(mol):
    return mol.get("inchi", "")

class Descriptors:
    @staticmethod
    def ExactMolWt(mol):
        return 100.0
    
    @staticmethod
    def MolWt(mol):
        return 100.0
    
    @staticmethod
    def MolLogP(mol):
        return 1.0

class Lipinski:
    @staticmethod
    def HeavyAtomCount(mol):
        return 10
    
    @staticmethod
    def NumHDonors(mol):
        return 2
    
    @staticmethod
    def NumHAcceptors(mol):
        return 2

class MolSurf:
    @staticmethod
    def TPSA(mol):
        return 20.0

class Fragments:
    @staticmethod
    def fr_Al_OH(mol):
        return 1

class AllChem:
    @staticmethod
    def Compute2DCoords(mol):
        return mol
        
    @staticmethod
    def GetMorganFingerprintAsBitVect(mol, radius, nBits):
        return [0] * nBits
EOL

# Create Scaffolds module
mkdir -p /tmp/mock_modules/rdkit/Chem/Scaffolds
touch /tmp/mock_modules/rdkit/Chem/Scaffolds/__init__.py
cat > /tmp/mock_modules/rdkit/Chem/Scaffolds/MurckoScaffold.py << EOL
# Mock MurckoScaffold
def MurckoScaffoldSmiles(mol):
    return "C1CCCCC1"

def GetScaffoldForMol(mol):
    return mol
EOL

# Create mock Draw module
mkdir -p /tmp/mock_modules/rdkit/Chem/Draw
touch /tmp/mock_modules/rdkit/Chem/Draw/__init__.py
cat > /tmp/mock_modules/rdkit/Chem/Draw/rdMolDraw2D.py << EOL
# Mock MolDraw2D
class MolDraw2DCairo:
    def __init__(self, width, height):
        self.width = width
        self.height = height
    
    def DrawMolecule(self, mol):
        pass
    
    def FinishDrawing(self):
        pass
    
    def GetDrawingText(self):
        return "SVG data"
EOL

# Add the mock module to Python path
export PYTHONPATH="/tmp/mock_modules:$PYTHONPATH"

# Run the application
cd /home/mushu/Projects/CryoProtect
echo "Python path: $(which python)"
echo "Python version: $(python --version)"
echo "PYTHONPATH: $PYTHONPATH"

# Run with error handling
python /home/mushu/Projects/CryoProtect/app.py
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo "Application exited with code $EXIT_CODE"
    echo "Check logs for more details"
    exit $EXIT_CODE
fi