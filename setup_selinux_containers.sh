#!/bin/bash
#
# SELinux Configuration for CryoProtect Containers
# This script sets up SELinux to work properly with CryoProtect containers

set -e

echo "Setting up SELinux for CryoProtect containers..."

# Check if running as root
if [ "$EUID" -ne 0 ]; then
  echo "Please run as root or with sudo"
  exit 1
fi

# Check SELinux status
echo "Checking SELinux status..."
SELINUX_STATUS=$(getenforce)
if [ "$SELINUX_STATUS" == "Disabled" ]; then
  echo "SELinux is disabled. No configuration needed."
  exit 0
fi

echo "SELinux is $SELINUX_STATUS"

# Create the cryoprotect policy module
echo "Creating CryoProtect SELinux policy module..."
cat > cryoprotect.te << EOF
module cryoprotect 1.0;

require {
    type container_t;
    type container_file_t;
    type container_runtime_t;
    class file { read write getattr open };
    class dir { read search open };
}

#============= Rules for CryoProtect containers ===============
allow container_t container_file_t:file { read write getattr open };
allow container_t container_file_t:dir { read search open };
EOF

# Compile and install the policy
echo "Compiling SELinux policy..."
checkmodule -M -m -o cryoprotect.mod cryoprotect.te
semodule_package -o cryoprotect.pp -m cryoprotect.mod

echo "Installing SELinux policy..."
semodule -i cryoprotect.pp

# Set up file contexts for CryoProtect directories
echo "Configuring file contexts for CryoProtect directories..."

# For application directory
APP_DIR="/home/mushu/Projects/CryoProtect"
echo "Setting context for $APP_DIR"
semanage fcontext -a -t container_file_t "$APP_DIR(/.*)?"
restorecon -Rv "$APP_DIR"

# Create dedicated data directories with proper contexts
echo "Creating data directories with proper contexts..."

DATA_DIR="/var/cryoprotect/data"
mkdir -p "$DATA_DIR/rdkit" "$DATA_DIR/app"
chown -R mushu:mushu "$DATA_DIR"
semanage fcontext -a -t container_file_t "$DATA_DIR(/.*)?"
restorecon -Rv "$DATA_DIR"

# Create a wrapper script
echo "Creating container wrapper script..."
cat > /usr/local/bin/run-cryoprotect-containers << EOF
#!/bin/bash
#
# Run CryoProtect containers with proper SELinux configuration
#

# Create the network if it doesn't exist
podman network create cryoprotect-net 2>/dev/null || true

# Stop existing containers
podman stop cryoprotect-rdkit cryoprotect-app 2>/dev/null || true
podman rm cryoprotect-rdkit cryoprotect-app 2>/dev/null || true

# RDKit service container
echo "Starting RDKit service container..."
podman run -d --name=cryoprotect-rdkit --replace \\
    --network=cryoprotect-net \\
    -p 5002:5000 \\
    -v "$DATA_DIR/rdkit:/data:Z" \\
    python:3.10-slim \\
    sh -c "cd /app && pip install flask && cat > mock_rdkit.py << 'EOT'
#!/usr/bin/env python3
from flask import Flask, jsonify, request

app = Flask(__name__)

@app.route('/health', methods=['GET'])
def health_check():
    return jsonify({
        'status': 'healthy',
        'rdkit_available': False,
        'rdkit_version': 'mock-2022.09.5',
        'rdkit_type': 'mock',
        'environment': 'production'
    })

@app.route('/calculate/<smiles>', methods=['GET'])
def calculate_properties(smiles):
    try:
        molecular_weight = len(smiles) * 10
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
    app.run(host='0.0.0.0', port=5000)
EOT
python mock_rdkit.py"

# App container
echo "Starting App container..."
podman run -d --name=cryoprotect-app --replace \\
    --network=cryoprotect-net \\
    -p 5001:5000 \\
    -e RDKIT_SERVICE_URL=http://cryoprotect-rdkit:5000 \\
    -v "$DATA_DIR/app:/data:Z" \\
    -v "$APP_DIR:/app:z" \\
    python:3.10-slim \\
    sh -c "cd /app && pip install -r requirements_essential.txt && python app.py"

echo "Containers started successfully."
echo "- App: http://localhost:5001"
echo "- RDKit service: http://localhost:5002"
EOF

chmod +x /usr/local/bin/run-cryoprotect-containers

echo "Creating systemd service files..."
cat > /etc/systemd/system/cryoprotect-containers.service << EOF
[Unit]
Description=CryoProtect Application Containers
After=network.target

[Service]
Type=oneshot
RemainAfterExit=yes
ExecStart=/usr/local/bin/run-cryoprotect-containers
ExecStop=podman stop cryoprotect-app cryoprotect-rdkit
User=mushu
Group=mushu

[Install]
WantedBy=multi-user.target
EOF

systemctl daemon-reload
systemctl enable cryoprotect-containers.service

echo "SELinux container configuration completed successfully."
echo "Use 'systemctl start cryoprotect-containers' to start the containers"
echo "Use 'systemctl stop cryoprotect-containers' to stop the containers"