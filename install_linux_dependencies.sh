#!/bin/bash
# CryoProtect - Complete Linux Dependencies Installation Script
# This script installs all required system packages and dependencies for CryoProtect

set -e  # Exit on error

echo "================================================"
echo "CryoProtect Linux Dependencies Installation"
echo "================================================"

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Detect Linux distribution
if [ -f /etc/os-release ]; then
    . /etc/os-release
    DISTRO=$ID
    echo "Detected Linux distribution: $DISTRO"
else
    echo "Unable to determine Linux distribution. Assuming Debian/Ubuntu compatible."
    DISTRO="ubuntu"
fi

# Update package lists
echo -e "\n[1/9] Updating package lists..."
case $DISTRO in
    ubuntu|debian|pop|mint|kali|elementary)
        sudo apt update
        ;;
    fedora|rhel|centos|rocky)
        sudo dnf check-update || true  # dnf returns non-zero if updates are available
        ;;
    arch|manjaro|endeavouros)
        sudo pacman -Sy
        ;;
    *)
        echo "Unsupported distribution. Please install packages manually."
        exit 1
        ;;
esac
echo "✅ Package lists updated"

# Install required system packages
echo -e "\n[2/9] Installing required system packages..."
case $DISTRO in
    ubuntu|debian|pop|mint|kali|elementary)
        sudo apt install -y git curl wget build-essential python3-dev postgresql postgresql-contrib libpq-dev nodejs npm
        ;;
    fedora|rhel|centos|rocky)
        sudo dnf install -y git curl wget gcc gcc-c++ make python3-devel postgresql postgresql-server postgresql-devel nodejs
        # Initialize PostgreSQL if not already done
        if [ ! -f /var/lib/pgsql/data/pg_hba.conf ]; then
            sudo postgresql-setup --initdb
            sudo systemctl enable postgresql
            sudo systemctl start postgresql
        fi
        ;;
    arch|manjaro|endeavouros)
        sudo pacman -S --noconfirm git curl wget base-devel python postgresql nodejs npm
        # Initialize PostgreSQL if not already done
        if [ ! -d /var/lib/postgres/data ]; then
            sudo -u postgres initdb -D /var/lib/postgres/data
            sudo systemctl enable postgresql
            sudo systemctl start postgresql
        fi
        ;;
esac
echo "✅ System packages installed"

# Start PostgreSQL service
echo -e "\n[3/9] Ensuring PostgreSQL service is running..."
if systemctl is-active --quiet postgresql; then
    echo "PostgreSQL is already running"
else
    echo "Starting PostgreSQL service..."
    sudo systemctl start postgresql
fi

# Enable PostgreSQL service to start at boot
sudo systemctl enable postgresql
echo "✅ PostgreSQL service running and enabled at boot"

# Create PostgreSQL user if it doesn't exist
echo -e "\n[4/9] Setting up PostgreSQL user and database..."
if ! sudo -u postgres psql -tAc "SELECT 1 FROM pg_roles WHERE rolname='$(whoami)'" | grep -q 1; then
    echo "Creating PostgreSQL user for $(whoami)..."
    sudo -u postgres createuser --superuser $(whoami)
    
    # Set password
    echo "Please set a password for your PostgreSQL user:"
    sudo -u postgres psql -c "ALTER USER $(whoami) WITH PASSWORD 'password';"
    echo "You can change this password later with: ALTER USER $(whoami) WITH PASSWORD 'newpassword';"
else
    echo "PostgreSQL user $(whoami) already exists"
fi

# Create CryoProtect database if it doesn't exist
if ! sudo -u postgres psql -lqt | cut -d \| -f 1 | grep -qw cryoprotect; then
    echo "Creating CryoProtect database..."
    sudo -u postgres createdb cryoprotect
else
    echo "CryoProtect database already exists"
fi
echo "✅ PostgreSQL user and database setup complete"

# Install Miniconda if not already installed
echo -e "\n[5/9] Setting up Miniconda..."
if ! command_exists conda; then
    echo "Installing Miniconda..."
    MINICONDA_PATH="$HOME/miniconda3"
    
    # Download and install Miniconda
    mkdir -p ~/Downloads
    cd ~/Downloads
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $MINICONDA_PATH
    
    # Add conda to path for current session
    export PATH="$MINICONDA_PATH/bin:$PATH"
    
    # Add conda to path permanently
    echo 'export PATH="$HOME/miniconda3/bin:$PATH"' >> ~/.bashrc
    
    # Initialize conda for bash
    $MINICONDA_PATH/bin/conda init bash
    
    # Source bashrc to get conda in the current session
    source ~/.bashrc
    
    # Update conda
    conda update -y conda
    
    echo "Miniconda installed at $MINICONDA_PATH"
else
    echo "Miniconda is already installed"
fi

# Ensure conda is in PATH
if ! command_exists conda; then
    echo "Conda installed but not in PATH. Please restart your terminal after the script completes."
    export PATH="$HOME/miniconda3/bin:$PATH"
fi
echo "✅ Miniconda setup complete"

# Create the conda environment if it doesn't exist
echo -e "\n[6/9] Creating conda environment for CryoProtect..."
if ! conda env list | grep -q "cryoprotect"; then
    echo "Creating new conda environment 'cryoprotect'..."
    conda create -y -n cryoprotect python=3.9
else
    echo "Conda environment 'cryoprotect' already exists"
fi
echo "✅ Conda environment created"

# Activate the conda environment and install packages
echo -e "\n[7/9] Installing Python dependencies in the conda environment..."
eval "$(conda shell.bash hook)"
conda activate cryoprotect

# Install RDKit with conda
echo "Installing RDKit through conda..."
conda install -y -c conda-forge rdkit=2023.9.1

# Install Python dependencies from requirements.txt
cd "$(dirname "$0")"  # Make sure we're in the project directory
echo "Installing Python packages from requirements.txt..."
if [ -f requirements.txt ]; then
    pip install -r requirements.txt
else
    echo "⚠️ requirements.txt not found. Please install requirements manually."
fi
echo "✅ Python dependencies installed"

# Set up .env file
echo -e "\n[8/9] Setting up environment file..."
if [ ! -f .env ] && [ -f .env.template ]; then
    cp .env.template .env
    echo "Created .env file from template"
    echo "⚠️ Please edit .env file to add your Supabase credentials"
elif [ ! -f .env ] && [ ! -f .env.template ]; then
    echo "⚠️ No .env.template found. Please create a .env file manually."
else
    echo ".env file already exists"
fi
echo "✅ Environment file setup complete"

# Make scripts executable
echo -e "\n[9/9] Making scripts executable..."
find . -name "*.sh" -exec chmod +x {} \;
echo "✅ All shell scripts are now executable"

echo -e "\n================================================"
echo "Installation Complete!"
echo "================================================"

echo -e "\nNext steps:"
echo "1. Edit the .env file with your Supabase credentials"
echo "2. Apply database migrations: node migrations/apply_migration.js"
echo "3. Run the application: ./run_app_with_fix.sh"
echo -e "\nIf you encounter any issues, please check LINUX_MIGRATION_GUIDE.md\n"

# Remind to restart terminal
echo "NOTE: Please restart your terminal or run 'source ~/.bashrc' to use conda commands"