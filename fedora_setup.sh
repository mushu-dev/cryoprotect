#!/bin/bash
# CryoProtect - Fedora-specific Setup Script
# This script installs all necessary components for CryoProtect on Fedora Linux

set -e  # Exit on error

BLUE='\033[0;34m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Function to print section headers
print_section() {
    echo -e "\n${BLUE}================================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}================================================${NC}\n"
}

# Function to print success messages
print_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

# Function to print warning messages
print_warning() {
    echo -e "${YELLOW}⚠️ $1${NC}"
}

# Function to print error messages
print_error() {
    echo -e "${RED}❌ $1${NC}"
}

# Function to check if a command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check Fedora version
check_fedora_version() {
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        if [ "$NAME" != "Fedora Linux" ]; then
            print_error "This script is designed for Fedora Linux. Detected: $NAME"
            print_warning "The script may not work properly on this distribution."
            read -p "Continue anyway? (y/n) " -n 1 -r
            echo
            if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                exit 1
            fi
        else
            print_success "Detected Fedora Linux $VERSION_ID"
            # Check for minimum version
            if [ "$VERSION_ID" -lt 36 ]; then
                print_warning "This script is tested on Fedora 36+. You're running Fedora $VERSION_ID."
                print_warning "Some features may not work as expected."
            fi
        fi
    else
        print_error "Cannot detect OS"
        exit 1
    fi
}

# Function to set up SELinux contexts
setup_selinux() {
    print_section "Setting up SELinux contexts"
    
    # Check if SELinux is enabled
    if command_exists getenforce; then
        selinux_status=$(getenforce)
        if [ "$selinux_status" = "Enforcing" ] || [ "$selinux_status" = "Permissive" ]; then
            print_success "SELinux is enabled: $selinux_status"
            
            # Install SELinux utilities if not already installed
            if ! command_exists semanage; then
                print_warning "SELinux policy management tools not found. Installing..."
                sudo dnf install -y policycoreutils-python-utils
            fi
            
            # Set SELinux contexts for project directories
            PROJECT_DIR=$(pwd)
            
            # Set context for application files
            sudo semanage fcontext -a -t httpd_sys_content_t "$PROJECT_DIR(/.*)?"
            
            # Set context for log directory
            mkdir -p "$PROJECT_DIR/logs"
            sudo semanage fcontext -a -t httpd_log_t "$PROJECT_DIR/logs(/.*)?"
            
            # Set context for cache directory
            mkdir -p "$PROJECT_DIR/cache"
            sudo semanage fcontext -a -t httpd_cache_t "$PROJECT_DIR/cache(/.*)?"
            
            # Set context for backup directory
            mkdir -p "$PROJECT_DIR/backup/data"
            sudo semanage fcontext -a -t httpd_sys_rw_content_t "$PROJECT_DIR/backup/data(/.*)?"
            
            # Apply contexts
            sudo restorecon -Rv "$PROJECT_DIR"
            
            print_success "SELinux contexts set"
            
            # Create SELinux policy module for PostgreSQL
            if [ -d "/var/lib/pgsql" ]; then
                echo "
module cryoprotect 1.0;

require {
    type httpd_t;
    type postgresql_t;
    type postgresql_db_t;
    class unix_stream_socket connectto;
    class dir search;
    class file { read write getattr open };
}

#============= httpd_t ==============
allow httpd_t postgresql_t:unix_stream_socket connectto;
allow httpd_t postgresql_db_t:dir search;
allow httpd_t postgresql_db_t:file { read write getattr open };
" > cryoprotect.te
                
                # Compile and install the policy module
                if command_exists checkmodule && command_exists semodule_package; then
                    checkmodule -M -m -o cryoprotect.mod cryoprotect.te
                    semodule_package -o cryoprotect.pp -m cryoprotect.mod
                    sudo semodule -i cryoprotect.pp
                    rm -f cryoprotect.te cryoprotect.mod cryoprotect.pp
                    print_success "SELinux policy module installed"
                else
                    print_warning "SELinux policy tools not found. Skipping policy module installation."
                    print_warning "You may need to run: sudo dnf install policycoreutils-devel"
                fi
            fi
        else
            print_warning "SELinux is disabled. Skipping SELinux configuration."
        fi
    else
        print_warning "SELinux tools not found. Skipping SELinux configuration."
    fi
}

# Function to set up firewall
setup_firewall() {
    print_section "Setting up firewall"
    
    if command_exists firewall-cmd; then
        # Check if firewalld is running
        if sudo systemctl is-active --quiet firewalld; then
            print_success "Firewalld is running"
            
            # Add ports for the application
            sudo firewall-cmd --add-port=5000/tcp --permanent
            
            # If using ELK stack, add those ports too
            sudo firewall-cmd --add-port=9200/tcp --permanent # Elasticsearch
            sudo firewall-cmd --add-port=5601/tcp --permanent # Kibana
            
            # If using monitoring, add those ports
            sudo firewall-cmd --add-port=9090/tcp --permanent # Prometheus
            sudo firewall-cmd --add-port=3000/tcp --permanent # Grafana
            
            # Reload firewall
            sudo firewall-cmd --reload
            
            print_success "Firewall configured"
        else
            print_warning "Firewalld is not running. Skipping firewall configuration."
        fi
    else
        print_warning "Firewall-cmd not found. Skipping firewall configuration."
    fi
}

# Main installation function
install_dependencies() {
    print_section "Installing system dependencies"
    
    # Update the system
    sudo dnf update -y
    print_success "System updated"
    
    # Install package groups
    sudo dnf group install -y "Development Tools" "Python Development"
    print_success "Development tools installed"
    
    # Install specific packages
    sudo dnf install -y \
        git \
        curl \
        wget \
        python3 \
        python3-devel \
        python3-pip \
        postgresql \
        postgresql-server \
        postgresql-devel \
        nodejs \
        npm \
        libXrender \
        libXext \
        libXcomposite \
        libXi \
        libXrandr \
        mesa-libGL \
        cairo-devel \
        pango-devel \
        fontconfig-devel \
        dejavu-sans-fonts \
        dejavu-serif-fonts \
        boost-devel \
        eigen3-devel \
        openssl-devel
    
    print_success "Required packages installed"
}

# Setup PostgreSQL
setup_postgresql() {
    print_section "Setting up PostgreSQL"
    
    # Check if PostgreSQL is initialized
    if [ ! -f /var/lib/pgsql/data/pg_hba.conf ]; then
        print_warning "PostgreSQL is not initialized. Initializing..."
        sudo postgresql-setup --initdb
    else
        print_success "PostgreSQL is already initialized"
    fi
    
    # Start and enable PostgreSQL service
    sudo systemctl enable postgresql
    sudo systemctl start postgresql
    
    if sudo systemctl is-active --quiet postgresql; then
        print_success "PostgreSQL service is running"
    else
        print_error "Failed to start PostgreSQL service"
        exit 1
    fi
    
    # Create database user and database
    echo "Creating database user and database..."
    read -p "Enter a username for PostgreSQL [$(whoami)]: " PG_USER
    PG_USER=${PG_USER:-$(whoami)}
    
    read -sp "Enter a password for PostgreSQL user $PG_USER: " PG_PASS
    echo
    
    # Create user and database
    sudo -u postgres psql -c "CREATE USER $PG_USER WITH SUPERUSER PASSWORD '$PG_PASS';" 2>/dev/null || echo "User $PG_USER may already exist"
    sudo -u postgres psql -c "CREATE DATABASE cryoprotect OWNER $PG_USER;" 2>/dev/null || echo "Database cryoprotect may already exist"
    sudo -u postgres psql -c "CREATE DATABASE cryoprotect_test OWNER $PG_USER;" 2>/dev/null || echo "Database cryoprotect_test may already exist"
    
    # Update pg_hba.conf to allow password authentication for local connections
    sudo sh -c "echo 'local   all             $PG_USER                                 md5' >> /var/lib/pgsql/data/pg_hba.conf"
    sudo sh -c "echo 'host    all             $PG_USER             127.0.0.1/32        md5' >> /var/lib/pgsql/data/pg_hba.conf"
    sudo sh -c "echo 'host    all             $PG_USER             ::1/128             md5' >> /var/lib/pgsql/data/pg_hba.conf"
    
    # Restart PostgreSQL to apply changes
    sudo systemctl restart postgresql
    
    print_success "PostgreSQL setup complete"
    
    # Update .env file with PostgreSQL credentials
    if [ -f .env ]; then
        sed -i "s/LOCAL_DB_USER=postgres/LOCAL_DB_USER=$PG_USER/g" .env
        sed -i "s/LOCAL_DB_PASSWORD=/LOCAL_DB_PASSWORD=$PG_PASS/g" .env
        print_success "Updated .env file with PostgreSQL credentials"
    elif [ -f .env.template ]; then
        cp .env.template .env
        sed -i "s/LOCAL_DB_USER=postgres/LOCAL_DB_USER=$PG_USER/g" .env
        sed -i "s/LOCAL_DB_PASSWORD=/LOCAL_DB_PASSWORD=$PG_PASS/g" .env
        print_success "Created .env file from template with PostgreSQL credentials"
    else
        print_warning ".env template not found. Please create manually."
    fi
}

# Setup conda environment
setup_conda() {
    print_section "Setting up conda environment"
    
    # Check if conda is installed
    if ! command_exists conda; then
        print_warning "Conda not found. Installing Miniconda..."
        
        # Download and install Miniconda
        wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
        bash miniconda.sh -b -p $HOME/miniconda3
        
        # Initialize conda
        eval "$($HOME/miniconda3/bin/conda shell.bash hook)"
        conda init bash
        
        print_success "Miniconda installed"
    else
        print_success "Conda is already installed"
    fi
    
    # Ensure conda is in PATH for this session
    if ! command_exists conda; then
        export PATH="$HOME/miniconda3/bin:$PATH"
    fi
    
    # Create conda environment if it doesn't exist
    if ! conda env list | grep -q "cryoprotect"; then
        print_warning "Creating conda environment..."
        
        if [ -f environment.yml ]; then
            conda env create -f environment.yml
        else
            conda create -n cryoprotect python=3.10 -y
            
            # Activate the environment
            conda activate cryoprotect
            
            # Install RDKit
            conda install -c conda-forge rdkit=2024.03.4 -y
            
            # Install other dependencies
            if [ -f requirements.txt ]; then
                pip install -r requirements.txt
            fi
        fi
        
        print_success "Conda environment created"
    else
        print_success "Conda environment already exists"
    fi
    
    # Activate the environment and verify installation
    eval "$(conda shell.bash hook)"
    conda activate cryoprotect
    
    # Verify RDKit installation
    python -c "from rdkit import Chem; print('RDKit installation successful!')" || {
        print_error "RDKit verification failed"
        exit 1
    }
    
    print_success "Conda environment setup complete"
}

# Setup Docker (optional)
setup_docker() {
    print_section "Setting up Docker (Optional)"
    
    read -p "Do you want to install Docker? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        # Install Docker
        if ! command_exists docker; then
            sudo dnf -y install dnf-plugins-core
            sudo dnf config-manager --add-repo https://download.docker.com/linux/fedora/docker-ce.repo
            sudo dnf install -y docker-ce docker-ce-cli containerd.io docker-compose-plugin
            
            # Start and enable Docker service
            sudo systemctl enable docker
            sudo systemctl start docker
            
            # Add current user to docker group
            sudo usermod -aG docker $(whoami)
            print_warning "You need to log out and log back in for docker group changes to take effect"
        else
            print_success "Docker is already installed"
        fi
        
        # Install Docker Compose if not already installed
        if ! command_exists docker-compose; then
            # Check if we're using docker-compose plugin
            if ! command_exists docker && docker compose version >/dev/null 2>&1; then
                print_success "Docker Compose plugin is already installed"
            else
                # Install standalone docker-compose
                DOCKER_COMPOSE_VERSION=$(curl -s https://api.github.com/repos/docker/compose/releases/latest | grep 'tag_name' | cut -d '"' -f 4)
                sudo curl -L "https://github.com/docker/compose/releases/download/${DOCKER_COMPOSE_VERSION}/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
                sudo chmod +x /usr/local/bin/docker-compose
            fi
        fi
        print_success "Docker setup complete"
    else
        print_warning "Skipping Docker installation"
    fi
}

# Setup systemd service
setup_systemd() {
    print_section "Setting up systemd service (Optional)"
    
    read -p "Do you want to set up a systemd service for CryoProtect? (y/n): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        # Create systemd service file
        echo "[Unit]
Description=CryoProtect Flask Application
After=network.target postgresql.service

[Service]
User=$(whoami)
WorkingDirectory=$(pwd)
Environment=PATH=$HOME/miniconda3/bin:$PATH
ExecStart=$HOME/miniconda3/envs/cryoprotect/bin/python app.py
Restart=always
RestartSec=10
StandardOutput=journal
StandardError=journal
SyslogIdentifier=cryoprotect
Environment=FLASK_ENV=production

[Install]
WantedBy=multi-user.target" > cryoprotect.service
        
        # Install service
        sudo mv cryoprotect.service /etc/systemd/system/
        sudo systemctl daemon-reload
        sudo systemctl enable cryoprotect.service
        
        print_success "Systemd service installed"
        print_warning "You can start the service with: sudo systemctl start cryoprotect.service"
    else
        print_warning "Skipping systemd service setup"
    fi
}

# Apply database migrations
apply_migrations() {
    print_section "Applying database migrations"
    
    # Check if the migrations directory exists
    if [ -d migrations ]; then
        # Check if we have node.js installed for migrations
        if command_exists node && [ -f migrations/apply_migration.js ]; then
            node migrations/apply_migration.js
            print_success "Database migrations applied using Node.js"
        # Check if we have Python migration scripts
        elif [ -f database/init_local_db.py ]; then
            python database/init_local_db.py
            print_success "Database initialized using Python script"
        else
            print_warning "No migration scripts found. Skipping database migrations."
        fi
    else
        print_warning "Migrations directory not found. Skipping database migrations."
    fi
}

# Final verification
verify_installation() {
    print_section "Verifying installation"
    
    # Check if conda environment is active
    if [ -z "$CONDA_PREFIX" ] || [[ "$CONDA_PREFIX" != *"cryoprotect"* ]]; then
        print_warning "Conda environment is not active. Activating..."
        eval "$(conda shell.bash hook)"
        conda activate cryoprotect
    fi
    
    # Check if Python is working
    python --version || {
        print_error "Python is not working properly"
        exit 1
    }
    
    # Check if RDKit is working
    python -c "from rdkit import Chem; print('RDKit is working')" || {
        print_error "RDKit is not working properly"
        exit 1
    }
    
    # Check if PostgreSQL is running
    if systemctl is-active --quiet postgresql; then
        print_success "PostgreSQL is running"
    else
        print_error "PostgreSQL is not running"
        exit 1
    }
    
    # Check database connection
    if command_exists psql; then
        PGPASSWORD=$PG_PASS psql -h localhost -U $PG_USER -d cryoprotect -c "SELECT version();" > /dev/null 2>&1 && {
            print_success "Database connection successful"
        } || {
            print_error "Database connection failed"
        }
    fi
    
    # Check if application runs
    if [ -f app.py ]; then
        print_warning "Testing application startup (will stop after 5 seconds)..."
        timeout 5 python app.py > /dev/null 2>&1 && {
            print_success "Application starts correctly"
        } || {
            print_warning "Application startup may have issues. This is not always a problem."
        }
    fi
    
    print_success "Verification complete"
}

# Main function
main() {
    print_section "CryoProtect Fedora Setup Script"
    
    # Check if running on Fedora
    check_fedora_version
    
    # Install dependencies
    install_dependencies
    
    # Setup PostgreSQL
    setup_postgresql
    
    # Setup conda environment
    setup_conda
    
    # Setup SELinux
    setup_selinux
    
    # Setup firewall
    setup_firewall
    
    # Setup Docker (optional)
    setup_docker
    
    # Apply database migrations
    apply_migrations
    
    # Setup systemd service (optional)
    setup_systemd
    
    # Verify installation
    verify_installation
    
    print_section "Setup Complete!"
    echo "You can now run the application with:"
    echo "  conda activate cryoprotect"
    echo "  python app.py"
    echo ""
    echo "Or use the provided scripts:"
    echo "  ./run_app.sh or ./run_app_with_fix.sh"
    echo ""
    echo "For more information, see the documentation in the docs/ directory."
}

# Run the main function
main