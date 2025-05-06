#!/bin/bash
# CryoProtect - Complete Fresh Linux Installation Script
# This script runs all necessary steps to get CryoProtect running on a fresh Linux installation

set -e  # Exit on error

echo "================================================"
echo "CryoProtect Fresh Linux Installation"
echo "================================================"

# Ask user if they want to use Docker or native installation
echo "How would you like to install CryoProtect?"
echo "1) Native installation (recommended for development)"
echo "2) Docker installation (easier, but less flexible)"
read -p "Enter your choice (1/2): " INSTALL_CHOICE

case $INSTALL_CHOICE in
    1)
        echo -e "\nProceeding with native installation..."
        
        # Run the dependency installation script
        echo -e "\n[1/4] Installing Linux dependencies..."
        ./install_linux_dependencies.sh
        
        # Set up environment after reloading bash profile
        echo -e "\n[2/4] Activating environment..."
        source ~/.bashrc
        eval "$(conda shell.bash hook)"
        conda activate cryoprotect
        
        # Apply database migrations
        echo -e "\n[3/4] Applying database migrations..."
        if [ -f "migrations/apply_migration.js" ]; then
            node migrations/apply_migration.js
        else
            echo "⚠️ Migration script not found. Please apply migrations manually."
        fi
        
        # Verify the installation
        echo -e "\n[4/4] Verifying installation..."
        ./verify_linux_migration.sh
        
        echo -e "\n================================================"
        echo "Native Installation Complete!"
        echo "================================================"
        
        echo -e "\nTo run CryoProtect:"
        echo "1. Activate the conda environment: conda activate cryoprotect"
        echo "2. Run the application: ./run_app_with_fix.sh"
        echo -e "\nAccess the application at http://localhost:5000"
        ;;
        
    2)
        echo -e "\nProceeding with Docker installation..."
        
        # Run the Docker setup script
        echo -e "\n[1/2] Setting up Docker..."
        ./docker_setup.sh
        
        # Build and start the containers
        echo -e "\n[2/2] Building and starting Docker containers..."
        docker-compose up -d
        
        echo -e "\n================================================"
        echo "Docker Installation Complete!"
        echo "================================================"
        
        echo -e "\nCryoProtect is now running via Docker!"
        echo "Access the application at http://localhost:5000"
        echo -e "\nTo view logs:"
        echo "   docker-compose logs -f"
        echo -e "\nTo stop the containers:"
        echo "   docker-compose down"
        ;;
        
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
esac

echo -e "\nThank you for installing CryoProtect!"
echo "For more information, see the documentation in the docs/ directory."
echo "If you encounter any issues, please check LINUX_MIGRATION_GUIDE.md"