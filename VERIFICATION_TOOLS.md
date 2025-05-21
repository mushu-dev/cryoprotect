# CryoProtect Verification Tools

This document describes the verification tools available to test your CryoProtect setup on Fedora Linux.

## Quick Start

To verify your Fedora environment and CryoProtect setup:

```bash
# Run the main verification script
./verify_fedora_setup.sh
```

This script will check all required components and display a summary of your setup status.

## Available Tools

### 1. verify_fedora_setup.sh

The main verification script that checks:
- Fedora version
- SELinux status
- Python installation
- Environment configuration
- Supabase connectivity
- RDKit installation
- Simplified app availability

### 2. check_rdkit.py

A simple Python script to verify RDKit installation and version:

```bash
# Activate your environment first
source quick_env/bin/activate

# Run the RDKit check
python check_rdkit.py

# Exit code 0 means RDKit is installed
# Exit code 1 means RDKit is not installed
# Exit code 2 means an error occurred during checking
```

### 3. simplified_app.py

A minimal Flask application that demonstrates core functionality:
- Basic Flask server
- Health check endpoint
- Supabase connectivity test
- Dependencies verification
- Sample data endpoint

### 4. run_simplified_app.sh

A script to run the simplified app with proper environment setup:

```bash
# Run in foreground
./run_simplified_app.sh

# Run in background
./run_simplified_app_background.sh
```

### 5. DATABASE_SCHEMA_SUMMARY.md

Documentation of the database schema with:
- Table descriptions
- Record counts
- Key relationships
- Sample queries

## Troubleshooting

### Common Issues

1. **SELinux Blocks Connections**
   - Check SELinux denials: `sudo ausearch -m avc --start recent`
   - Create custom SELinux policy if needed

2. **RDKit Not Found**
   - Ensure conda environment is activated
   - Install with: `conda install -c conda-forge rdkit`

3. **Supabase Connection Fails**
   - Verify .env file has correct credentials
   - Check network connectivity with: `curl https://your-project.supabase.co/rest/v1/`
   - Ensure firewalld allows outgoing connections

## Next Steps

After successful verification:
1. Run the full CryoProtect application
2. Set up development environment
3. Configure SELinux contexts for production
4. Implement Podman containerization

See [FEDORA_SETUP_SUMMARY.md](FEDORA_SETUP_SUMMARY.md) for more details on next steps.