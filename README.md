# CryoProtect Analyzer

A tool for analyzing cryoprotectant molecules using RDKit and Supabase.

## Overview

CryoProtect Analyzer is a Flask-based web application that allows users to analyze cryoprotectant molecules, calculate their properties, and store the results in a Supabase database. The application uses RDKit for molecular property calculations and visualization.

## Features

- Molecular property calculation using RDKit
- Visualization of molecules
- Substructure search
- Similarity calculation
- Database storage using Supabase
- Web interface for easy access

## Installation

### Prerequisites

- Python 3.9+
- Conda
- Docker (optional)

### Setup

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/cryoprotect-analyzer.git
   cd cryoprotect-analyzer
   ```

2. Set up the environment:
   ```
   ./setup_environment.sh  # Linux/Mac
   setup_environment.bat   # Windows
   ```

3. Configure Supabase:
   - Create a `.env` file with your Supabase credentials
   - Apply the database migration using `node migrations/apply_migration.js`

4. Run the application:
   ```
   ./run_app.sh  # Linux/Mac
   run_app.bat   # Windows
   ```

### Docker

You can also run the application using Docker:

```
docker-compose up
```

## Usage

1. Open your browser and navigate to `http://localhost:5000`
2. Use the web interface to analyze molecules, create mixtures, and view predictions

## Documentation

- [API Documentation](README_API.md)
- [RDKit Integration](README_RDKit.md)
- [Supabase Integration](README_Supabase.md)
- [Web Interface](README_Web.md)
- [RDKit Troubleshooting](README_RDKit_Troubleshooting.md)

## License

MIT