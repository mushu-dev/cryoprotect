# CryoProtect Analyzer System Documentation

## System Overview

CryoProtect Analyzer is a Flask-based web application for analyzing cryoprotectant molecules, calculating their properties, and storing results in a Supabase database. The system allows researchers to analyze, compare, and design cryoprotective mixtures for biological sample preservation.

## Core Objectives

1. **Molecular Data Management:** Catalog and store molecular data from multiple sources (PubChem, ChEMBL)
2. **Property Calculation:** Calculate molecular properties using RDKit and custom algorithms
3. **Mixture Analysis:** Create, analyze, and optimize cryoprotectant mixtures
4. **Predictive Modeling:** Predict properties of novel mixtures
5. **Experiment Tracking:** Document and analyze experimental results
6. **Scientific Collaboration:** Support team-based scientific workflows

## Technology Stack

### Backend
- **Python 3.10+**: Primary programming language
- **Flask 3.0.2**: Web framework
- **Flask-RESTful 0.3.10**: REST API framework
- **Flask-CORS 4.0.0**: Cross-Origin Resource Sharing support
- **Marshmallow 3.21.0**: Object serialization/deserialization
- **Supabase 2.1.0**: PostgreSQL database with serverless API
- **RDKit 2024.03.4**: Cheminformatics toolkit for molecular analysis
- **Psycopg2**: Direct PostgreSQL connection
- **JWT**: Authentication with JSON Web Tokens
- **Prometheus**: Monitoring metrics collection

### Scientific Computing
- **NumPy/SciPy**: Numerical computations
- **Pandas**: Data analysis
- **Scikit-Learn**: Machine learning models
- **Matplotlib/Seaborn**: Data visualization

### Infrastructure
- **Podman**: Container management for Fedora Linux
- **Prometheus/Grafana**: Monitoring and alerting
- **ELK Stack**: Centralized logging
- **Scheduled Backups**: Automated backup system with integrity verification

### Security
- **Row-Level Security (RLS)**: Database access control
- **JWT-based Authentication**: Secure token-based authentication
- **Role-Based Access Control (RBAC)**: Permission-based authorization
- **API Rate Limiting**: Protection against abuse
- **SELinux Integration**: Enhanced security on Fedora

## Architecture

CryoProtect follows a layered architecture with clear separation of concerns:

### 1. Presentation Layer
- **REST API**: Primary interface for frontend and external clients
- **OpenAPI Documentation**: Self-documenting API using apispec
- **Rate Limiting**: Traffic control for API stability

### 2. Application Layer
- **Resource Controllers**: Implement API endpoints
- **Middleware**: Common functionality (authentication, logging, monitoring)
- **Service Modules**: Business logic implementation

### 3. Domain Layer
- **Models**: Domain entities (molecules, mixtures, experiments)
- **Services**: Core business logic
- **Scientific Calculations**: RDKit integration and custom algorithms

### 4. Data Access Layer
- **Supabase Adapter**: Primary database access
- **Direct Connection**: Performance-optimized PostgreSQL access
- **Cache**: Multi-level caching for improved performance

### 5. External Integration Layer
- **PubChem Client**: Resilient client for PubChem API
- **ChEMBL Client**: Resilient client for ChEMBL API
- **RDKit Bridge**: Integration with RDKit calculation engine

## Authentication and Security

### Authentication System
- JWT-based authentication with token refresh
- Role-based access control (RBAC)
- Service role for system operations
- Session management with secure cookies

### Security Features
- HTTPS-only cookies
- CSRF protection
- HTTP-only JWT cookies
- Content Security Policy (CSP)
- Fine-grained Row Level Security (RLS) policies
- SELinux context management for sensitive operations

## Key Components

### 1. Molecular Database
- Sources molecules from PubChem and ChEMBL
- Stores comprehensive molecular structures and properties
- Maintains versioning and modification history

### 2. Mixture Analysis System
- Allows creation of molecule mixtures
- Calculates predicted properties
- Compares mixtures based on properties

### 3. Scientific Calculation Engine
- Integrates with RDKit for molecular property calculations
- Implements custom algorithms for cryoprotection analysis
- Supports batch processing for large datasets

### 4. Experiment Tracking
- Records experimental conditions and results
- Compares experimental data with predictions
- Supports laboratory verification workflows

### 5. Monitoring and Observability
- Prometheus metrics for performance tracking
- Structured logging with ELK integration
- Health checks and alerting

### 6. Data Protection
- Automated backup system with retention policies
- Verification of backup integrity
- Restoration procedures

## Molecular Filtering and Search System

The system provides advanced molecular filtering capabilities critical for cryoprotectant research:

### Structure-Based Similarity Search

Multiple complementary methods are implemented to find structurally similar molecules:

1. **Morgan/ECFP Fingerprints**
   - Extended-Connectivity Fingerprints capture local molecular environments
   - Configurable radius (ECFP4, ECFP6) for different levels of structural detail
   - Fast comparison using Tanimoto similarity coefficient
   - Example API: `GET /api/molecules/{id}/similar?method=morgan&radius=2&limit=10`

2. **Topological Fingerprints**
   - Path-based structural features capturing molecular scaffolds
   - Effective for detecting molecules with similar backbones
   - Example API: `GET /api/molecules/{id}/similar?method=topological&limit=10`

3. **MACCS Keys**
   - 166 predefined structural features for focused similarity
   - Specialized for functional group matching
   - Example API: `GET /api/molecules/{id}/similar?method=maccs&limit=10`

4. **Pharmacophore-Based Search**
   - 3D arrangement of functional groups (hydrogen bond donors/acceptors, etc.)
   - More biologically relevant for cryoprotection mechanisms
   - Example API: `GET /api/molecules/pharmacophore?query=jsonConfig&limit=10`

### Property-Based Filtering

Molecules can be filtered based on scientific criteria relevant for cryoprotection:

1. **Physicochemical Properties**
   - Molecular weight (optimized range: 76-150 Da for cellular penetration)
   - LogP (water solubility indicator, preferred: -3 to 1)
   - Hydrogen bond donors (≥2) and acceptors (≥3) for water interaction
   - Polar surface area (optimal range: 40-90 Å²)
   - Rotatable bonds (<5 for conformational stability)

2. **Cryoprotection-Specific Criteria**
   - Glass transition temperature: Ability to form vitreous state
   - Ice nucleation inhibition: Ability to prevent ice crystal formation
   - Cell membrane permeability: Transport across biological membranes
   - Protein stabilization potential: Interaction with biomolecules
   - Cytotoxicity profiles: Safety for biological systems
   - Viscosity at low temperatures: Flow characteristics during freezing

3. **Implementation:**
   - Flexible JSON-based filter configuration
   - Cascading filter operations for multi-parameter refinement
   - Dynamic property calculation for unmeasured attributes
   - Batch processing for high-throughput screening

Example API request for property-based filtering:
```json
POST /api/molecules/filter
{
  "filters": [
    {"property": "molecular_weight", "min": 76, "max": 150},
    {"property": "h_bond_acceptors", "min": 3},
    {"property": "h_bond_donors", "min": 2},
    {"property": "logp", "min": -3, "max": 1},
    {"property": "glass_transition_temp", "min": -30}
  ],
  "sort_by": "predicted_cryoprotection_score",
  "limit": 50
}
```

### Combined Search Approaches

The system enables sophisticated search strategies combining multiple methods:

1. **Similarity-Based Filtering:**
   - Find molecules similar to known cryoprotectants, then apply property filters
   - Example: Similar to glycerol, but with improved glass transition temperature

2. **Property-First Screening:**
   - Apply rigorous property filters, then cluster by structural similarity
   - Useful for identifying novel chemical scaffolds with desired properties

3. **Machine Learning Enhanced Search:**
   - Train models on known cryoprotectants to predict effectiveness
   - Apply model-based scoring to filter results by predicted performance

4. **Experiment-Guided Refinement:**
   - Use experimental results to adjust property weights
   - Iteratively refine search parameters based on laboratory feedback

## API Structure

The API is organized into several logical groups:

1. **Molecules**: Endpoints for molecular data management
2. **Properties**: Molecular property calculation and retrieval
3. **Mixtures**: Mixture creation and analysis
4. **Experiments**: Experimental data management
5. **Predictions**: Property prediction endpoints
6. **Projects**: Project management endpoints
7. **Teams**: Team collaboration endpoints
8. **Users**: User profile and settings
9. **Analysis**: Data analysis and visualization

## Data Flow

1. **Data Import**: 
   - Molecules are imported from PubChem/ChEMBL via resilient clients
   - Properties are calculated using RDKit and stored

2. **User Workflow**:
   - Users create and manage projects and teams
   - Users search and select molecules for analysis
   - Users create mixtures of molecules
   - System predicts properties of mixtures
   - Users record experimental results
   - System compares predictions with experiments

3. **Analysis Pipeline**:
   - Raw molecular data → Property calculation → Mixture analysis → Prediction models → Visualization

## Performance Optimizations

1. **Connection Pooling**: Optimized database connection management
2. **Multi-level Caching**: In-memory and disk-based caching
3. **Batch Processing**: Efficient bulk operations
4. **Query Optimization**: Indexed searches and optimized queries
5. **Asynchronous Processing**: Background jobs for long-running operations

## Fedora-Specific Environment

### System Requirements
- Fedora 42+ (Workstation or Server edition)
- At least 4GB RAM, 20GB free disk space
- Python 3.10+ (included with Fedora)
- PostgreSQL 15+ (available via DNF)

### Setup Process
- Installation script: `fedora_setup.sh`
- Automated dependency installation
- SELinux configuration for Flask application
- Firewall configuration for required ports
- Systemd service creation (optional)

### Container Management
CryoProtect is designed for deployment using Podman containers, the default container engine in Fedora Linux:

1. **Development**: Local environment with hot reloading
   ```bash
   podman-compose -f podman-compose.dev.yml up
   ```

2. **Staging**: Testing environment with staging database
   ```bash
   podman-compose up
   ```

3. **Production**: Optimized for reliability and performance
   ```bash
   podman-compose -f podman-compose.prod.yml up -d
   ```

### Network Connectivity Configuration

CryoProtect's connection to Supabase requires specific network considerations, especially on Fedora systems that prefer IPv6 by default:

#### IPv4/IPv6 Configuration

1. **Supabase Connection Configuration**:
   - In `config.py`, explicitly force IPv4 for Supabase connections:
   ```python
   SUPABASE_OPTIONS = {
       "postgrest_client_options": {
           "fetch_parameters": {
               "family": 4  # Force IPv4
           }
       }
   }
   ```

2. **Connection Testing**:
   - Use the included `scripts/test_connectivity.sh` script to verify proper IPv4 connectivity
   - Diagnose and troubleshoot dual-stack behavior

3. **System-level Configuration**:
   - If needed, modify `/etc/gai.conf` to prioritize IPv4 over IPv6:
   ```
   precedence ::ffff:0:0/96  100
   ```

4. **Fallback Mechanisms**:
   - Implement connection retry logic with protocol switching
   - Log connection attempts for troubleshooting

### Podman-Docker Compatibility
- The project maintains Docker Compose files that are compatible with Podman
- `podman-compose` is used instead of `docker-compose`
- Container images are built and managed with Podman
- Root-less container execution for enhanced security

## Implementation Priorities

For efficient implementation, we'll focus on the following priority areas:

1. **Core Database Functionality**: Ensure proper connection and data flow
2. **API Endpoints for Critical Features**: Focus on molecules, properties, and mixtures first
3. **Authentication and Security**: Implement JWT and RLS policies
4. **Basic UI**: Simple interface for core operations
5. **Testing Framework**: Comprehensive test coverage

## Project Phases

The project is currently divided into four phases:

1. **Phase 1: Technical Foundation** (Completed)
   - Database Architecture
   - Authentication System

2. **Phase 2: Feature Completion** (Current)
   - API Layer Completion
   - Core Functionality
   - User Interface

3. **Phase 3: Production Readiness** (Planned)
   - Deployment Infrastructure with Podman
   - Monitoring and Maintenance
   - Security with SELinux integration

4. **Phase 4: Documentation and Knowledge Transfer** (Planned)
   - Documentation
   - Knowledge Transfer

## Known Limitations and Gaps

1. **SELinux Integration**: Need to create proper security contexts for application
2. **Podman Configuration**: Need to adapt Docker files for Podman
3. **Fedora-specific Paths**: Some path references need updates for Fedora filesystem layout
4. **Systemd Service Integration**: Create service files for automatic startup

## Future Directions

1. **Enhanced Machine Learning Models**: More advanced predictive algorithms
2. **Toxicity Prediction**: Improved toxicity assessment
3. **Interactive Molecular Visualizations**: 3D visualization of molecular structures
4. **Integration with Lab Equipment**: Direct data collection from lab instruments
5. **Expanded API Ecosystem**: Third-party integration support