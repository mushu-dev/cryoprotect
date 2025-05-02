# ToxCast/Tox21 Implementation Prompt for ROO Agent

You are responsible for implementing the remaining components of the toxicity data integration for CryoProtect v2. This system already has a significant foundation for Tox21 toxicity data integration, but several key components still need to be implemented.

## Understanding the Current Implementation

The project has already implemented:

1. **Database Schema** in `migrations/012_toxicity_schema.sql`
2. **Chemical Data Module** in `chemical_data/toxicity/` directory
3. **API Resources** in `api/toxicity_resources.py`
4. **Unified Scoring System** in `api/unified_scoring.py`
5. **Data Models** in `api/models.py`
6. **Migration Scripts** like `apply_toxicity_schema.py`

## Your Implementation Tasks

Your task is to complete the toxicity data implementation by focusing on these key areas:

### 1. Test Suite Development

Create comprehensive test files for the toxicity data integration:

```python
# tests/test_toxicity_data.py
"""Unit tests for the ToxicityData class and related functionality."""
import unittest
from unittest.mock import patch, MagicMock
from api.models import ToxicityData

class TestToxicityData(unittest.TestCase):
    """Test cases for ToxicityData class."""

    def setUp(self):
        """Set up test case."""
        # Your setup code here
        pass

    @patch('api.models.ToxicityData.get_supabase')
    def test_get_for_molecule(self, mock_get_supabase):
        """Test get_for_molecule method."""
        # Mock setup and assertions
        pass

    # Additional test methods
```

```python
# tests/test_toxicity_scorer.py
"""Unit tests for the ToxicityScorer class."""
import unittest
from unittest.mock import patch, MagicMock
from chemical_data.toxicity.toxicity_scorer import ToxicityScorer

class TestToxicityScorer(unittest.TestCase):
    """Test cases for ToxicityScorer class."""

    def setUp(self):
        """Set up test case."""
        # Your setup code here
        pass

    def test_calculate_endpoint_score(self):
        """Test _calculate_endpoint_score method."""
        # Setup and assertions
        pass

    # Additional test methods
```

```python
# tests/test_toxicity_resources.py
"""Integration tests for toxicity API resources."""
import unittest
from unittest.mock import patch, MagicMock
from app import app

class TestToxicityResources(unittest.TestCase):
    """Test cases for toxicity API resources."""

    def setUp(self):
        """Set up test case."""
        self.app = app.test_client()
        # Additional setup
        pass

    def test_get_toxicity_data(self):
        """Test GET /api/toxicity/molecule/{id} endpoint."""
        # Test implementation
        pass

    # Additional test methods
```

### 2. Data Population Script

Create a script to populate toxicity data:

```python
# populate_toxicity_data.py
#!/usr/bin/env python3
"""
Populate toxicity data from Tox21.

This script:
1. Downloads and processes Tox21 assay and chemical data
2. Maps compounds to existing molecules in the database
3. Imports toxicity data and calculates scores
"""

import os
import sys
import logging
import argparse
from supabase import create_client
from config import Config
from chemical_data.toxicity.tox21_client import Tox21Client
from chemical_data.toxicity.toxicity_scorer import ToxicityScorer

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Populate toxicity data from Tox21')
    parser.add_argument('--limit', type=int, default=None, help='Limit the number of compounds to process')
    parser.add_argument('--force-refresh', action='store_true', help='Force refresh of cached data')
    parser.add_argument('--dry-run', action='store_true', help='Perform a dry run without modifying the database')
    return parser.parse_args()

def main():
    """Main function to populate toxicity data."""
    args = parse_arguments()
    
    # Initialize Supabase client
    config = Config()
    supabase_url = os.environ.get("SUPABASE_URL") or config.SUPABASE_URL
    supabase_key = os.environ.get("SUPABASE_KEY") or config.SUPABASE_KEY
    
    if not supabase_url or not supabase_key:
        logger.error("Supabase URL or key not found")
        return False
    
    supabase = create_client(supabase_url, supabase_key)
    
    # Initialize Tox21 client and toxicity scorer
    tox21_client = Tox21Client(supabase)
    toxicity_scorer = ToxicityScorer(supabase)
    
    # Download and import Tox21 assay data
    logger.info("Downloading and importing Tox21 assay data")
    if not args.dry_run:
        assay_file = tox21_client.download_assay_data(force_refresh=args.force_refresh)
        assays_imported = tox21_client.import_assays(assay_file)
        logger.info(f"Imported {assays_imported} Tox21 assays")
    else:
        logger.info("Dry run - would import Tox21 assay data")
    
    # Download and import Tox21 chemical data
    logger.info("Downloading and importing Tox21 chemical data")
    if not args.dry_run:
        chemical_file = tox21_client.download_chemical_data(force_refresh=args.force_refresh)
        stats = tox21_client.import_chemical_data(chemical_file, max_compounds=args.limit)
        logger.info(f"Chemical data import stats: {stats}")
    else:
        logger.info("Dry run - would import Tox21 chemical data")
    
    # Calculate toxicity scores
    logger.info("Calculating toxicity scores for molecules")
    if not args.dry_run:
        scores_calculated = toxicity_scorer.calculate_toxicity_scores()
        logger.info(f"Calculated toxicity scores for {scores_calculated} molecules")
    else:
        logger.info("Dry run - would calculate toxicity scores")
    
    logger.info("Toxicity data population completed successfully")
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
```

### 3. Visualization Components

Create JavaScript visualization for toxicity data:

```javascript
// static/js/toxicity-visualization.js
/**
 * Toxicity Visualization Module
 * 
 * This module provides visualization components for toxicity data,
 * including toxicity scores, endpoint breakdowns, and unified scoring.
 */

class ToxicityVisualization {
    /**
     * Initialize the toxicity visualization
     * @param {string} containerId - ID of the container element
     */
    constructor(containerId) {
        this.container = document.getElementById(containerId);
        if (!this.container) {
            console.error(`Container element with ID "${containerId}" not found`);
            return;
        }
        
        this.colors = {
            efficacy: '#4CAF50',
            toxicity: '#F44336',
            glassTransition: '#2196F3',
            overall: '#9C27B0'
        };
        
        this.initialize();
    }
    
    /**
     * Initialize the visualization components
     */
    initialize() {
        // Create container elements
        this.scoreContainer = document.createElement('div');
        this.scoreContainer.className = 'toxicity-score-container';
        this.container.appendChild(this.scoreContainer);
        
        this.endpointContainer = document.createElement('div');
        this.endpointContainer.className = 'toxicity-endpoint-container';
        this.container.appendChild(this.endpointContainer);
        
        this.unifiedContainer = document.createElement('div');
        this.unifiedContainer.className = 'unified-score-container';
        this.container.appendChild(this.unifiedContainer);
    }
    
    /**
     * Load and display toxicity data for a molecule
     * @param {string} moleculeId - ID of the molecule
     */
    loadMoleculeToxicity(moleculeId) {
        // Fetch toxicity data from API
        fetch(`/api/toxicity/molecule/${moleculeId}`)
            .then(response => response.json())
            .then(data => {
                this.renderToxicityData(data);
            })
            .catch(error => {
                console.error('Error loading toxicity data:', error);
                this.renderError('Failed to load toxicity data');
            });
        
        // Fetch toxicity scores
        fetch(`/api/toxicity/scores/molecule/${moleculeId}`)
            .then(response => response.json())
            .then(data => {
                this.renderToxicityScores(data);
            })
            .catch(error => {
                console.error('Error loading toxicity scores:', error);
            });
        
        // Fetch unified score
        fetch(`/api/toxicity/unified/molecule/${moleculeId}`)
            .then(response => response.json())
            .then(data => {
                this.renderUnifiedScore(data);
            })
            .catch(error => {
                console.error('Error loading unified score:', error);
            });
    }
    
    /**
     * Render toxicity data
     * @param {Object} data - Toxicity data
     */
    renderToxicityData(data) {
        // Implementation details
    }
    
    /**
     * Render toxicity scores
     * @param {Object} data - Toxicity score data
     */
    renderToxicityScores(data) {
        // Implementation details
    }
    
    /**
     * Render unified score
     * @param {Object} data - Unified score data
     */
    renderUnifiedScore(data) {
        // Implementation details
    }
    
    /**
     * Render error message
     * @param {string} message - Error message
     */
    renderError(message) {
        this.container.innerHTML = `<div class="error-message">${message}</div>`;
    }
}

// Initialize when document is ready
document.addEventListener('DOMContentLoaded', () => {
    // Check if we're on a molecule details page
    const moleculeIdElement = document.getElementById('molecule-id');
    if (moleculeIdElement) {
        const moleculeId = moleculeIdElement.value;
        const toxicityViz = new ToxicityVisualization('toxicity-viz-container');
        toxicityViz.loadMoleculeToxicity(moleculeId);
    }
});
```

### 4. Documentation

Update API documentation:

```markdown
# Toxicity Data API Documentation

## Overview

The Toxicity Data API provides access to toxicity data from the Tox21 program. This data includes assay results, toxicity scores, and unified scores that combine toxicity with efficacy and glass transition temperature (Tg) data.

## Endpoints

### Get Toxicity Data for a Molecule

`GET /api/toxicity/molecule/{molecule_id}`

Retrieves toxicity data for a specific molecule.

**Parameters:**
- `molecule_id` (path parameter): UUID of the molecule
- `assay_id` (query parameter, optional): Filter by assay ID
- `endpoint` (query parameter, optional): Filter by toxicological endpoint
- `active_only` (query parameter, optional): Return only active results (default: false)

**Response:**
```json
{
  "molecule_id": "123e4567-e89b-12d3-a456-426614174000",
  "toxicity_data": [
    {
      "id": "789e0123-c45d-67e8-f901-234567890123",
      "assay_id": "456e7890-a12b-34c5-d67e-890123456789",
      "activity_value": 10.5,
      "activity_unit": "µM",
      "hit_call": true,
      "reliability_score": 0.85
    }
  ],
  "count": 1
}
```

### Get Toxicity Endpoints

`GET /api/toxicity/endpoints`

Retrieves available toxicity endpoints.

**Response:**
```json
{
  "endpoints": [
    {
      "id": "123e4567-e89b-12d3-a456-426614174000",
      "name": "nuclear_receptor",
      "description": "Nuclear receptor activity",
      "category": "endocrine"
    }
  ],
  "count": 1
}
```

### [Additional endpoint documentation...]
```

### 5. Performance Optimization

```python
# optimize_toxicity_queries.py
#!/usr/bin/env python3
"""
Optimize database queries for toxicity data.

This script:
1. Adds additional indexes for commonly queried toxicity data
2. Creates materialized views for frequently accessed data
3. Verifies query performance with explain analyze
"""

import os
import sys
import logging
import time
from supabase import create_client
from config import Config

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def add_indexes(supabase):
    """Add performance indexes for toxicity data."""
    logger.info("Adding performance indexes for toxicity data...")
    
    indexes = [
        # Compound indexes for commonly used query patterns
        "CREATE INDEX IF NOT EXISTS idx_toxicity_data_molecule_hit_call ON toxicity_data (molecule_id, hit_call)",
        "CREATE INDEX IF NOT EXISTS idx_toxicity_data_assay_activity ON toxicity_data (assay_id, activity_value)",
        "CREATE INDEX IF NOT EXISTS idx_toxicity_score_molecule_type ON toxicity_score (molecule_id, score_type)",
        "CREATE INDEX IF NOT EXISTS idx_assay_endpoint_relevance ON assay_endpoint_mapping (assay_id, endpoint_id, relevance_score)",
        "CREATE INDEX IF NOT EXISTS idx_toxicity_assay_endpoint ON toxicity_assay (toxicological_endpoint)"
    ]
    
    success_count = 0
    for index_sql in indexes:
        try:
            logger.info(f"Creating index: {index_sql}")
            # Execute the index creation
            supabase.postgrest.rpc('exec_sql', {'query': index_sql}).execute()
            success_count += 1
        except Exception as e:
            logger.error(f"Error creating index: {str(e)}")
    
    logger.info(f"Successfully created {success_count}/{len(indexes)} indexes")
    return success_count == len(indexes)

def create_materialized_views(supabase):
    """Create materialized views for frequently accessed toxicity data."""
    logger.info("Creating materialized views for toxicity data...")
    
    views = [
        """
        CREATE MATERIALIZED VIEW IF NOT EXISTS toxicity_molecule_summary AS
        SELECT 
            m.id AS molecule_id,
            m.name AS molecule_name,
            m.toxicity_score,
            COUNT(td.id) AS total_assays,
            SUM(CASE WHEN td.hit_call = true THEN 1 ELSE 0 END) AS active_assays,
            AVG(td.reliability_score) AS avg_reliability
        FROM 
            molecule m
        LEFT JOIN 
            toxicity_data td ON m.id = td.molecule_id
        GROUP BY 
            m.id, m.name, m.toxicity_score
        """,
        
        """
        CREATE MATERIALIZED VIEW IF NOT EXISTS toxicity_endpoint_summary AS
        SELECT 
            te.id AS endpoint_id,
            te.name AS endpoint_name,
            te.category,
            COUNT(DISTINCT ta.id) AS assay_count,
            COUNT(DISTINCT td.molecule_id) AS molecule_count
        FROM 
            toxicity_endpoint te
        LEFT JOIN 
            assay_endpoint_mapping aem ON te.id = aem.endpoint_id
        LEFT JOIN 
            toxicity_assay ta ON aem.assay_id = ta.id
        LEFT JOIN 
            toxicity_data td ON ta.id = td.assay_id
        GROUP BY 
            te.id, te.name, te.category
        """
    ]
    
    success_count = 0
    for view_sql in views:
        try:
            logger.info(f"Creating materialized view: {view_sql.split('AS')[0]}")
            # Execute the view creation
            supabase.postgrest.rpc('exec_sql', {'query': view_sql}).execute()
            success_count += 1
        except Exception as e:
            logger.error(f"Error creating materialized view: {str(e)}")
    
    logger.info(f"Successfully created {success_count}/{len(views)} materialized views")
    return success_count == len(views)

def test_query_performance(supabase):
    """Test the performance of common toxicity data queries."""
    logger.info("Testing query performance...")
    
    queries = [
        # Common query patterns
        "SELECT * FROM toxicity_data WHERE molecule_id = '00000000-0000-0000-0000-000000000000' LIMIT 1",
        "SELECT * FROM toxicity_score WHERE molecule_id = '00000000-0000-0000-0000-000000000000' LIMIT 1",
        "SELECT * FROM toxicity_assay WHERE toxicological_endpoint = 'nuclear_receptor' LIMIT 10"
    ]
    
    # Replace the dummy UUID with a real molecule ID
    try:
        molecule_response = supabase.table("molecule").select("id").limit(1).execute()
        if molecule_response.data and len(molecule_response.data) > 0:
            real_id = molecule_response.data[0]["id"]
            queries = [q.replace('00000000-0000-0000-0000-000000000000', real_id) for q in queries]
    except Exception as e:
        logger.warning(f"Couldn't get real molecule ID: {str(e)}")
    
    for query in queries:
        try:
            logger.info(f"Testing query: {query}")
            # Execute EXPLAIN ANALYZE
            explain_query = f"EXPLAIN ANALYZE {query}"
            start_time = time.time()
            supabase.postgrest.rpc('exec_sql', {'query': explain_query}).execute()
            duration = time.time() - start_time
            logger.info(f"Query executed in {duration:.4f} seconds")
        except Exception as e:
            logger.error(f"Error testing query: {str(e)}")
    
    return True

def main():
    """Main function to optimize toxicity queries."""
    # Initialize Supabase client
    config = Config()
    supabase_url = os.environ.get("SUPABASE_URL") or config.SUPABASE_URL
    supabase_key = os.environ.get("SUPABASE_KEY") or config.SUPABASE_KEY
    
    if not supabase_url or not supabase_key:
        logger.error("Supabase URL or key not found")
        return False
    
    supabase = create_client(supabase_url, supabase_key)
    
    # Add indexes
    if not add_indexes(supabase):
        logger.warning("Failed to add all indexes")
    
    # Create materialized views
    if not create_materialized_views(supabase):
        logger.warning("Failed to create all materialized views")
    
    # Test query performance
    test_query_performance(supabase)
    
    logger.info("Toxicity query optimization completed")
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
```

## Implementation Guidelines

1. **Follow existing code patterns:** Maintain consistency with the existing codebase.
2. **Add comprehensive tests:** Ensure at least 70% test coverage for all new code.
3. **Document thoroughly:** Include docstrings in code and update API documentation.
4. **Focus on performance:** Optimize database queries and batch operations.
5. **Ensure security:** Follow the existing RLS policies pattern for new features.

## Success Criteria

Your implementation will be considered complete when:

1. All toxicity-related tests pass with >70% code coverage
2. Tox21 data can be successfully imported for at least 1000 compounds
3. Toxicity scores are calculated for all molecules in the database
4. Visualization components effectively display the toxicity data
5. Documentation covers all aspects of the toxicity integration
6. Query performance meets requirements (API response <200ms)

## Resources

1. Existing code in the `chemical_data/toxicity/` directory
2. Database schema in `migrations/012_toxicity_schema.sql`
3. API resources in `api/toxicity_resources.py`
4. Models in `api/models.py`
5. Unified scoring in `api/unified_scoring.py`

I've provided implementation snippets as guidance, but you should adapt these to match the specific requirements and patterns of the CryoProtect v2 codebase.