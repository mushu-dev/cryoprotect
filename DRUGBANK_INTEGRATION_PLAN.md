# DrugBank Integration Implementation Plan

## Overview

This document outlines the technical implementation plan for integrating DrugBank as a data source into the CryoProtect v2 system. DrugBank will provide rich pharmaceutical data that significantly enhances our ability to identify and evaluate cryoprotectant candidates by leveraging known drug properties, interactions, and pharmacological data.

## Business Value

1. **Enhanced Compound Screening**: Access to pharmaceutically-validated compounds with known safety profiles
2. **Rich Property Data**: Detailed pharmacokinetic and pharmacodynamic data to improve ML model training
3. **Drug Repurposing**: Identify existing approved drugs that may have cryoprotectant properties
4. **Mechanism Insights**: Leverage known mechanism of action data to understand cryoprotection mechanisms
5. **Interaction Data**: Drug-target and drug-drug interaction data to improve mixture modeling

## Technical Implementation Steps

### 1. DrugBank Client Module (Week 1)

```python
# drugbank/client.py
from dataclasses import dataclass
import requests
import time
import logging
from typing import Dict, List, Optional, Any, Union
from xml.etree import ElementTree as ET

@dataclass
class DrugBankCredentials:
    """DrugBank API credentials"""
    api_key: str
    app_id: str

class DrugBankClient:
    """Client for interacting with DrugBank API and data"""
    
    def __init__(self, credentials: DrugBankCredentials, 
                 rate_limit: float = 1.0,
                 cache_dir: str = "./cache/drugbank"):
        """Initialize the DrugBank client
        
        Args:
            credentials: DrugBank API credentials
            rate_limit: Maximum requests per second
            cache_dir: Directory for caching responses
        """
        self.credentials = credentials
        self.base_url = "https://go.drugbank.com/api/v1"
        self.rate_limit = rate_limit
        self.last_request_time = 0
        self.cache_dir = cache_dir
        
        # Create cache directory
        os.makedirs(cache_dir, exist_ok=True)
        
        # Set up logging
        self.logger = logging.getLogger(__name__)
    
    def _apply_rate_limit(self):
        """Apply rate limiting to API requests"""
        current_time = time.time()
        elapsed = current_time - self.last_request_time
        
        if elapsed < (1.0 / self.rate_limit):
            sleep_time = (1.0 / self.rate_limit) - elapsed
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()
    
    def _get_auth_headers(self):
        """Get authentication headers for API requests"""
        return {
            "Authorization": f"Bearer {self.credentials.api_key}",
            "X-DrugBank-App-ID": self.credentials.app_id
        }
    
    def get_drug_by_id(self, drugbank_id: str) -> Dict:
        """Get drug information by DrugBank ID
        
        Args:
            drugbank_id: DrugBank ID (e.g., 'DB00001')
            
        Returns:
            Dictionary containing drug information
        """
        self._apply_rate_limit()
        
        url = f"{self.base_url}/drugs/{drugbank_id}"
        response = requests.get(url, headers=self._get_auth_headers())
        
        if response.status_code == 200:
            return response.json()
        else:
            self.logger.error(f"Error fetching drug {drugbank_id}: {response.status_code}")
            response.raise_for_status()
    
    def search_drugs(self, query: str, page: int = 1, per_page: int = 10) -> Dict:
        """Search for drugs by name, synonym, or identifier
        
        Args:
            query: Search query
            page: Page number
            per_page: Results per page
            
        Returns:
            Dictionary containing search results
        """
        self._apply_rate_limit()
        
        url = f"{self.base_url}/drugs/search"
        params = {
            "q": query,
            "page": page,
            "per_page": per_page
        }
        
        response = requests.get(url, headers=self._get_auth_headers(), params=params)
        
        if response.status_code == 200:
            return response.json()
        else:
            self.logger.error(f"Error searching drugs: {response.status_code}")
            response.raise_for_status()
    
    def get_drug_interactions(self, drugbank_id: str) -> List[Dict]:
        """Get drug-drug interactions for a specific drug
        
        Args:
            drugbank_id: DrugBank ID
            
        Returns:
            List of drug interactions
        """
        self._apply_rate_limit()
        
        url = f"{self.base_url}/drugs/{drugbank_id}/interactions"
        response = requests.get(url, headers=self._get_auth_headers())
        
        if response.status_code == 200:
            return response.json()
        else:
            self.logger.error(f"Error fetching interactions for drug {drugbank_id}: {response.status_code}")
            response.raise_for_status()
    
    def extract_structure_info(self, drug_data: Dict) -> Dict:
        """Extract structure information from drug data
        
        Args:
            drug_data: Drug data from API
            
        Returns:
            Dictionary with structure information
        """
        structure_info = {
            "smiles": drug_data.get("smiles", None),
            "inchi": drug_data.get("inchi", None),
            "inchikey": drug_data.get("inchikey", None),
            "formula": drug_data.get("formula", None),
            "molecular_weight": drug_data.get("molecular_weight", None),
            "drugbank_id": drug_data.get("drugbank_id", None),
            "name": drug_data.get("name", None)
        }
        
        return structure_info
    
    def extract_properties(self, drug_data: Dict) -> List[Dict]:
        """Extract property information from drug data
        
        Args:
            drug_data: Drug data from API
            
        Returns:
            List of properties as dictionaries
        """
        properties = []
        
        # Basic properties
        if "molecular_weight" in drug_data:
            properties.append({
                "name": "Molecular Weight",
                "value": drug_data["molecular_weight"],
                "type": "numeric"
            })
        
        # Pharmacokinetic properties
        if "pharmacokinetics" in drug_data:
            pk = drug_data["pharmacokinetics"]
            
            for pk_property in ["absorption", "distribution", "metabolism", "elimination"]:
                if pk_property in pk:
                    properties.append({
                        "name": f"Pharmacokinetics - {pk_property.capitalize()}",
                        "value": pk[pk_property],
                        "type": "text"
                    })
        
        # ADMET properties
        if "admet" in drug_data:
            admet = drug_data["admet"]
            
            for admet_property, value in admet.items():
                properties.append({
                    "name": f"ADMET - {admet_property.capitalize()}",
                    "value": value,
                    "type": "text" if isinstance(value, str) else "numeric"
                })
        
        return properties
    
    def extract_categories(self, drug_data: Dict) -> List[str]:
        """Extract drug categories
        
        Args:
            drug_data: Drug data from API
            
        Returns:
            List of category names
        """
        if "categories" in drug_data and isinstance(drug_data["categories"], list):
            return [category.get("name") for category in drug_data["categories"] if "name" in category]
        return []
```

### 2. Database Schema Updates (Week 1)

```sql
-- Create tables for DrugBank integration

-- Drug categories table
CREATE TABLE IF NOT EXISTS drug_categories (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT NOT NULL UNIQUE,
    description TEXT,
    parent_category_id UUID REFERENCES drug_categories(id),
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- Drug-category mapping table
CREATE TABLE IF NOT EXISTS molecule_drug_categories (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id UUID NOT NULL REFERENCES molecule(id),
    category_id UUID NOT NULL REFERENCES drug_categories(id),
    source TEXT DEFAULT 'DrugBank',
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(molecule_id, category_id)
);

-- Drug targets table
CREATE TABLE IF NOT EXISTS drug_targets (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    name TEXT NOT NULL,
    uniprot_id TEXT,
    gene_name TEXT,
    organism TEXT,
    target_type TEXT,
    description TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- Drug-target mapping table
CREATE TABLE IF NOT EXISTS molecule_drug_targets (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id UUID NOT NULL REFERENCES molecule(id),
    target_id UUID NOT NULL REFERENCES drug_targets(id),
    action TEXT,
    affinity_value FLOAT,
    affinity_unit TEXT,
    source TEXT DEFAULT 'DrugBank',
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(molecule_id, target_id)
);

-- Drug interactions table
CREATE TABLE IF NOT EXISTS drug_interactions (
    id UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    molecule_id_1 UUID NOT NULL REFERENCES molecule(id),
    molecule_id_2 UUID NOT NULL REFERENCES molecule(id),
    description TEXT,
    severity TEXT,
    interaction_type TEXT,
    source TEXT DEFAULT 'DrugBank',
    created_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP,
    UNIQUE(molecule_id_1, molecule_id_2)
);

-- Add fields to molecule table for DrugBank information
ALTER TABLE molecule 
ADD COLUMN IF NOT EXISTS drugbank_id TEXT,
ADD COLUMN IF NOT EXISTS drug_type TEXT,
ADD COLUMN IF NOT EXISTS approved BOOLEAN DEFAULT FALSE;

-- Add new property types
INSERT INTO property_types (name, data_type, description, units)
VALUES 
('LogD', 'numeric', 'Distribution coefficient at pH 7.4', NULL),
('pKa', 'numeric', 'Acid dissociation constant', NULL),
('Bioavailability', 'numeric', 'Oral bioavailability percentage', '%'),
('Half Life', 'numeric', 'Elimination half-life', 'hours'),
('Drug Class', 'text', 'Pharmacological classification', NULL),
('Mechanism of Action', 'text', 'Primary mechanism of action', NULL),
('Blood-Brain Barrier', 'boolean', 'Ability to cross the blood-brain barrier', NULL),
('Protein Binding', 'numeric', 'Percentage bound to plasma proteins', '%'),
('Cmax', 'numeric', 'Maximum serum concentration', 'ng/mL'),
('AUC', 'numeric', 'Area under curve', 'ng·h/mL'),
('Volume of Distribution', 'numeric', 'Apparent volume of distribution', 'L/kg'),
('Clearance', 'numeric', 'Systemic clearance', 'mL/min/kg'),
('Cryoprotection Evidence', 'text', 'Evidence of cryoprotection activity', NULL),
('Vitrification Potential', 'numeric', 'Potential for vitrification', NULL),
('Cell Preservation Score', 'numeric', 'Effectiveness in cell preservation', NULL);

-- Add RLS policies
ALTER TABLE drug_categories ENABLE ROW LEVEL SECURITY;
ALTER TABLE molecule_drug_categories ENABLE ROW LEVEL SECURITY;
ALTER TABLE drug_targets ENABLE ROW LEVEL SECURITY;
ALTER TABLE molecule_drug_targets ENABLE ROW LEVEL SECURITY;
ALTER TABLE drug_interactions ENABLE ROW LEVEL SECURITY;

-- Create policies
CREATE POLICY "Public drug_categories are viewable by everyone"
    ON drug_categories FOR SELECT
    USING (true);

CREATE POLICY "Users can view their molecule category mappings"
    ON molecule_drug_categories FOR SELECT
    USING (
        molecule_id IN (
            SELECT id FROM molecule WHERE created_by = auth.uid()
        )
    );

-- Similar policies for other tables...
```

### 3. Identifier Mapping (Week 2)

```python
# drugbank/mapper.py
import logging
from typing import Dict, Optional, List, Tuple
from dataclasses import dataclass

@dataclass
class MappingResult:
    """Result of mapping a DrugBank entity to a CryoProtect molecule"""
    molecule_id: Optional[str]
    confidence: float  # 0.0 to 1.0
    mapping_method: str
    external_id: str
    external_id_type: str = 'DrugBank'

class DrugBankMapper:
    """Maps DrugBank entities to CryoProtect molecules"""
    
    def __init__(self, db_connection):
        """Initialize with database connection
        
        Args:
            db_connection: Database connection object
        """
        self.db = db_connection
        self.logger = logging.getLogger(__name__)
    
    def map_by_structure(self, structure_info: Dict) -> MappingResult:
        """Map DrugBank drug to molecule based on structure information
        
        Args:
            structure_info: Structure information dictionary
            
        Returns:
            MappingResult with mapping details
        """
        # Try InChIKey first (most reliable)
        if structure_info.get("inchikey"):
            molecule_id = self._find_by_inchikey(structure_info["inchikey"])
            if molecule_id:
                return MappingResult(
                    molecule_id=molecule_id,
                    confidence=1.0,
                    mapping_method="inchikey",
                    external_id=structure_info["drugbank_id"]
                )
        
        # Try InChI next
        if structure_info.get("inchi"):
            molecule_id = self._find_by_inchi(structure_info["inchi"])
            if molecule_id:
                return MappingResult(
                    molecule_id=molecule_id,
                    confidence=0.95,
                    mapping_method="inchi",
                    external_id=structure_info["drugbank_id"]
                )
        
        # Try SMILES next
        if structure_info.get("smiles"):
            molecule_id = self._find_by_smiles(structure_info["smiles"])
            if molecule_id:
                return MappingResult(
                    molecule_id=molecule_id,
                    confidence=0.9,
                    mapping_method="smiles",
                    external_id=structure_info["drugbank_id"]
                )
        
        # Try name (least reliable)
        if structure_info.get("name"):
            molecule_id = self._find_by_name(structure_info["name"])
            if molecule_id:
                return MappingResult(
                    molecule_id=molecule_id,
                    confidence=0.7,
                    mapping_method="name",
                    external_id=structure_info["drugbank_id"]
                )
        
        # No mapping found
        return MappingResult(
            molecule_id=None,
            confidence=0.0,
            mapping_method="none",
            external_id=structure_info["drugbank_id"]
        )
    
    def _find_by_inchikey(self, inchikey: str) -> Optional[str]:
        """Find molecule by InChIKey
        
        Args:
            inchikey: InChIKey to search for
            
        Returns:
            Molecule ID if found, None otherwise
        """
        try:
            result = self.db.execute_query(
                "SELECT id FROM molecule WHERE inchikey = %s",
                (inchikey,)
            )
            if result and len(result) > 0:
                return result[0]["id"]
        except Exception as e:
            self.logger.error(f"Error finding molecule by InChIKey: {e}")
        return None
    
    def _find_by_inchi(self, inchi: str) -> Optional[str]:
        """Find molecule by InChI
        
        Args:
            inchi: InChI to search for
            
        Returns:
            Molecule ID if found, None otherwise
        """
        try:
            result = self.db.execute_query(
                "SELECT id FROM molecule WHERE inchi = %s",
                (inchi,)
            )
            if result and len(result) > 0:
                return result[0]["id"]
        except Exception as e:
            self.logger.error(f"Error finding molecule by InChI: {e}")
        return None
    
    def _find_by_smiles(self, smiles: str) -> Optional[str]:
        """Find molecule by SMILES
        
        Args:
            smiles: SMILES to search for
            
        Returns:
            Molecule ID if found, None otherwise
        """
        try:
            result = self.db.execute_query(
                "SELECT id FROM molecule WHERE smiles = %s",
                (smiles,)
            )
            if result and len(result) > 0:
                return result[0]["id"]
        except Exception as e:
            self.logger.error(f"Error finding molecule by SMILES: {e}")
        return None
    
    def _find_by_name(self, name: str) -> Optional[str]:
        """Find molecule by name (fuzzy match)
        
        Args:
            name: Name to search for
            
        Returns:
            Molecule ID if found, None otherwise
        """
        try:
            # Use trigram similarity for fuzzy matching
            result = self.db.execute_query(
                """
                SELECT id, similarity(name, %s) as sim
                FROM molecule
                WHERE similarity(name, %s) > 0.6
                ORDER BY sim DESC
                LIMIT 1
                """,
                (name, name)
            )
            if result and len(result) > 0:
                return result[0]["id"]
        except Exception as e:
            self.logger.error(f"Error finding molecule by name: {e}")
        return None
    
    def create_mapping(self, mapping_result: MappingResult, user_id: str) -> bool:
        """Create mapping record in database
        
        Args:
            mapping_result: Mapping result
            user_id: User ID creating the mapping
            
        Returns:
            True if successful, False otherwise
        """
        if not mapping_result.molecule_id:
            return False
        
        try:
            result = self.db.execute_query(
                """
                INSERT INTO molecule_identifier_mapping
                (molecule_id, external_id, external_id_type, confidence_score, 
                mapping_method, created_by)
                VALUES (%s, %s, %s, %s, %s, %s)
                RETURNING id
                """,
                (
                    mapping_result.molecule_id,
                    mapping_result.external_id,
                    mapping_result.external_id_type,
                    mapping_result.confidence,
                    mapping_result.mapping_method,
                    user_id
                )
            )
            return result and len(result) > 0
        except Exception as e:
            self.logger.error(f"Error creating mapping: {e}")
            return False
```

### 4. Importer Implementation (Week 2-3)

```python
# drugbank/importer.py
import os
import json
import time
import uuid
import logging
from datetime import datetime
from typing import Dict, List, Optional, Any, Tuple

from .client import DrugBankClient, DrugBankCredentials
from .mapper import DrugBankMapper, MappingResult
from database.connection_pool import DatabaseConnectionPool

class DrugBankImporter:
    """Imports data from DrugBank to CryoProtect database"""
    
    def __init__(self, credentials: DrugBankCredentials, 
                checkpoint_dir: str = "./checkpoints",
                batch_size: int = 50):
        """Initialize the DrugBank importer
        
        Args:
            credentials: DrugBank API credentials
            checkpoint_dir: Directory for checkpoints
            batch_size: Batch size for processing
        """
        self.client = DrugBankClient(credentials)
        self.checkpoint_dir = checkpoint_dir
        self.batch_size = batch_size
        self.logger = logging.getLogger(__name__)
        
        # Create checkpoint directory
        os.makedirs(checkpoint_dir, exist_ok=True)
        
        # Get database connection
        self.db = DatabaseConnectionPool.get_instance()
        
        # Create mapper
        self.mapper = DrugBankMapper(self.db)
        
        # Initialize statistics
        self.stats = {
            "total_processed": 0,
            "total_imported": 0,
            "total_mapped": 0,
            "total_new": 0,
            "total_categories": 0,
            "total_targets": 0,
            "total_interactions": 0,
            "total_errors": 0
        }
    
    def _load_checkpoint(self, job_id: str) -> Dict:
        """Load checkpoint for import job
        
        Args:
            job_id: Import job ID
            
        Returns:
            Checkpoint data or empty dict if not found
        """
        checkpoint_path = os.path.join(self.checkpoint_dir, f"drugbank_import_{job_id}.json")
        if os.path.exists(checkpoint_path):
            with open(checkpoint_path, "r") as f:
                return json.load(f)
        return {}
    
    def _save_checkpoint(self, job_id: str, data: Dict) -> None:
        """Save checkpoint for import job
        
        Args:
            job_id: Import job ID
            data: Checkpoint data
        """
        checkpoint_path = os.path.join(self.checkpoint_dir, f"drugbank_import_{job_id}.json")
        data["timestamp"] = datetime.now().isoformat()
        with open(checkpoint_path, "w") as f:
            json.dump(data, f, indent=2)
    
    def _create_import_job(self, description: str, user_id: str) -> str:
        """Create import job record
        
        Args:
            description: Job description
            user_id: User ID
            
        Returns:
            Job ID
        """
        try:
            result = self.db.execute_query(
                """
                INSERT INTO import_job
                (description, status, created_by)
                VALUES (%s, %s, %s)
                RETURNING id
                """,
                (description, "pending", user_id)
            )
            if result and len(result) > 0:
                return result[0]["id"]
        except Exception as e:
            self.logger.error(f"Error creating import job: {e}")
            raise
    
    def _update_import_job(self, job_id: str, status: str, stats: Dict = None) -> None:
        """Update import job status
        
        Args:
            job_id: Job ID
            status: New status
            stats: Optional statistics
        """
        try:
            params = [status]
            sql = "UPDATE import_job SET status = %s"
            
            if status in ["completed", "failed", "cancelled"]:
                sql += ", completed_at = %s"
                params.append(datetime.now().isoformat())
            
            if stats:
                sql += ", metadata = %s"
                params.append(json.dumps(stats))
            
            sql += " WHERE id = %s"
            params.append(job_id)
            
            self.db.execute_query(sql, params)
        except Exception as e:
            self.logger.error(f"Error updating import job: {e}")
    
    def import_drug_list(self, drug_ids: List[str], user_id: str) -> Dict:
        """Import a list of drugs from DrugBank
        
        Args:
            drug_ids: List of DrugBank IDs
            user_id: User ID
            
        Returns:
            Dictionary with import statistics
        """
        # Create import job
        job_id = self._create_import_job(
            f"DrugBank import - {len(drug_ids)} drugs",
            user_id
        )
        
        # Update job status
        self._update_import_job(job_id, "in_progress")
        
        # Load checkpoint if exists
        checkpoint = self._load_checkpoint(job_id)
        processed_ids = set(checkpoint.get("processed_ids", []))
        
        # Process drugs in batches
        try:
            total_batches = (len(drug_ids) + self.batch_size - 1) // self.batch_size
            
            for batch_idx in range(total_batches):
                batch_start = batch_idx * self.batch_size
                batch_end = min(batch_start + self.batch_size, len(drug_ids))
                batch = drug_ids[batch_start:batch_end]
                
                # Skip already processed drugs
                batch = [drug_id for drug_id in batch if drug_id not in processed_ids]
                
                if not batch:
                    self.logger.info(f"Skipping batch {batch_idx+1} - all drugs already processed")
                    continue
                
                self.logger.info(f"Processing batch {batch_idx+1}/{total_batches} ({len(batch)} drugs)")
                
                # Process batch
                self._process_drug_batch(batch, user_id, processed_ids)
                
                # Save checkpoint
                checkpoint["processed_ids"] = list(processed_ids)
                checkpoint["stats"] = self.stats
                self._save_checkpoint(job_id, checkpoint)
                
                # Update job status
                self._update_import_job(
                    job_id, 
                    "in_progress", 
                    {
                        "batch": batch_idx + 1,
                        "total_batches": total_batches,
                        "processed": self.stats["total_processed"],
                        "imported": self.stats["total_imported"]
                    }
                )
            
            # Complete job
            self._update_import_job(job_id, "completed", self.stats)
            
            return self.stats
            
        except Exception as e:
            self.logger.error(f"Error importing drugs: {e}")
            self._update_import_job(job_id, "failed", {"error": str(e)})
            raise
    
    def _process_drug_batch(self, drug_ids: List[str], user_id: str, processed_ids: set) -> None:
        """Process a batch of drugs
        
        Args:
            drug_ids: List of DrugBank IDs
            user_id: User ID
            processed_ids: Set of already processed IDs
        """
        for drug_id in drug_ids:
            try:
                # Fetch drug data
                drug_data = self.client.get_drug_by_id(drug_id)
                if not drug_data:
                    self.logger.warning(f"No data found for drug {drug_id}")
                    self.stats["total_errors"] += 1
                    processed_ids.add(drug_id)
                    continue
                
                # Extract structure info
                structure_info = self.client.extract_structure_info(drug_data)
                if not structure_info.get("smiles") or not structure_info.get("inchi"):
                    self.logger.warning(f"Incomplete structure data for drug {drug_id}")
                    self.stats["total_errors"] += 1
                    processed_ids.add(drug_id)
                    continue
                
                # Map to existing molecule or create new one
                mapping_result = self.mapper.map_by_structure(structure_info)
                
                if mapping_result.molecule_id:
                    # Existing molecule - update with DrugBank info
                    self._update_existing_molecule(mapping_result.molecule_id, drug_data, user_id)
                    self.mapper.create_mapping(mapping_result, user_id)
                    self.stats["total_mapped"] += 1
                else:
                    # New molecule - create from DrugBank data
                    molecule_id = self._create_new_molecule(structure_info, drug_data, user_id)
                    if molecule_id:
                        # Create mapping for the new molecule
                        mapping_result.molecule_id = molecule_id
                        self.mapper.create_mapping(mapping_result, user_id)
                        self.stats["total_new"] += 1
                    else:
                        self.logger.error(f"Failed to create molecule for drug {drug_id}")
                        self.stats["total_errors"] += 1
                        continue
                
                # Import additional drug data
                if mapping_result.molecule_id:
                    self._import_drug_properties(mapping_result.molecule_id, drug_data, user_id)
                    self._import_drug_categories(mapping_result.molecule_id, drug_data, user_id)
                    self._import_drug_targets(mapping_result.molecule_id, drug_data, user_id)
                    self._import_drug_interactions(mapping_result.molecule_id, drug_id, user_id)
                    self.stats["total_imported"] += 1
                
                self.stats["total_processed"] += 1
                processed_ids.add(drug_id)
                
            except Exception as e:
                self.logger.error(f"Error processing drug {drug_id}: {e}")
                self.stats["total_errors"] += 1
                processed_ids.add(drug_id)
                continue
            
            # Rate limiting
            time.sleep(1.0 / self.client.rate_limit)
    
    def _update_existing_molecule(self, molecule_id: str, drug_data: Dict, user_id: str) -> None:
        """Update existing molecule with DrugBank information
        
        Args:
            molecule_id: Molecule ID
            drug_data: DrugBank data
            user_id: User ID
        """
        try:
            # Update molecule record
            self.db.execute_query(
                """
                UPDATE molecule SET
                    drugbank_id = %s,
                    drug_type = %s,
                    approved = %s,
                    updated_at = %s,
                    modification_history = COALESCE(modification_history, '[]'::jsonb) || %s::jsonb
                WHERE id = %s
                """,
                (
                    drug_data.get("drugbank_id"),
                    drug_data.get("type"),
                    drug_data.get("approved") == "true",
                    datetime.now().isoformat(),
                    json.dumps([{
                        "timestamp": datetime.now().isoformat(),
                        "action": "updated_from_drugbank",
                        "user_id": user_id
                    }]),
                    molecule_id
                )
            )
        except Exception as e:
            self.logger.error(f"Error updating molecule {molecule_id}: {e}")
            raise
    
    def _create_new_molecule(self, structure_info: Dict, drug_data: Dict, user_id: str) -> Optional[str]:
        """Create new molecule from DrugBank data
        
        Args:
            structure_info: Structure information
            drug_data: DrugBank data
            user_id: User ID
            
        Returns:
            Molecule ID if created successfully
        """
        try:
            # Create molecule record
            molecule_id = str(uuid.uuid4())
            result = self.db.execute_query(
                """
                INSERT INTO molecule
                (id, name, smiles, inchi, inchikey, formula, molecular_weight,
                drugbank_id, drug_type, approved, created_by, data_source,
                version, modification_history)
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                RETURNING id
                """,
                (
                    molecule_id,
                    structure_info.get("name"),
                    structure_info.get("smiles"),
                    structure_info.get("inchi"),
                    structure_info.get("inchikey"),
                    structure_info.get("formula"),
                    structure_info.get("molecular_weight"),
                    drug_data.get("drugbank_id"),
                    drug_data.get("type"),
                    drug_data.get("approved") == "true",
                    user_id,
                    "DrugBank",
                    1,
                    json.dumps([{
                        "timestamp": datetime.now().isoformat(),
                        "action": "created_from_drugbank",
                        "user_id": user_id
                    }])
                )
            )
            
            if result and len(result) > 0:
                return result[0]["id"]
            return None
            
        except Exception as e:
            self.logger.error(f"Error creating molecule: {e}")
            raise
    
    def _import_drug_properties(self, molecule_id: str, drug_data: Dict, user_id: str) -> None:
        """Import drug properties
        
        Args:
            molecule_id: Molecule ID
            drug_data: DrugBank data
            user_id: User ID
        """
        try:
            # Extract properties
            properties = self.client.extract_properties(drug_data)
            
            # Get property types
            property_types = {}
            result = self.db.execute_query("SELECT id, name, data_type FROM property_types")
            for row in result:
                property_types[row["name"].lower()] = {
                    "id": row["id"],
                    "data_type": row["data_type"]
                }
            
            # Insert properties
            for prop in properties:
                prop_type = property_types.get(prop["name"].lower())
                if not prop_type:
                    self.logger.warning(f"Property type {prop['name']} not found")
                    continue
                
                prop_id = str(uuid.uuid4())
                
                # Determine value column based on data type
                value_col = None
                value_param = None
                
                if prop["type"] == "numeric" or prop_type["data_type"] == "numeric":
                    value_col = "numeric_value"
                    try:
                        value_param = float(prop["value"])
                    except (ValueError, TypeError):
                        continue
                elif prop["type"] == "boolean" or prop_type["data_type"] == "boolean":
                    value_col = "boolean_value"
                    if isinstance(prop["value"], bool):
                        value_param = prop["value"]
                    elif isinstance(prop["value"], str):
                        value_param = prop["value"].lower() == "true"
                    else:
                        continue
                else:
                    value_col = "text_value"
                    value_param = str(prop["value"])
                
                # Insert property
                sql = f"""
                INSERT INTO molecular_property
                (id, molecule_id, property_type_id, {value_col}, created_by, data_source)
                VALUES (%s, %s, %s, %s, %s, %s)
                """
                
                self.db.execute_query(
                    sql,
                    (
                        prop_id,
                        molecule_id,
                        prop_type["id"],
                        value_param,
                        user_id,
                        "DrugBank"
                    )
                )
                
        except Exception as e:
            self.logger.error(f"Error importing properties for molecule {molecule_id}: {e}")
            raise
    
    def _import_drug_categories(self, molecule_id: str, drug_data: Dict, user_id: str) -> None:
        """Import drug categories
        
        Args:
            molecule_id: Molecule ID
            drug_data: DrugBank data
            user_id: User ID
        """
        if "categories" not in drug_data or not isinstance(drug_data["categories"], list):
            return
        
        try:
            for category in drug_data["categories"]:
                # Get or create category
                cat_name = category.get("name")
                if not cat_name:
                    continue
                
                # Check if category exists
                result = self.db.execute_query(
                    "SELECT id FROM drug_categories WHERE name = %s",
                    (cat_name,)
                )
                
                if result and len(result) > 0:
                    cat_id = result[0]["id"]
                else:
                    # Create new category
                    cat_id = str(uuid.uuid4())
                    self.db.execute_query(
                        """
                        INSERT INTO drug_categories
                        (id, name, description)
                        VALUES (%s, %s, %s)
                        """,
                        (
                            cat_id,
                            cat_name,
                            category.get("description")
                        )
                    )
                    self.stats["total_categories"] += 1
                
                # Add molecule-category mapping
                self.db.execute_query(
                    """
                    INSERT INTO molecule_drug_categories
                    (molecule_id, category_id)
                    VALUES (%s, %s)
                    ON CONFLICT (molecule_id, category_id) DO NOTHING
                    """,
                    (molecule_id, cat_id)
                )
                
        except Exception as e:
            self.logger.error(f"Error importing categories for molecule {molecule_id}: {e}")
            raise
    
    def _import_drug_targets(self, molecule_id: str, drug_data: Dict, user_id: str) -> None:
        """Import drug targets
        
        Args:
            molecule_id: Molecule ID
            drug_data: DrugBank data
            user_id: User ID
        """
        if "targets" not in drug_data or not isinstance(drug_data["targets"], list):
            return
        
        try:
            for target in drug_data["targets"]:
                # Get or create target
                target_name = target.get("name")
                if not target_name:
                    continue
                
                # Check if target exists
                result = self.db.execute_query(
                    "SELECT id FROM drug_targets WHERE name = %s",
                    (target_name,)
                )
                
                if result and len(result) > 0:
                    target_id = result[0]["id"]
                else:
                    # Create new target
                    target_id = str(uuid.uuid4())
                    self.db.execute_query(
                        """
                        INSERT INTO drug_targets
                        (id, name, uniprot_id, gene_name, organism, target_type, description)
                        VALUES (%s, %s, %s, %s, %s, %s, %s)
                        """,
                        (
                            target_id,
                            target_name,
                            target.get("uniprot_id"),
                            target.get("gene_name"),
                            target.get("organism"),
                            target.get("target_type"),
                            target.get("description")
                        )
                    )
                    self.stats["total_targets"] += 1
                
                # Add molecule-target mapping
                self.db.execute_query(
                    """
                    INSERT INTO molecule_drug_targets
                    (molecule_id, target_id, action, affinity_value, affinity_unit)
                    VALUES (%s, %s, %s, %s, %s)
                    ON CONFLICT (molecule_id, target_id) DO NOTHING
                    """,
                    (
                        molecule_id,
                        target_id,
                        target.get("action"),
                        target.get("affinity_value"),
                        target.get("affinity_unit")
                    )
                )
                
        except Exception as e:
            self.logger.error(f"Error importing targets for molecule {molecule_id}: {e}")
            raise
    
    def _import_drug_interactions(self, molecule_id: str, drug_id: str, user_id: str) -> None:
        """Import drug interactions
        
        Args:
            molecule_id: Molecule ID
            drug_id: DrugBank ID
            user_id: User ID
        """
        try:
            # Get interactions from API
            interactions = self.client.get_drug_interactions(drug_id)
            if not interactions:
                return
            
            for interaction in interactions:
                other_drug_id = interaction.get("drugbank_id")
                if not other_drug_id:
                    continue
                
                # Find molecule for the other drug
                result = self.db.execute_query(
                    "SELECT id FROM molecule WHERE drugbank_id = %s",
                    (other_drug_id,)
                )
                
                if not result or len(result) == 0:
                    continue
                
                other_molecule_id = result[0]["id"]
                
                # Add interaction
                self.db.execute_query(
                    """
                    INSERT INTO drug_interactions
                    (molecule_id_1, molecule_id_2, description, severity, interaction_type)
                    VALUES (%s, %s, %s, %s, %s)
                    ON CONFLICT (molecule_id_1, molecule_id_2) DO NOTHING
                    """,
                    (
                        molecule_id,
                        other_molecule_id,
                        interaction.get("description"),
                        interaction.get("severity"),
                        interaction.get("type")
                    )
                )
                self.stats["total_interactions"] += 1
                
        except Exception as e:
            self.logger.error(f"Error importing interactions for molecule {molecule_id}: {e}")
            raise
```

### 5. Integration with ML Models (Week 3)

```python
# drugbank/features.py
import numpy as np
from typing import Dict, List, Any, Optional
from rdkit import Chem
from rdkit.Chem import Descriptors

class DrugBankFeatureExtractor:
    """Extract features from DrugBank data for ML models"""
    
    def __init__(self, db_connection):
        """Initialize with database connection
        
        Args:
            db_connection: Database connection
        """
        self.db = db_connection
    
    def extract_drug_class_features(self, molecule_id: str) -> Dict[str, float]:
        """Extract features based on drug classification
        
        Args:
            molecule_id: Molecule ID
            
        Returns:
            Dictionary of features
        """
        features = {}
        
        # Get drug categories
        result = self.db.execute_query(
            """
            SELECT dc.name
            FROM molecule_drug_categories mdc
            JOIN drug_categories dc ON mdc.category_id = dc.id
            WHERE mdc.molecule_id = %s
            """,
            (molecule_id,)
        )
        
        if not result:
            return features
        
        # Define category groups
        category_groups = {
            "antifreeze": ["antifreeze", "cryoprotect"],
            "membrane_stabilizer": ["membrane stabiliz", "cell membrane", "lipid"],
            "protein_stabilizer": ["protein stabiliz", "enzyme stabiliz", "chaperone"],
            "sugar_derivative": ["sugar", "saccharide", "polysaccharide", "disaccharide"],
            "alcohol": ["alcohol", "polyol", "glycol"],
            "dmso_like": ["dmso", "dimethyl sulfoxide", "sulfoxide"],
            "antioxidant": ["antioxidant", "radical scaveng"]
        }
        
        # Initialize all category features to 0
        for group in category_groups:
            features[f"drug_class_{group}"] = 0.0
        
        # Set features based on categories
        categories = [row["name"].lower() for row in result]
        
        for group, keywords in category_groups.items():
            if any(any(keyword in category for keyword in keywords) for category in categories):
                features[f"drug_class_{group}"] = 1.0
        
        return features
    
    def extract_target_features(self, molecule_id: str) -> Dict[str, float]:
        """Extract features based on drug targets
        
        Args:
            molecule_id: Molecule ID
            
        Returns:
            Dictionary of features
        """
        features = {}
        
        # Get drug targets
        result = self.db.execute_query(
            """
            SELECT dt.name, dt.target_type, mdt.action
            FROM molecule_drug_targets mdt
            JOIN drug_targets dt ON mdt.target_id = dt.id
            WHERE mdt.molecule_id = %s
            """,
            (molecule_id,)
        )
        
        if not result:
            return features
        
        # Define target groups
        target_groups = {
            "membrane_protein": ["membrane", "channel", "transporter"],
            "enzyme": ["enzyme", "hydrolase", "oxidase", "reductase"],
            "receptor": ["receptor"],
            "carrier": ["carrier", "binding protein"],
            "structural_protein": ["structural", "cytoskeletal"]
        }
        
        # Initialize all target features to 0
        for group in target_groups:
            features[f"target_{group}"] = 0.0
        
        # Count targets by type
        target_types = {}
        for row in result:
            target_type = row["target_type"].lower() if row["target_type"] else ""
            target_types[target_type] = target_types.get(target_type, 0) + 1
        
        # Set features based on target types
        for target_type in target_types:
            for group, keywords in target_groups.items():
                if any(keyword in target_type for keyword in keywords):
                    features[f"target_{group}"] = target_types[target_type] / len(result)
        
        return features
    
    def extract_pharmacokinetic_features(self, molecule_id: str) -> Dict[str, float]:
        """Extract pharmacokinetic features
        
        Args:
            molecule_id: Molecule ID
            
        Returns:
            Dictionary of features
        """
        features = {}
        
        # Get pharmacokinetic properties
        properties = [
            "Blood-Brain Barrier",
            "Protein Binding",
            "Bioavailability",
            "Half Life",
            "Volume of Distribution",
            "Clearance"
        ]
        
        placeholders = ", ".join(["%s"] * len(properties))
        result = self.db.execute_query(
            f"""
            SELECT pt.name, mp.numeric_value, mp.boolean_value
            FROM molecular_property mp
            JOIN property_types pt ON mp.property_type_id = pt.id
            WHERE mp.molecule_id = %s AND pt.name IN ({placeholders})
            """,
            (molecule_id, *properties)
        )
        
        if not result:
            return features
        
        # Set features based on properties
        for row in result:
            name = row["name"].lower().replace(" ", "_").replace("-", "_")
            if row["numeric_value"] is not None:
                features[f"pk_{name}"] = row["numeric_value"]
            elif row["boolean_value"] is not None:
                features[f"pk_{name}"] = 1.0 if row["boolean_value"] else 0.0
        
        return features
    
    def extract_all_features(self, molecule_id: str) -> Dict[str, float]:
        """Extract all DrugBank features for a molecule
        
        Args:
            molecule_id: Molecule ID
            
        Returns:
            Dictionary of all features
        """
        # Get molecule SMILES
        result = self.db.execute_query(
            "SELECT smiles FROM molecule WHERE id = %s",
            (molecule_id,)
        )
        
        if not result or not result[0]["smiles"]:
            return {}
        
        smiles = result[0]["smiles"]
        
        # Initialize features with basic RDKit descriptors
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return {}
        
        features = {
            "mw": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "tpsa": Descriptors.TPSA(mol),
            "hba": Descriptors.NumHAcceptors(mol),
            "hbd": Descriptors.NumHDonors(mol),
            "rtb": Descriptors.NumRotatableBonds(mol),
            "rings": Descriptors.RingCount(mol)
        }
        
        # Add DrugBank specific features
        drug_class_features = self.extract_drug_class_features(molecule_id)
        target_features = self.extract_target_features(molecule_id)
        pk_features = self.extract_pharmacokinetic_features(molecule_id)
        
        features.update(drug_class_features)
        features.update(target_features)
        features.update(pk_features)
        
        # Get if the molecule is a known drug
        result = self.db.execute_query(
            "SELECT approved FROM molecule WHERE id = %s",
            (molecule_id,)
        )
        
        if result:
            features["is_approved_drug"] = 1.0 if result[0]["approved"] else 0.0
        
        return features
    
    def extract_batch_features(self, molecule_ids: List[str]) -> Dict[str, Dict[str, float]]:
        """Extract features for multiple molecules
        
        Args:
            molecule_ids: List of molecule IDs
            
        Returns:
            Dictionary mapping molecule IDs to feature dictionaries
        """
        features = {}
        
        for molecule_id in molecule_ids:
            features[molecule_id] = self.extract_all_features(molecule_id)
        
        return features
```

### 6. API Integration (Week 4)

```python
# api/drugbank_resources.py
from flask import request, jsonify, Blueprint
from flask_restx import Api, Resource, fields, Namespace

from drugbank.client import DrugBankClient, DrugBankCredentials
from drugbank.importer import DrugBankImporter
from drugbank.features import DrugBankFeatureExtractor
from api.utils import get_user_id, login_required, admin_required
from config import Config

# Create blueprint
drugbank_bp = Blueprint('drugbank', __name__, url_prefix='/api/v1/drugbank')
drugbank_ns = Namespace('drugbank', description='DrugBank operations')

# Initialize client
config = Config()
credentials = DrugBankCredentials(
    api_key=config.DRUGBANK_API_KEY,
    app_id=config.DRUGBANK_APP_ID
)
client = DrugBankClient(credentials)

# Models
drug_model = drugbank_ns.model('Drug', {
    'drugbank_id': fields.String(required=True, description='DrugBank ID'),
    'name': fields.String(required=True, description='Drug name'),
    'smiles': fields.String(description='SMILES'),
    'inchi': fields.String(description='InChI'),
    'inchikey': fields.String(description='InChIKey'),
    'formula': fields.String(description='Molecular formula'),
    'molecular_weight': fields.Float(description='Molecular weight'),
    'approved': fields.Boolean(description='Approval status')
})

import_model = drugbank_ns.model('ImportRequest', {
    'drug_ids': fields.List(fields.String, required=True, description='List of DrugBank IDs to import'),
    'description': fields.String(description='Import job description')
})

import_result_model = drugbank_ns.model('ImportResult', {
    'job_id': fields.String(description='Import job ID'),
    'status': fields.String(description='Job status'),
    'total_processed': fields.Integer(description='Total drugs processed'),
    'total_imported': fields.Integer(description='Total drugs imported'),
    'total_mapped': fields.Integer(description='Total drugs mapped to existing molecules'),
    'total_new': fields.Integer(description='Total new molecules created'),
    'total_errors': fields.Integer(description='Total errors encountered')
})

# Routes
@drugbank_ns.route('/drug/<string:drugbank_id>')
class DrugResource(Resource):
    @drugbank_ns.marshal_with(drug_model)
    @login_required
    def get(self, drugbank_id):
        """Get drug information from DrugBank"""
        try:
            drug_data = client.get_drug_by_id(drugbank_id)
            return drug_data
        except Exception as e:
            drugbank_ns.abort(500, f"Error fetching drug: {str(e)}")

@drugbank_ns.route('/search')
class DrugSearchResource(Resource):
    @drugbank_ns.doc(params={'q': 'Search query', 'page': 'Page number', 'per_page': 'Results per page'})
    @login_required
    def get(self):
        """Search for drugs in DrugBank"""
        query = request.args.get('q', '')
        page = int(request.args.get('page', 1))
        per_page = int(request.args.get('per_page', 10))
        
        try:
            results = client.search_drugs(query, page, per_page)
            return results
        except Exception as e:
            drugbank_ns.abort(500, f"Error searching drugs: {str(e)}")

@drugbank_ns.route('/import')
class DrugImportResource(Resource):
    @drugbank_ns.expect(import_model)
    @drugbank_ns.marshal_with(import_result_model)
    @admin_required
    def post(self):
        """Import drugs from DrugBank"""
        data = request.json
        drug_ids = data.get('drug_ids', [])
        description = data.get('description', f"Import {len(drug_ids)} drugs")
        
        if not drug_ids:
            drugbank_ns.abort(400, "No drug IDs provided")
        
        user_id = get_user_id()
        
        try:
            importer = DrugBankImporter(credentials)
            result = importer.import_drug_list(drug_ids, user_id)
            
            return {
                'job_id': result.get('job_id', ''),
                'status': 'completed',
                'total_processed': result.get('total_processed', 0),
                'total_imported': result.get('total_imported', 0),
                'total_mapped': result.get('total_mapped', 0),
                'total_new': result.get('total_new', 0),
                'total_errors': result.get('total_errors', 0)
            }
        except Exception as e:
            drugbank_ns.abort(500, f"Error importing drugs: {str(e)}")

@drugbank_ns.route('/features/<string:molecule_id>')
class DrugFeaturesResource(Resource):
    @login_required
    def get(self, molecule_id):
        """Get DrugBank features for a molecule"""
        from database.connection_pool import DatabaseConnectionPool
        db = DatabaseConnectionPool.get_instance()
        
        try:
            extractor = DrugBankFeatureExtractor(db)
            features = extractor.extract_all_features(molecule_id)
            return {'molecule_id': molecule_id, 'features': features}
        except Exception as e:
            drugbank_ns.abort(500, f"Error extracting features: {str(e)}")

# Register namespace with blueprint
api = Api(drugbank_bp, doc='/doc')
api.add_namespace(drugbank_ns)
```

### 7. ML Model Extension (Week 4)

```python
# Update predictive_models.py to include DrugBank features

def _extract_features(self, molecule_data: Dict[str, Any]) -> List[float]:
    """
    Extracts scientifically relevant features from molecule data for model input.
    
    Now enhanced with DrugBank features.
    """
    features = []
    
    # Standard features (hydrogen bonding, LogP, etc.)
    # [Original code here...]
    
    # DrugBank features (if available)
    drugbank_features = molecule_data.get('drugbank_features', {})
    
    if drugbank_features:
        # Drug class features
        for drug_class in ['antifreeze', 'membrane_stabilizer', 'protein_stabilizer', 
                          'sugar_derivative', 'alcohol', 'dmso_like', 'antioxidant']:
            feature_name = f'drug_class_{drug_class}'
            features.append(drugbank_features.get(feature_name, 0.0))
        
        # Target features
        for target_type in ['membrane_protein', 'enzyme', 'receptor', 
                          'carrier', 'structural_protein']:
            feature_name = f'target_{target_type}'
            features.append(drugbank_features.get(feature_name, 0.0))
        
        # Pharmacokinetic features
        for pk_prop in ['blood_brain_barrier', 'protein_binding', 'bioavailability', 
                       'half_life', 'volume_of_distribution', 'clearance']:
            feature_name = f'pk_{pk_prop}'
            features.append(drugbank_features.get(feature_name, 0.0))
        
        # Approval status
        features.append(drugbank_features.get('is_approved_drug', 0.0))
    else:
        # Add zeros for all DrugBank features if not available
        features.extend([0.0] * 16)  # 7 drug classes + 5 target types + 3 PK props + 1 approval
    
    return features

def _get_feature_names(self) -> List[str]:
    """Get the names of the features used by the model."""
    standard_features = [
        'h_bond_donors', 'h_bond_acceptors', 'h_bond_total',
        'logp', 'tpsa',
        'molecular_weight', 'heavy_atom_count', 'rotatable_bond_count', 'ring_count', 'fraction_csp3',
        'hydroxyl_count', 'alcohol_count', 'ether_count', 'amine_count', 'amide_count',
        'rule_of_5_violations', 'veber_violations', 'bbb_permeant', 'intestinal_absorption', 'estimated_log_papp'
    ]
    
    # Add DrugBank features
    drugbank_features = []
    
    # Drug class features
    for drug_class in ['antifreeze', 'membrane_stabilizer', 'protein_stabilizer', 
                      'sugar_derivative', 'alcohol', 'dmso_like', 'antioxidant']:
        drugbank_features.append(f'drug_class_{drug_class}')
    
    # Target features
    for target_type in ['membrane_protein', 'enzyme', 'receptor', 
                       'carrier', 'structural_protein']:
        drugbank_features.append(f'target_{target_type}')
    
    # Pharmacokinetic features
    for pk_prop in ['blood_brain_barrier', 'protein_binding', 'bioavailability', 
                   'half_life', 'volume_of_distribution', 'clearance']:
        drugbank_features.append(f'pk_{pk_prop}')
    
    # Approval status
    drugbank_features.append('is_approved_drug')
    
    return standard_features + drugbank_features
```

## Testing Strategy

1. **Unit Tests**: Create comprehensive unit tests for each component:
   - `test_drugbank_client.py`: Test API interactions
   - `test_drugbank_mapper.py`: Test identifier mapping
   - `test_drugbank_importer.py`: Test import functionality
   - `test_drugbank_features.py`: Test feature extraction

2. **Integration Tests**: Create tests for data flow:
   - Test DrugBank import process with mock API responses
   - Test integration with existing molecules
   - Test feature extraction with real data
   - Test API endpoints

3. **Performance Tests**:
   - Test batch processing with large data sets
   - Test connection pooling and cache performance
   - Test feature extraction speed

## Deployment Plan

1. **Database Migration**:
   - Apply schema updates using Supabase migrations
   - Add property types and categories

2. **API Updates**:
   - Deploy new API endpoints for DrugBank integration
   - Update ML endpoints to include DrugBank features

3. **Initial Data Import**:
   - Import curated list of known cryoprotectants from DrugBank
   - Import related compounds with similar structures

4. **Documentation**:
   - Update API documentation with new endpoints
   - Create user guide for DrugBank data integration

## Monitoring and Maintenance

1. **Logging**:
   - Track API usage and rate limits
   - Monitor import job progress and success rates
   - Track feature usage in ML models

2. **Alerts**:
   - Set up alerts for API rate limit approaching
   - Monitor for failed imports

3. **Performance**:
   - Track query performance for DrugBank-enhanced features
   - Monitor memory usage during batch imports

## Success Criteria

The DrugBank integration will be considered successful when:

1. Users can search and import drugs from DrugBank interface
2. ML models show improved prediction accuracy with DrugBank features
3. Imported drug structures are properly mapped to existing molecules
4. Drug categories and targets provide meaningful insights in the UI
5. Batch import performance meets targets (<2 minutes for 100 drugs)

## References

1. DrugBank API Documentation: https://docs.drugbank.com/v1/
2. DrugBank Data Model: https://go.drugbank.com/documentation/data-model
3. Mapping Guidelines: https://go.drugbank.com/documentation/mapping
4. DrugBank Features for ML: https://go.drugbank.com/documentation/features