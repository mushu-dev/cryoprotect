"""
Tox21 Client for CryoProtect v2.

This module provides functionality to fetch, process, and import toxicity data from the Tox21 program.
Tox21 is a federal collaboration among EPA, NIH, and FDA that uses high-throughput screening
approaches to identify environmental chemicals that could potentially affect human health.

References:
- Tox21 Data: https://ntp.niehs.nih.gov/whatwestudy/tox21/index.html
- Tox21 Dashboard: https://tripod.nih.gov/tox21/assays/
"""

import os
import csv
import json
import logging
import requests
import pandas as pd
import numpy as np
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Any, Union
from pathlib import Path
import time
import zipfile
import io

# Database imports
from supabase import create_client, Client

# Local imports
from config import Config

# Set up logging
logger = logging.getLogger(__name__)

# Tox21 API and data constants
TOX21_DATA_URL = "https://tripod.nih.gov/tox21/assays/download"
TOX21_ASSAY_SUMMARY_FILE = "Tox21_Assay_Summary.csv"
TOX21_CHEMICAL_DATA_FILE = "Tox21_Chemical_Data.csv"

# Cache directory
CACHE_DIR = Path("cache/tox21")


class Tox21Client:
    """Client for interacting with Tox21 data."""

    def __init__(self, supabase_client: Optional[Client] = None, config: Optional[Config] = None):
        """Initialize the Tox21 client."""
        self.config = config or Config()
        
        # Initialize Supabase client if not provided
        if supabase_client:
            self.supabase = supabase_client
        else:
            supabase_url = os.environ.get("SUPABASE_URL") or self.config.SUPABASE_URL
            supabase_key = os.environ.get("SUPABASE_KEY") or self.config.SUPABASE_KEY
            self.supabase = create_client(supabase_url, supabase_key)
        
        # Create cache directory if it doesn't exist
        CACHE_DIR.mkdir(parents=True, exist_ok=True)
        
        # Get Tox21 source ID from database
        self.tox21_source_id = self._get_tox21_source_id()
        
        # Initialize rate limiting parameters
        self.request_count = 0
        self.last_request_time = 0
        self.rate_limit = 5  # requests per second
        
    def _get_tox21_source_id(self) -> str:
        """Get the Tox21 source ID from the database."""
        try:
            response = self.supabase.table("toxicity_data_source").select("id").eq("name", "Tox21").execute()
            if response.data and len(response.data) > 0:
                return response.data[0]["id"]
            else:
                logger.error("Tox21 data source not found in database")
                raise ValueError("Tox21 data source not found in database. Run migrations first.")
        except Exception as e:
            logger.error(f"Error getting Tox21 source ID: {str(e)}")
            raise
    
    def _rate_limit_request(self):
        """Apply rate limiting to API requests."""
        current_time = time.time()
        elapsed = current_time - self.last_request_time
        
        # If we've made a request recently, sleep to respect rate limit
        if elapsed < (1.0 / self.rate_limit):
            sleep_time = (1.0 / self.rate_limit) - elapsed
            time.sleep(sleep_time)
        
        self.last_request_time = time.time()
        self.request_count += 1
    
    def _get_cache_path(self, filename: str) -> Path:
        """Get the path to a cached file."""
        return CACHE_DIR / filename
    
    def _is_cache_valid(self, cache_path: Path, max_age_days: int = 30) -> bool:
        """Check if a cached file is valid (exists and is not too old)."""
        if not cache_path.exists():
            return False
        
        # Check if file is older than max_age_days
        file_age = datetime.now().timestamp() - cache_path.stat().st_mtime
        max_age_seconds = max_age_days * 24 * 60 * 60
        
        return file_age < max_age_seconds
    
    def download_assay_data(self, force_refresh: bool = False) -> Path:
        """Download Tox21 assay summary data."""
        cache_path = self._get_cache_path(TOX21_ASSAY_SUMMARY_FILE)
        
        # Use cached file if it exists and is valid, unless force_refresh is True
        if not force_refresh and self._is_cache_valid(cache_path):
            logger.info(f"Using cached Tox21 assay data: {cache_path}")
            return cache_path
        
        # Download the file
        logger.info("Downloading Tox21 assay summary data...")
        url = f"{TOX21_DATA_URL}/{TOX21_ASSAY_SUMMARY_FILE}"
        
        try:
            self._rate_limit_request()
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            with open(cache_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    f.write(chunk)
            
            logger.info(f"Downloaded Tox21 assay data to {cache_path}")
            return cache_path
        
        except requests.exceptions.RequestException as e:
            logger.error(f"Error downloading Tox21 assay data: {str(e)}")
            raise
    
    def download_chemical_data(self, force_refresh: bool = False) -> Path:
        """Download Tox21 chemical data."""
        cache_path = self._get_cache_path(TOX21_CHEMICAL_DATA_FILE)
        
        # Use cached file if it exists and is valid, unless force_refresh is True
        if not force_refresh and self._is_cache_valid(cache_path):
            logger.info(f"Using cached Tox21 chemical data: {cache_path}")
            return cache_path
        
        # Download the file
        logger.info("Downloading Tox21 chemical data...")
        url = f"{TOX21_DATA_URL}/{TOX21_CHEMICAL_DATA_FILE}"
        
        try:
            self._rate_limit_request()
            response = requests.get(url, stream=True)
            response.raise_for_status()
            
            # Check if it's a zip file
            if url.endswith('.zip') or response.headers.get('Content-Type') == 'application/zip':
                # Extract the CSV file from the zip
                with zipfile.ZipFile(io.BytesIO(response.content)) as z:
                    csv_files = [f for f in z.namelist() if f.endswith('.csv')]
                    if not csv_files:
                        raise ValueError("No CSV files found in the downloaded zip file")
                    
                    # Extract the first CSV file
                    with z.open(csv_files[0]) as csv_file, open(cache_path, 'wb') as f:
                        f.write(csv_file.read())
            else:
                # Direct CSV download
                with open(cache_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
            
            logger.info(f"Downloaded Tox21 chemical data to {cache_path}")
            return cache_path
        
        except requests.exceptions.RequestException as e:
            logger.error(f"Error downloading Tox21 chemical data: {str(e)}")
            raise
    
    def import_assays(self, assay_file: Optional[Path] = None) -> int:
        """Import Tox21 assay data into the database."""
        # Download assay data if file not provided
        if assay_file is None:
            assay_file = self.download_assay_data()
        
        logger.info(f"Importing Tox21 assay data from {assay_file}...")
        
        # Read assay data
        try:
            assay_df = pd.read_csv(assay_file)
            logger.info(f"Read {len(assay_df)} assays from file")
        except Exception as e:
            logger.error(f"Error reading assay data: {str(e)}")
            raise
        
        # Process and import assays
        assays_imported = 0
        batch_size = 100
        
        for i in range(0, len(assay_df), batch_size):
            batch = assay_df.iloc[i:i+batch_size]
            
            # Prepare assay data for import
            assay_data = []
            for _, row in batch.iterrows():
                assay_data.append({
                    "source_id": self.tox21_source_id,
                    "assay_name": row.get("assay_name", ""),
                    "assay_id": row.get("assay_id", ""),
                    "description": row.get("assay_description", ""),
                    "assay_type": row.get("assay_type", ""),
                    "organism": row.get("organism", ""),
                    "tissue": row.get("tissue", ""),
                    "cell_line": row.get("cell_line", ""),
                    "assay_target": row.get("intended_target", ""),
                    "assay_function_type": row.get("assay_function_type", ""),
                    "detection_technology": row.get("detection_technology", ""),
                    "key_assay_reagent": row.get("key_assay_reagent", ""),
                    "assay_footprint": row.get("assay_footprint", ""),
                    "assay_format": row.get("assay_format", ""),
                    "biological_process_target": row.get("biological_process_target", ""),
                    "detection_technology_type": row.get("detection_technology_type", ""),
                    "signal_direction": row.get("signal_direction", "")
                })
            
            # Import assay data
            try:
                response = self.supabase.table("toxicity_assay").upsert(assay_data).execute()
                assays_imported += len(response.data)
                logger.info(f"Imported {len(response.data)} assays (batch {i//batch_size + 1})")
            except Exception as e:
                logger.error(f"Error importing assays: {str(e)}")
                raise
        
        logger.info(f"Successfully imported {assays_imported} Tox21 assays")
        return assays_imported
    
    def import_chemical_data(self, chemical_file: Optional[Path] = None, 
                            batch_size: int = 1000, 
                            max_compounds: Optional[int] = None) -> Dict[str, int]:
        """Import Tox21 chemical data into the database."""
        # Download chemical data if file not provided
        if chemical_file is None:
            chemical_file = self.download_chemical_data()
        
        logger.info(f"Importing Tox21 chemical data from {chemical_file}...")
        
        # Create import job record
        import_job_id = self._create_import_job()
        
        # Read chemical data
        try:
            # Use pandas to read the file efficiently
            chem_df = pd.read_csv(chemical_file)
            
            # Limit to max_compounds if specified
            if max_compounds is not None:
                chem_df = chem_df.head(max_compounds)
                
            logger.info(f"Read {len(chem_df)} chemical records from file")
        except Exception as e:
            logger.error(f"Error reading chemical data: {str(e)}")
            self._update_import_job(import_job_id, "failed", error=str(e))
            raise
        
        # Get assay mapping from database
        assay_mapping = self._get_assay_mapping()
        
        # Process and import chemical data
        stats = {
            "total_records": len(chem_df),
            "compounds_processed": 0,
            "mappings_created": 0,
            "toxicity_data_imported": 0,
            "failed_records": 0
        }
        
        # Group by compound to process one compound at a time
        # Note: Tox21 might use different compound identifiers (e.g., NCGC ID)
        compound_id_column = self._identify_compound_id_column(chem_df)
        compound_groups = chem_df.groupby(compound_id_column)
        
        for i, (compound_id, compound_df) in enumerate(compound_groups):
            if i % 100 == 0:
                logger.info(f"Processing compound {i+1}/{len(compound_groups)} ({compound_id_column}: {compound_id})")
                # Update import job status
                self._update_import_job(import_job_id, "in_progress", 
                                       records_processed=stats["compounds_processed"],
                                       records_imported=stats["toxicity_data_imported"])
            
            try:
                # Map compound to internal molecule ID
                mapping_result = self._map_compound(compound_id, compound_id_column, compound_df)
                
                if mapping_result["molecule_id"]:
                    # Import toxicity data for this compound
                    import_result = self._import_compound_data(
                        mapping_result["molecule_id"],
                        mapping_result["mapping_id"],
                        compound_df,
                        assay_mapping
                    )
                    
                    stats["mappings_created"] += mapping_result["created"]
                    stats["toxicity_data_imported"] += import_result["imported"]
                else:
                    logger.warning(f"Could not map compound {compound_id} to internal molecule")
                
                stats["compounds_processed"] += 1
                
            except Exception as e:
                logger.error(f"Error processing compound {compound_id}: {str(e)}")
                stats["failed_records"] += 1
        
        # Update import job status
        status = "completed" if stats["failed_records"] == 0 else "completed_with_errors"
        self._update_import_job(import_job_id, status, 
                               records_processed=stats["compounds_processed"],
                               records_imported=stats["toxicity_data_imported"],
                               records_failed=stats["failed_records"])
        
        logger.info(f"Tox21 import complete: {stats}")
        return stats
    
    def _identify_compound_id_column(self, df: pd.DataFrame) -> str:
        """Identify the column containing compound identifiers in the Tox21 data."""
        # Common identifier columns in Tox21 data
        possible_columns = ["NCGC_ID", "TOX21_ID", "compound_id", "cid", "casrn"]
        
        for col in possible_columns:
            if col in df.columns:
                return col
        
        # If none of the expected columns are found, use the first column
        logger.warning(f"Could not identify compound ID column, using {df.columns[0]}")
        return df.columns[0]
    
    def _create_import_job(self) -> str:
        """Create a new import job record."""
        try:
            response = self.supabase.table("toxicity_import_job").insert({
                "source_id": self.tox21_source_id,
                "status": "pending",
                "started_at": datetime.now().isoformat()
            }).execute()
            
            if response.data and len(response.data) > 0:
                return response.data[0]["id"]
            else:
                logger.error("Failed to create import job")
                raise ValueError("Failed to create import job")
        except Exception as e:
            logger.error(f"Error creating import job: {str(e)}")
            raise
    
    def _update_import_job(self, job_id: str, status: str, 
                          records_processed: int = None, 
                          records_imported: int = None,
                          records_failed: int = None,
                          error: str = None) -> None:
        """Update an import job record."""
        update_data = {"status": status, "updated_at": datetime.now().isoformat()}
        
        if status in ["completed", "completed_with_errors", "failed"]:
            update_data["completed_at"] = datetime.now().isoformat()
        
        if records_processed is not None:
            update_data["records_processed"] = records_processed
        
        if records_imported is not None:
            update_data["records_imported"] = records_imported
            
        if records_failed is not None:
            update_data["records_failed"] = records_failed
        
        if error:
            # Fetch current error log
            try:
                response = self.supabase.table("toxicity_import_job").select("error_log").eq("id", job_id).execute()
                if response.data and len(response.data) > 0:
                    current_log = response.data[0].get("error_log", [])
                    if not isinstance(current_log, list):
                        current_log = []
                else:
                    current_log = []
            except Exception as e:
                logger.error(f"Error fetching error log: {str(e)}")
                current_log = []
            
            # Add new error to log
            current_log.append({
                "timestamp": datetime.now().isoformat(),
                "error": error
            })
            update_data["error_log"] = current_log
        
        try:
            self.supabase.table("toxicity_import_job").update(update_data).eq("id", job_id).execute()
        except Exception as e:
            logger.error(f"Error updating import job: {str(e)}")
    
    def _get_assay_mapping(self) -> Dict[str, str]:
        """Get mapping of assay_id to UUID from database."""
        try:
            response = self.supabase.table("toxicity_assay").select("id, assay_id").eq("source_id", self.tox21_source_id).execute()
            
            assay_mapping = {}
            if response.data:
                for assay in response.data:
                    assay_mapping[assay["assay_id"]] = assay["id"]
            
            logger.info(f"Loaded {len(assay_mapping)} assay mappings from database")
            return assay_mapping
        except Exception as e:
            logger.error(f"Error getting assay mapping: {str(e)}")
            raise
    
    def _map_compound(self, compound_id: str, id_type: str, compound_df: pd.DataFrame) -> Dict[str, Any]:
        """Map a Tox21 compound to an internal molecule ID."""
        # Extract compound identifiers
        identifiers = {
            id_type.lower(): compound_id
        }
        
        # Get the first row for basic compound info
        first_row = compound_df.iloc[0]
        
        # Add other identifiers if available
        for col, id_type_name in [
            ("casrn", "casrn"),
            ("name", "name"),
            ("preferred_name", "name"),
            ("smiles", "smiles"),
            ("inchi", "inchi"),
            ("inchikey", "inchikey"),
            ("pubchem_cid", "pubchem_cid")
        ]:
            if col in first_row and pd.notna(first_row[col]):
                identifiers[id_type_name] = first_row[col]
        
        # Check if mapping already exists
        try:
            response = self.supabase.table("molecule_identifier_mapping").select("id, molecule_id").eq("source_id", self.tox21_source_id).eq("external_id_type", id_type.upper()).eq("external_id", compound_id).execute()
            
            if response.data and len(response.data) > 0:
                # Mapping exists
                return {
                    "molecule_id": response.data[0]["molecule_id"],
                    "mapping_id": response.data[0]["id"],
                    "created": 0
                }
        except Exception as e:
            logger.error(f"Error checking existing mapping: {str(e)}")
        
        # No mapping exists, try to map using available identifiers
        molecule_id = None
        confidence = 0.0
        mapping_method = None
        
        # Try mapping by InChIKey (most reliable)
        if "inchikey" in identifiers:
            try:
                response = self.supabase.table("molecule").select("id").eq("inchikey", identifiers["inchikey"]).execute()
                if response.data and len(response.data) > 0:
                    molecule_id = response.data[0]["id"]
                    confidence = 1.0
                    mapping_method = "inchikey_match"
            except Exception as e:
                logger.error(f"Error mapping by InChIKey: {str(e)}")
        
        # Try mapping by InChI if InChIKey failed
        if not molecule_id and "inchi" in identifiers:
            try:
                response = self.supabase.table("molecule").select("id").eq("inchi", identifiers["inchi"]).execute()
                if response.data and len(response.data) > 0:
                    molecule_id = response.data[0]["id"]
                    confidence = 0.95
                    mapping_method = "inchi_match"
            except Exception as e:
                logger.error(f"Error mapping by InChI: {str(e)}")
        
        # Try mapping by SMILES if InChI failed
        if not molecule_id and "smiles" in identifiers:
            try:
                response = self.supabase.table("molecule").select("id").eq("smiles", identifiers["smiles"]).execute()
                if response.data and len(response.data) > 0:
                    molecule_id = response.data[0]["id"]
                    confidence = 0.9
                    mapping_method = "smiles_match"
            except Exception as e:
                logger.error(f"Error mapping by SMILES: {str(e)}")
        
        # If mapping successful, create mapping record
        if molecule_id:
            try:
                mapping_data = {
                    "molecule_id": molecule_id,
                    "source_id": self.tox21_source_id,
                    "external_id": compound_id,
                    "external_id_type": id_type.upper(),
                    "confidence_score": confidence,
                    "mapping_method": mapping_method
                }
                
                response = self.supabase.table("molecule_identifier_mapping").insert(mapping_data).execute()
                
                if response.data and len(response.data) > 0:
                    return {
                        "molecule_id": molecule_id,
                        "mapping_id": response.data[0]["id"],
                        "created": 1
                    }
                else:
                    logger.error(f"Failed to create mapping for {compound_id}")
                    return {
                        "molecule_id": molecule_id,
                        "mapping_id": None,
                        "created": 0
                    }
            except Exception as e:
                logger.error(f"Error creating mapping: {str(e)}")
                return {
                    "molecule_id": molecule_id,
                    "mapping_id": None,
                    "created": 0
                }
        
        # No mapping found
        return {
            "molecule_id": None,
            "mapping_id": None,
            "created": 0
        }
    
    def _import_compound_data(self, molecule_id: str, mapping_id: str, 
                             compound_df: pd.DataFrame, 
                             assay_mapping: Dict[str, str]) -> Dict[str, int]:
        """Import toxicity data for a compound."""
        stats = {"imported": 0}
        
        # Prepare toxicity data for import
        toxicity_data = []
        
        for _, row in compound_df.iterrows():
            assay_id = row.get("assay_id")
            if assay_id not in assay_mapping:
                logger.warning(f"Assay {assay_id} not found in mapping, skipping")
                continue
            
            # Extract activity data
            activity_value = row.get("activity_value")
            if pd.isna(activity_value):
                activity_value = None
            
            hit_call = row.get("hit_call")
            if isinstance(hit_call, str):
                hit_call = hit_call.lower() == "true" or hit_call == "1" or hit_call.lower() == "active"
            elif isinstance(hit_call, (int, float)):
                hit_call = hit_call == 1
            else:
                hit_call = None
            
            # Create toxicity data record
            toxicity_data.append({
                "molecule_id": molecule_id,
                "assay_id": assay_mapping[assay_id],
                "mapping_id": mapping_id,
                "activity_value": activity_value,
                "activity_unit": row.get("activity_unit"),
                "activity_type": row.get("activity_type"),
                "hit_call": hit_call,
                "significance": row.get("significance"),
                "reliability_score": 0.8,  # Default reliability score
                "data_source": "Tox21",
                "version": 1
            })
        
        # Import toxicity data in batches
        batch_size = 100
        for i in range(0, len(toxicity_data), batch_size):
            batch = toxicity_data[i:i+batch_size]
            
            try:
                response = self.supabase.table("toxicity_data").upsert(batch).execute()
                stats["imported"] += len(response.data)
            except Exception as e:
                logger.error(f"Error importing toxicity data batch: {str(e)}")
        
        # Update toxicity_data_available flag on molecule
        try:
            self.supabase.table("molecule").update({"toxicity_data_available": True}).eq("id", molecule_id).execute()
        except Exception as e:
            logger.error(f"Error updating toxicity_data_available flag: {str(e)}")
        
        return stats
    
    def calculate_toxicity_scores(self, molecule_ids: Optional[List[str]] = None) -> int:
        """Calculate toxicity scores for molecules."""
        # Get calculation method ID for toxicity scoring
        method_id = self._get_toxicity_calculation_method_id()
        
        # Get molecules to calculate scores for
        if molecule_ids:
            try:
                response = self.supabase.table("molecule").select("id").in_("id", molecule_ids).eq("toxicity_data_available", True).execute()
                molecules = [m["id"] for m in response.data]
            except Exception as e:
                logger.error(f"Error getting molecules: {str(e)}")
                raise
        else:
            try:
                response = self.supabase.table("molecule").select("id").eq("toxicity_data_available", True).execute()
                molecules = [m["id"] for m in response.data]
            except Exception as e:
                logger.error(f"Error getting molecules: {str(e)}")
                raise
        
        logger.info(f"Calculating toxicity scores for {len(molecules)} molecules")
        
        scores_calculated = 0
        
        for molecule_id in molecules:
            try:
                # Get toxicity data for molecule
                response = self.supabase.table("toxicity_data").select("*").eq("molecule_id", molecule_id).eq("data_source", "Tox21").execute()
                
                if not response.data:
                    logger.warning(f"No Tox21 data found for molecule {molecule_id}")
                    continue
                
                # Calculate overall toxicity score
                toxicity_data = response.data
                overall_score = self._calculate_overall_toxicity_score(toxicity_data)
                
                # Save overall score to molecule
                self.supabase.table("molecule").update({"toxicity_score": overall_score}).eq("id", molecule_id).execute()
                
                scores_calculated += 1
                
                if scores_calculated % 100 == 0:
                    logger.info(f"Calculated scores for {scores_calculated}/{len(molecules)} molecules")
                
            except Exception as e:
                logger.error(f"Error calculating toxicity score for molecule {molecule_id}: {str(e)}")
        
        logger.info(f"Successfully calculated toxicity scores for {scores_calculated} molecules")
        return scores_calculated
    
    def _get_toxicity_calculation_method_id(self) -> str:
        """Get the toxicity calculation method ID from the database."""
        try:
            response = self.supabase.table("calculation_method").select("id").eq("name", "Tox21 Score").execute()
            if response.data and len(response.data) > 0:
                return response.data[0]["id"]
            else:
                logger.error("Toxicity calculation method not found in database")
                raise ValueError("Toxicity calculation method not found in database. Run migrations first.")
        except Exception as e:
            logger.error(f"Error getting toxicity calculation method ID: {str(e)}")
            raise
    
    def _calculate_overall_toxicity_score(self, toxicity_data: List[Dict]) -> float:
        """Calculate an overall toxicity score based on toxicity data."""
        if not toxicity_data:
            return 0.0
        
        # Count positive hits
        hit_count = sum(1 for record in toxicity_data if record.get("hit_call") is True)
        
        # Calculate hit ratio
        hit_ratio = hit_count / len(toxicity_data)
        
        # Scale to 0-100
        toxicity_score = hit_ratio * 100
        
        return toxicity_score
