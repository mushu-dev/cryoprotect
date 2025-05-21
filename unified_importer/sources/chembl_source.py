"""
ChEMBL data source implementation.

This module provides the ChEMBLDataSource class for importing compound
data from the ChEMBL database using both the REST API and the official
ChEMBL Web Resources client.
"""

import asyncio
import json
import time
import logging
import aiohttp
import re
import uuid
from typing import Dict, List, Any, Optional, Tuple, Union, Set, AsyncIterator

# Import the official ChEMBL client for more efficient access
try:
    from chembl_webresource_client.new_client import new_client
    CHEMBL_CLIENT_AVAILABLE = True
except ImportError:
    CHEMBL_CLIENT_AVAILABLE = False

from .source_base import MolecularDataSource
from ..core.database import DatabaseOperations
from ..core.checkpoint import CheckpointManager
from ..core.progress import ProgressTracker
from ..core.validation import MoleculeValidator


class ChEMBLDataSource(MolecularDataSource):
    """
    ChEMBL data source for molecular data.
    
    This class provides methods to interact with the ChEMBL database
    via its REST API (https://www.ebi.ac.uk/chembl/api/data).
    """
    
    def __init__(
        self,
        db_operations: DatabaseOperations,
        checkpoint_manager: Optional[CheckpointManager] = None,
        progress_tracker: Optional[ProgressTracker] = None,
        validator: Optional[MoleculeValidator] = None,
        config: Optional[Dict[str, Any]] = None,
        logger: Optional[logging.Logger] = None
    ):
        """
        Initialize the ChEMBL data source.
        
        Args:
            db_operations: Database operations instance
            checkpoint_manager: Checkpoint manager for resumable imports
            progress_tracker: Progress tracker for monitoring import progress
            validator: Validator for molecular data
            config: Configuration dictionary
            logger: Logger instance
        """
        super().__init__(
            db_operations=db_operations,
            checkpoint_manager=checkpoint_manager,
            progress_tracker=progress_tracker,
            validator=validator,
            config=config,
            logger=logger
        )
        
        # ChEMBL-specific configuration
        self.api_base_url = self.config.get(
            'chembl_api_url',
            'https://www.ebi.ac.uk/chembl/api/data'
        )
        self.timeout = self.config.get('timeout', 30)

        # HTTP session
        self._session = None

        # Regex for ChEMBL ID validation
        self.chembl_id_pattern = re.compile(r'^CHEMBL\d+$')

        # Initialize ChEMBL client if available
        if CHEMBL_CLIENT_AVAILABLE:
            self.molecule_client = new_client.molecule
            self.logger.info("Using ChEMBL Web Resources client for improved performance")
        else:
            self.molecule_client = None
            self.logger.warning("ChEMBL Web Resources client not available. Install with: pip install chembl_webresource_client")

        # Molecular property filters for identifying potential cryoprotectants
        self.property_filters = self.config.get('property_filters', [
            # Small molecules with hydrogen bonding capability
            {
                'full_mwt__lte': 200.0,  # Molecular weight less than or equal to 200
                'hba__gte': 2,           # At least 2 hydrogen bond acceptors
                'hbd__gte': 1            # At least 1 hydrogen bond donor
            },
            # Medium-sized molecules with high hydrogen bonding capability
            {
                'full_mwt__lte': 350.0,  # Molecular weight less than or equal to 350
                'full_mwt__gte': 150.0,  # Molecular weight greater than or equal to 150
                'hba__gte': 3,           # At least 3 hydrogen bond acceptors
                'hbd__gte': 2            # At least 2 hydrogen bond donors
            },
            # Polyols and similar compounds
            {
                'full_mwt__lte': 400.0,  # Molecular weight less than or equal to 400
                'hba_lipinski__gte': 4,  # At least 4 hydrogen bond acceptors (Lipinski)
                'hbd_lipinski__gte': 3   # At least 3 hydrogen bond donors (Lipinski)
            },
            # Compounds with good water solubility (based on LogP)
            {
                'alogp__lte': 1.0,       # LogP less than or equal to 1.0
                'full_mwt__lte': 300.0,  # Molecular weight less than or equal to 300
                'hba__gte': 2            # At least 2 hydrogen bond acceptors
            }
        ])

        # Standard reference compounds to ensure they are always imported
        self.reference_compounds = self.config.get('reference_compounds', [
            "CHEMBL25",    # Aspirin
            "CHEMBL1118",  # Caffeine
            "CHEMBL1234",  # Glycerol (common cryoprotectant)
            "CHEMBL444",   # Glucose
            "CHEMBL230130", # Ethylene glycol (common cryoprotectant)
            "CHEMBL9335",  # Dimethyl sulfoxide (DMSO, common cryoprotectant)
            "CHEMBL15151"  # Trehalose (common cryoprotectant)
        ])

        # Cache for client responses
        self.cache = {}
        self.last_request_time = time.time()
    
    async def _get_session(self) -> aiohttp.ClientSession:
        """
        Get or create the HTTP session.
        
        Returns:
            aiohttp.ClientSession instance
        """
        if self._session is None or self._session.closed:
            self._session = aiohttp.ClientSession(
                timeout=aiohttp.ClientTimeout(total=self.timeout)
            )
        return self._session
    
    async def _close_session(self) -> None:
        """Close the HTTP session if open."""
        if self._session and not self._session.closed:
            await self._session.close()
            self._session = None
    
    async def _make_api_request(
        self,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        retries: int = 3
    ) -> Any:
        """
        Make a request to the ChEMBL API.

        Args:
            endpoint: API endpoint path
            params: Query parameters
            retries: Number of retry attempts

        Returns:
            Parsed response data
        """
        session = await self._get_session()
        url = f"{self.api_base_url}/{endpoint}"

        # Add headers to request JSON content
        headers = {
            'Accept': 'application/json',
            'Content-Type': 'application/json'
        }

        last_error = None
        for attempt in range(1, retries + 1):
            try:
                async with session.get(url, params=params, headers=headers) as response:
                    if response.status == 200:
                        # Check content type
                        content_type = response.headers.get('Content-Type', '')
                        if 'application/json' in content_type or 'json' in content_type.lower():
                            return await response.json()
                        elif 'application/xml' in content_type or 'xml' in content_type.lower():
                            # Try to extract JSON from XML response - for backwards compatibility
                            text_response = await response.text()
                            self.logger.warning(
                                f"Received XML instead of JSON from ChEMBL API. "
                                f"Will try to extract JSON content."
                            )

                            # Try to find JSON content embedded in the XML
                            try:
                                # Look for JSON-like structures
                                json_start = text_response.find('{')
                                json_end = text_response.rfind('}')

                                if json_start >= 0 and json_end > json_start:
                                    potential_json = text_response[json_start:json_end+1]
                                    return json.loads(potential_json)
                                else:
                                    raise ValueError("No JSON content found in XML response")
                            except Exception as json_error:
                                self.logger.error(f"Could not extract JSON from XML: {str(json_error)}")
                                return None
                        else:
                            # Try to parse as JSON anyway
                            try:
                                return await response.json()
                            except Exception:
                                # If JSON parsing fails, return text
                                self.logger.warning(
                                    f"Unexpected content type: {content_type}. "
                                    f"Trying to parse as text."
                                )
                                return None
                    elif response.status == 404:
                        return None
                    elif response.status == 429:
                        # Rate limited, wait longer
                        retry_after = int(response.headers.get('Retry-After', 1))
                        self.logger.warning(
                            f"Rate limited by ChEMBL API, waiting {retry_after} seconds"
                        )
                        await asyncio.sleep(retry_after)
                        continue
                    else:
                        response.raise_for_status()
            except (aiohttp.ClientError, asyncio.TimeoutError) as e:
                last_error = e
                if attempt < retries:
                    retry_delay = self.retry_delay * (2 ** (attempt - 1))  # Exponential backoff
                    self.logger.warning(
                        f"API request failed, retrying in {retry_delay:.2f}s "
                        f"(attempt {attempt}/{retries}): {str(e)}"
                    )
                    await asyncio.sleep(retry_delay)
                else:
                    self.logger.error(f"API request failed after {retries} attempts: {str(e)}")
                    raise

        if last_error:
            raise last_error

        return None
    
    def get_source_name(self) -> str:
        """
        Get the name of this data source.

        Returns:
            Source name
        """
        return "ChEMBL"
    
    async def _apply_rate_limiting(self) -> None:
        """Apply rate limiting to avoid overwhelming the ChEMBL API."""
        now = time.time()
        time_since_last_request = now - self.last_request_time

        if time_since_last_request < self.api_delay:
            await asyncio.sleep(self.api_delay - time_since_last_request)

        self.last_request_time = time.time()

    async def fetch_compound(self, identifier: str) -> Optional[Dict[str, Any]]:
        """
        Fetch a single compound by ChEMBL ID.

        Args:
            identifier: ChEMBL ID (e.g., CHEMBL123)

        Returns:
            Compound data dictionary or None if not found
        """
        try:
            # Validate ChEMBL ID format
            if not self.chembl_id_pattern.match(identifier):
                self.logger.error(f"Invalid ChEMBL ID format: {identifier}")
                return None

            # Check cache first
            cache_key = f"chembl_{identifier}"
            if cache_key in self.cache:
                return self.cache[cache_key]

            # Use the client if available for improved performance
            if CHEMBL_CLIENT_AVAILABLE and self.molecule_client:
                try:
                    # Apply rate limiting
                    await self._apply_rate_limiting()

                    # Use asyncio.to_thread to run the blocking client in a thread
                    data = await asyncio.to_thread(self.molecule_client.get, identifier)

                    if data:
                        # Cache the result
                        self.cache[cache_key] = data
                        return data
                except Exception as client_error:
                    self.logger.warning(f"ChEMBL client error for {identifier}, falling back to API: {str(client_error)}")

            # Fall back to API request
            endpoint = f"molecule/{identifier}"
            data = await self._make_api_request(endpoint)

            if data:
                # Cache the result
                self.cache[cache_key] = data

                # Report finding to logger
                if 'molecule_chembl_id' in data:
                    self.logger.info(f"Found compound {data['molecule_chembl_id']} - {data.get('pref_name', 'Unnamed')}")

            return data
        except Exception as e:
            self.logger.error(f"Error fetching compound {identifier}: {str(e)}")
            return None
    
    async def search_compounds(
        self,
        query: str,
        max_results: Optional[int] = None
    ) -> List[str]:
        """
        Search for compounds matching a query.

        Args:
            query: Search query string
            max_results: Maximum number of results to return

        Returns:
            List of ChEMBL IDs matching the query
        """
        try:
            if query == "*":
                # Special case for "all compounds"
                return await self._get_all_compounds(max_results)

            # Try to use the client first if available
            if CHEMBL_CLIENT_AVAILABLE and self.molecule_client:
                try:
                    # Apply rate limiting
                    await self._apply_rate_limiting()

                    # Use asyncio.to_thread to run the blocking client in a thread
                    # Use pref_name__icontains for broader matching
                    results = await asyncio.to_thread(
                        self.molecule_client.filter,
                        pref_name__icontains=query
                    )

                    # Extract ChEMBL IDs
                    chembl_ids = [result['molecule_chembl_id'] for result in results
                                if 'molecule_chembl_id' in result]

                    # Apply limit if specified
                    if max_results and len(chembl_ids) > max_results:
                        chembl_ids = chembl_ids[:max_results]

                    self.logger.info(f"Found {len(chembl_ids)} compounds matching query '{query}' using client")
                    return chembl_ids

                except Exception as client_error:
                    self.logger.warning(f"ChEMBL client search error, falling back to API: {str(client_error)}")

            # Fall back to API
            # Prepare query parameters
            params = {
                'search': query,
                'limit': min(max_results or 1000, 1000)  # ChEMBL API limit
            }

            # Fetch the search results
            endpoint = "molecule"
            data = await self._make_api_request(endpoint, params)

            if not data or 'molecules' not in data:
                return []

            # Extract ChEMBL IDs
            chembl_ids = [molecule['molecule_chembl_id']
                         for molecule in data['molecules']]

            # Apply limit if specified
            if max_results and len(chembl_ids) > max_results:
                chembl_ids = chembl_ids[:max_results]

            self.logger.info(f"Found {len(chembl_ids)} compounds matching query '{query}' using API")
            return chembl_ids
        except Exception as e:
            self.logger.error(f"Error searching compounds: {str(e)}")
            return []
    
    async def _get_all_compounds(self, max_results: Optional[int] = None) -> List[str]:
        """
        Get all compounds from ChEMBL.

        This method is optimized to use the property filters and reference compounds
        to fetch a reasonable subset of compounds rather than the entire database.

        Args:
            max_results: Maximum number of results to return

        Returns:
            List of ChEMBL IDs
        """
        try:
            # Start with reference compounds
            chembl_ids = set(self.reference_compounds)
            self.logger.info(f"Added {len(chembl_ids)} reference compounds")

            if CHEMBL_CLIENT_AVAILABLE and self.molecule_client:
                # Use property-based filtering with the client
                for i, prop_filter in enumerate(self.property_filters):
                    if max_results and len(chembl_ids) >= max_results:
                        break

                    try:
                        self.logger.info(f"Applying property filter #{i+1}: {prop_filter}")

                        # Apply rate limiting
                        await self._apply_rate_limiting()

                        # Use asyncio.to_thread to run the blocking client in a thread
                        results = await asyncio.to_thread(
                            self.molecule_client.filter,
                            **prop_filter
                        )

                        # Extract ChEMBL IDs
                        for result in results:
                            if 'molecule_chembl_id' in result:
                                chembl_ids.add(result['molecule_chembl_id'])

                                # Check limit
                                if max_results and len(chembl_ids) >= max_results:
                                    break

                        self.logger.info(f"Found {len(chembl_ids)} compounds after filter #{i+1}")

                    except Exception as e:
                        self.logger.error(f"Error applying property filter #{i+1}: {str(e)}")
            else:
                # Use API-based paging as fallback
                # Get the total count
                total_count = await self.get_compound_count()

                if total_count == 0:
                    return []

                # Apply limit if specified
                limit = min(max_results or total_count, total_count)

                # Collect all IDs
                page_size = 1000  # ChEMBL API limit

                for offset in range(0, limit, page_size):
                    # Calculate current page size
                    current_limit = min(page_size, limit - offset)

                    # Fetch the page
                    params = {
                        'limit': current_limit,
                        'offset': offset
                    }

                    endpoint = "molecule"
                    data = await self._make_api_request(endpoint, params)

                    if not data or 'molecules' not in data:
                        break

                    # Extract ChEMBL IDs
                    page_ids = [molecule['molecule_chembl_id']
                               for molecule in data['molecules']]
                    chembl_ids.update(page_ids)

                    if len(page_ids) < current_limit:
                        # Fewer results than requested, we've reached the end
                        break

                    # Rate limiting
                    await asyncio.sleep(self.api_delay)

            # Convert set to list
            id_list = list(chembl_ids)

            # Apply limit if specified
            if max_results and len(id_list) > max_results:
                id_list = id_list[:max_results]

            self.logger.info(f"Returning {len(id_list)} ChEMBL IDs")
            return id_list
        except Exception as e:
            self.logger.error(f"Error getting all compounds: {str(e)}")
            return []
    
    async def get_compound_batch(
        self,
        identifiers: List[str]
    ) -> List[Dict[str, Any]]:
        """
        Fetch multiple compounds by ChEMBL ID.
        
        Args:
            identifiers: List of ChEMBL IDs
            
        Returns:
            List of compound data dictionaries
        """
        try:
            # ChEMBL API doesn't support batch fetching directly
            # We'll fetch compounds individually but concurrently
            
            # Create tasks for each compound
            tasks = []
            for chembl_id in identifiers:
                tasks.append(self.fetch_compound(chembl_id))
            
            # Fetch all compounds concurrently
            results = await asyncio.gather(*tasks, return_exceptions=True)
            
            # Process results
            compounds = []
            for i, result in enumerate(results):
                if isinstance(result, Exception):
                    self.logger.error(
                        f"Error fetching compound {identifiers[i]}: {str(result)}"
                    )
                elif result is not None:
                    compounds.append(result)
            
            return compounds
        except Exception as e:
            self.logger.error(f"Error fetching compound batch: {str(e)}")
            return []
    
    async def get_property_data(
        self,
        compound_data: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """
        Fetch property data for a compound.

        Args:
            compound_data: Compound data dictionary

        Returns:
            List of property data dictionaries
        """
        properties = []

        try:
            # Extract ChEMBL ID
            chembl_id = compound_data.get('molecule_chembl_id')

            if not chembl_id:
                return properties

            # Extract basic properties from compound data
            if 'molecule_properties' in compound_data:
                mol_props = compound_data['molecule_properties']

                # Process each property
                for prop_name, value in mol_props.items():
                    if value is not None:
                        # Skip null values

                        # Determine value type and property type
                        value_type = self._determine_value_type(value)
                        prop_type = self._determine_property_type(prop_name)

                        # Create property record
                        property_data = {
                            'property_name': prop_name,
                            'property_type': prop_type,
                            f'{value_type}_value': value,
                            'source': 'ChEMBL'
                        }

                        properties.append(property_data)

            # Fetch additional properties
            await self._add_compound_properties(chembl_id, properties)

            return properties
        except Exception as e:
            self.logger.error(f"Error getting property data: {str(e)}")
            return []
    
    async def _add_compound_properties(
        self,
        chembl_id: str,
        properties: List[Dict[str, Any]]
    ) -> None:
        """
        Add additional properties for a compound.
        
        Args:
            chembl_id: ChEMBL ID
            properties: List to add properties to
        """
        try:
            # Fetch drug properties
            endpoint = f"mechanism?molecule_chembl_id={chembl_id}"
            mechanism_data = await self._make_api_request(endpoint)
            
            if mechanism_data and 'mechanisms' in mechanism_data:
                for mechanism in mechanism_data['mechanisms']:
                    # Extract mechanism of action
                    if 'mechanism_of_action' in mechanism and mechanism['mechanism_of_action']:
                        properties.append({
                            'property_name': 'mechanism_of_action',
                            'property_type': 'biological',
                            'text_value': mechanism['mechanism_of_action'],
                            'source': 'ChEMBL'
                        })
                    
                    # Extract target information
                    if 'target_chembl_id' in mechanism and mechanism['target_chembl_id']:
                        properties.append({
                            'property_name': 'target_chembl_id',
                            'property_type': 'biological',
                            'text_value': mechanism['target_chembl_id'],
                            'source': 'ChEMBL'
                        })
            
            # Fetch activity data (limited to a few examples)
            endpoint = f"activity?molecule_chembl_id={chembl_id}&limit=5"
            activity_data = await self._make_api_request(endpoint)
            
            if activity_data and 'activities' in activity_data:
                for activity in activity_data['activities']:
                    # Extract activity value
                    if ('standard_value' in activity and 
                        'standard_units' in activity and 
                        'standard_type' in activity and
                        activity['standard_value'] is not None):
                        
                        value_type = self._determine_value_type(activity['standard_value'])
                        property_data = {
                            'property_name': activity['standard_type'],
                            'property_type': 'biological',
                            f'{value_type}_value': activity['standard_value'],
                            'source': 'ChEMBL'
                        }
                        
                        if activity['standard_units']:
                            property_data['unit'] = activity['standard_units']
                            
                        properties.append(property_data)
        except Exception as e:
            self.logger.warning(f"Error adding compound properties for {chembl_id}: {str(e)}")
    
    def _determine_property_type(self, property_name: str) -> str:
        """
        Determine the type of a property based on its name.
        
        Args:
            property_name: Name of the property
            
        Returns:
            Property type category
        """
        physicochemical_props = {
            'molecular_weight', 'alogp', 'psa', 'rtb', 'hba', 'hbd',
            'aromatic_rings', 'heavy_atoms', 'qed_weighted',
            'cx_logp', 'cx_logd', 'full_mwt', 'mw_freebase', 'mw_monoisotopic',
            'num_lipinski_ro5_violations', 'ro3_pass'
        }
        
        structural_props = {
            'molecule_type', 'full_molformula', 'structure_type',
            'inorganic_flag', 'natural_product', 'polymer_flag',
            'chirality', 'num_ro5_violations', 'med_chem_friendly',
            'molecule_structures', 'atc_classifications'
        }
        
        biological_props = {
            'activity', 'potency', 'ic50', 'ec50', 'inhibition',
            'ki', 'target', 'mechanism', 'pathway', 'indication'
        }
        
        # Check property name against categories
        property_name_clean = property_name.replace('_', '').lower()
        
        for prop in physicochemical_props:
            if prop.replace('_', '') in property_name_clean:
                return 'physicochemical'
                
        for prop in structural_props:
            if prop.replace('_', '') in property_name_clean:
                return 'structural'
                
        for prop in biological_props:
            if prop.replace('_', '') in property_name_clean:
                return 'biological'
        
        # Default category
        return 'other'
    
    def _determine_value_type(self, value: Any) -> str:
        """
        Determine the type of a property value.
        
        Args:
            value: Property value
            
        Returns:
            Value type (numeric, text, boolean)
        """
        if isinstance(value, (int, float)):
            return 'numeric'
        elif isinstance(value, bool):
            return 'boolean'
        else:
            # Try to convert to number
            try:
                float(value)
                return 'numeric'
            except (ValueError, TypeError):
                return 'text'
    
    async def get_compound_count(self, query: Optional[str] = None) -> int:
        """
        Get the total count of compounds matching a query.
        
        Args:
            query: Optional search query to filter compounds
            
        Returns:
            Total count of matched compounds
        """
        try:
            if query and query != "*":
                # Prepare query parameters
                params = {
                    'search': query,
                    'limit': 1  # We only need the count
                }
            else:
                # No query, get total count
                params = {
                    'limit': 1  # We only need the count
                }
            
            # Fetch the search results
            endpoint = "molecule"
            data = await self._make_api_request(endpoint, params)
            
            if data and 'page_meta' in data and 'total_count' in data['page_meta']:
                return data['page_meta']['total_count']
                
            return 0
        except Exception as e:
            self.logger.error(f"Error getting compound count: {str(e)}")
            return 0
    
    async def transform_compound_data(
        self,
        raw_data: Dict[str, Any]
    ) -> Dict[str, Any]:
        """
        Transform raw ChEMBL compound data to the unified format.

        Args:
            raw_data: Raw ChEMBL compound data

        Returns:
            Transformed compound data in unified format
        """
        # Extract structures
        structures = raw_data.get('molecule_structures', {}) or {}

        # Extract molecule properties
        mol_props = raw_data.get('molecule_properties', {}) or {}

        # Get ChEMBL version if available
        chembl_version = getattr(raw_data, 'chembl_version', None) or "Unknown"

        # Get the preferred name or a fallback
        name = raw_data.get('pref_name')
        if not name:
            # Try alternative name fields
            name = (raw_data.get('molecule_synonyms', [{}])[0].get('synonym') if raw_data.get('molecule_synonyms')
                    else raw_data.get('molecule_chembl_id', 'Unknown Compound'))

        # Calculate initial properties for JSONB field
        properties = {}
        for prop_key, prop_value in mol_props.items():
            if prop_value is not None:
                properties[prop_key] = prop_value

        # Ensure common properties are explicitly set if available in mol_props
        if mol_props.get('full_mwt'):
            properties['Molecular Weight'] = mol_props.get('full_mwt')
        if mol_props.get('alogp'):
            properties['LogP'] = mol_props.get('alogp')
        if mol_props.get('psa'):
            properties['TPSA'] = mol_props.get('psa')
        if mol_props.get('hba'):
            properties['Hydrogen Bond Acceptor Count'] = mol_props.get('hba')
        if mol_props.get('hbd'):
            properties['Hydrogen Bond Donor Count'] = mol_props.get('hbd')
        if mol_props.get('rtb'):
            properties['Rotatable Bond Count'] = mol_props.get('rtb')
        if mol_props.get('ro3_pass'):
            properties['Rule of Three Pass'] = mol_props.get('ro3_pass')
        if mol_props.get('num_ro5_violations'):
            properties['Number Rule of Five Violations'] = mol_props.get('num_ro5_violations')

        # Extract synonyms
        synonyms = []
        if raw_data.get('molecule_synonyms'):
            for synonym_data in raw_data.get('molecule_synonyms'):
                if synonym_data.get('synonym'):
                    synonyms.append(synonym_data.get('synonym'))

        # Filter and limit synonyms
        if synonyms:
            seen = set()
            filtered_synonyms = []

            for synonym in synonyms:
                if synonym and synonym not in seen:
                    filtered_synonyms.append(synonym)
                    seen.add(synonym)

            # Limit to a reasonable number
            max_synonyms = 20
            if len(filtered_synonyms) > max_synonyms:
                filtered_synonyms = filtered_synonyms[:max_synonyms]

            synonyms = filtered_synonyms

        # Transform to the unified format
        transformed = {
            "id": str(uuid.uuid4()),  # Generate a unique ID for the molecule
            "name": name,
            "smiles": structures.get('canonical_smiles'),
            "inchi": structures.get('standard_inchi'),
            "inchikey": structures.get('standard_inchi_key'),
            "formula": mol_props.get('full_molformula'),
            "molecular_weight": mol_props.get('full_mwt'),
            "pubchem_cid": None,  # Will be populated later by cross-reference resolution
            "chembl_id": raw_data.get('molecule_chembl_id'),
            "source": f"ChEMBL v{chembl_version}",
            "source_id": raw_data.get('molecule_chembl_id'),
            "properties": properties,
            "synonyms": synonyms,
            "metadata": {
                "chembl_version": chembl_version,
                "first_approval": raw_data.get('first_approval'),
                "oral": raw_data.get('oral'),
                "topical": raw_data.get('topical'),
                "parenteral": raw_data.get('parenteral'),
                "black_box_warning": raw_data.get('black_box_warning'),
                "availability_type": raw_data.get('availability_type'),
                "chirality": raw_data.get('chirality'),
                "max_phase": raw_data.get('max_phase')
            }
        }

        return transformed
    
    async def transform_property_data(
        self,
        raw_properties: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """
        Transform raw property data to the unified format.
        
        Args:
            raw_properties: List of raw property dictionaries
            
        Returns:
            List of transformed property dictionaries
        """
        # In this case, the properties are already in the unified format
        # Return as is, but ensure source is set
        for prop in raw_properties:
            if 'source' not in prop:
                prop['source'] = 'ChEMBL'
        
        return raw_properties
    
    async def import_all(
        self,
        max_results: Optional[int] = None
    ) -> Tuple[int, int, List[Tuple[str, str]]]:
        """
        Import all compounds from ChEMBL based on property filters and reference compounds.

        Instead of importing the entire ChEMBL database (which is very large),
        this method uses property filters to identify potential cryoprotectants
        and other compounds of interest.

        Args:
            max_results: Maximum number of results to import

        Returns:
            Tuple of (success_count, failure_count, failures_with_reasons)
        """
        self.logger.info(f"Importing compounds from ChEMBL using optimized property-based approach")

        # Get compounds based on property filters and reference compounds
        identifiers = await self._get_all_compounds(max_results)

        # Set total items in progress tracker
        if self.progress_tracker:
            self.progress_tracker.total_items = len(identifiers)

        # Import the compounds
        return await self.import_compounds(identifiers)

    async def stream_compound_identifiers(
        self,
        query: Optional[str] = None,
        batch_size: int = 1000
    ) -> AsyncIterator[List[str]]:
        """
        Stream compound identifiers in batches.

        Args:
            query: Optional query string to filter compounds
            batch_size: Number of identifiers to yield in each batch

        Yields:
            Batches of compound identifiers
        """
        try:
            # If no query or special "*" query, use optimized method
            if not query or query == "*":
                # First yield reference compounds
                if self.reference_compounds:
                    yield self.reference_compounds

                # Then yield results from property filters if client is available
                if CHEMBL_CLIENT_AVAILABLE and self.molecule_client:
                    for prop_filter in self.property_filters:
                        try:
                            self.logger.info(f"Streaming compounds with property filter: {prop_filter}")

                            # Apply rate limiting
                            await self._apply_rate_limiting()

                            # Use asyncio.to_thread to run the blocking client in a thread
                            results = await asyncio.to_thread(
                                self.molecule_client.filter,
                                **prop_filter
                            )

                            # Extract ChEMBL IDs
                            ids_batch = []
                            for result in results:
                                if 'molecule_chembl_id' in result:
                                    ids_batch.append(result['molecule_chembl_id'])

                                    # Yield when batch is full
                                    if len(ids_batch) >= batch_size:
                                        yield ids_batch
                                        ids_batch = []

                            # Yield any remaining IDs
                            if ids_batch:
                                yield ids_batch

                        except Exception as e:
                            self.logger.error(f"Error applying property filter: {str(e)}")
                else:
                    # Fall back to standard API pagination
                    # Prepare base parameters
                    params = {}

                    # Get total count
                    total_count = await self.get_compound_count()

                    if total_count == 0:
                        return

                    # Stream in batches
                    for offset in range(0, total_count, batch_size):
                        # Add pagination parameters
                        current_params = {**params, 'limit': min(batch_size, total_count - offset), 'offset': offset}

                        # Fetch the page
                        endpoint = "molecule"
                        data = await self._make_api_request(endpoint, current_params)

                        if not data or 'molecules' not in data:
                            break

                        # Extract ChEMBL IDs
                        chembl_ids = [molecule['molecule_chembl_id']
                                    for molecule in data['molecules']]

                        yield chembl_ids

                        if len(chembl_ids) < current_params['limit']:
                            # Fewer results than requested, we've reached the end
                            break

                        # Rate limiting
                        await asyncio.sleep(self.api_delay)
            else:
                # Specific query - use regular search
                # Try to use the client first if available
                if CHEMBL_CLIENT_AVAILABLE and self.molecule_client:
                    try:
                        # Apply rate limiting
                        await self._apply_rate_limiting()

                        # Use the client for search
                        results = await asyncio.to_thread(
                            self.molecule_client.filter,
                            pref_name__icontains=query
                        )

                        # Extract IDs and yield in batches
                        chembl_ids = []
                        for result in results:
                            if 'molecule_chembl_id' in result:
                                chembl_ids.append(result['molecule_chembl_id'])

                                if len(chembl_ids) >= batch_size:
                                    yield chembl_ids
                                    chembl_ids = []

                        # Yield any remaining IDs
                        if chembl_ids:
                            yield chembl_ids

                        return
                    except Exception as client_error:
                        self.logger.warning(f"Client search failed, using API: {str(client_error)}")

                # Fall back to API search
                # Prepare base parameters
                params = {'search': query}

                # Get total count
                total_count = await self.get_compound_count(query)

                if total_count == 0:
                    return

                # Stream in batches
                for offset in range(0, total_count, batch_size):
                    # Add pagination parameters
                    current_params = {**params, 'limit': min(batch_size, total_count - offset), 'offset': offset}

                    # Fetch the page
                    endpoint = "molecule"
                    data = await self._make_api_request(endpoint, current_params)

                    if not data or 'molecules' not in data:
                        break

                    # Extract ChEMBL IDs
                    chembl_ids = [molecule['molecule_chembl_id']
                                for molecule in data['molecules']]

                    yield chembl_ids

                    if len(chembl_ids) < current_params['limit']:
                        # Fewer results than requested, we've reached the end
                        break

                    # Rate limiting
                    await asyncio.sleep(self.api_delay)
        except Exception as e:
            self.logger.error(f"Error streaming compound identifiers: {str(e)}")
            yield []