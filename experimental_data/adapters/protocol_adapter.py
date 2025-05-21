#!/usr/bin/env python3
"""
Protocol Adapter for CryoProtect Enhanced Experimental Data System.

This module provides an adapter to integrate between Convex-based protocols 
and Flask API-based experimental data system.
"""

from typing import Dict, List, Any, Optional, Union, Tuple, TypedDict
import uuid
from datetime import datetime
import aiohttp
import json
import logging
import os
from pydantic import BaseModel, Field

from ..models import Protocol, ProtocolStep, ValidationError, JSONDict
from ..utils.authentication import get_auth_headers

# Configure logging
logger = logging.getLogger(__name__)

class ConvexProtocolDTO(BaseModel):
    """Data Transfer Object for Convex Protocol data."""
    id: str
    name: str
    description: Optional[str] = None
    steps: List[Dict[str, Any]]
    version: Optional[int] = None
    parent_id: Optional[str] = None
    category: Optional[str] = None
    is_template: Optional[bool] = None
    parameters: Optional[Dict[str, Any]] = None
    created_by: Optional[str] = None
    created_at: Optional[int] = None
    updated_at: Optional[int] = None
    public: Optional[bool] = None

class ProtocolAdapter:
    """Adapter for integrating between Convex protocols and experimental data system."""
    
    def __init__(self, convex_api_url: Optional[str] = None):
        """
        Initialize protocol adapter.
        
        Args:
            convex_api_url: URL of the Convex API. If not provided, it will be read
                           from the CONVEX_API_URL environment variable.
        """
        self.convex_api_url = convex_api_url or os.environ.get("CONVEX_API_URL")
        if not self.convex_api_url:
            raise ValueError("Convex API URL is required. Set CONVEX_API_URL environment variable or provide it to the constructor.")
        
        # Ensure URL ends with a slash
        if not self.convex_api_url.endswith("/"):
            self.convex_api_url += "/"
    
    async def get_protocol(self, protocol_id: str) -> Protocol:
        """
        Get a protocol from Convex API and convert to experimental data format.
        
        Args:
            protocol_id: ID of the protocol in Convex
            
        Returns:
            Protocol in experimental data format
            
        Raises:
            ValidationError: If the protocol is not found or invalid
        """
        try:
            # Get auth headers
            headers = await get_auth_headers()
            
            # Make API request to Convex
            async with aiohttp.ClientSession() as session:
                url = f"{self.convex_api_url}api/protocols/{protocol_id}"
                async with session.get(url, headers=headers) as response:
                    if response.status != 200:
                        error_text = await response.text()
                        raise ValidationError(f"Failed to get protocol: {error_text}")
                    
                    convex_protocol = await response.json()
            
            # Convert to experimental data format
            return self._convert_from_convex(convex_protocol)
        
        except aiohttp.ClientError as e:
            logger.error(f"HTTP error while getting protocol {protocol_id}: {str(e)}")
            raise ValidationError(f"Failed to get protocol: {str(e)}")
        
        except Exception as e:
            logger.error(f"Error while getting protocol {protocol_id}: {str(e)}")
            raise ValidationError(f"Failed to get protocol: {str(e)}")
    
    async def list_protocols(
        self, 
        filters: Optional[Dict[str, Any]] = None,
        page: int = 1,
        page_size: int = 20,
        sort_by: str = "updatedAt",
        sort_order: str = "desc"
    ) -> Tuple[List[Protocol], int]:
        """
        List protocols from Convex API and convert to experimental data format.
        
        Args:
            filters: Filters to apply (name, category, isTemplate, etc.)
            page: Page number (1-indexed)
            page_size: Number of items per page
            sort_by: Field to sort by
            sort_order: Sort order (asc or desc)
            
        Returns:
            Tuple of (protocols, total_count)
            
        Raises:
            ValidationError: If the request fails
        """
        try:
            # Get auth headers
            headers = await get_auth_headers()
            
            # Build query parameters
            params = {
                "page": page,
                "per_page": page_size,
                "sort_by": sort_by,
                "sort_order": sort_order
            }
            
            # Add filters if provided
            if filters:
                for key, value in filters.items():
                    params[key] = value
            
            # Make API request to Convex
            async with aiohttp.ClientSession() as session:
                url = f"{self.convex_api_url}api/protocols"
                async with session.get(url, params=params, headers=headers) as response:
                    if response.status != 200:
                        error_text = await response.text()
                        raise ValidationError(f"Failed to list protocols: {error_text}")
                    
                    result = await response.json()
                    convex_protocols = result.get("data", [])
                    total_count = result.get("total", 0)
            
            # Convert each protocol
            protocols = [self._convert_from_convex(p) for p in convex_protocols]
            
            return protocols, total_count
        
        except aiohttp.ClientError as e:
            logger.error(f"HTTP error while listing protocols: {str(e)}")
            raise ValidationError(f"Failed to list protocols: {str(e)}")
        
        except Exception as e:
            logger.error(f"Error while listing protocols: {str(e)}")
            raise ValidationError(f"Failed to list protocols: {str(e)}")
    
    async def create_protocol(self, protocol: Protocol) -> str:
        """
        Create a protocol in Convex API from experimental data format.
        
        Args:
            protocol: Protocol in experimental data format
            
        Returns:
            ID of the created protocol in Convex
            
        Raises:
            ValidationError: If the creation fails
        """
        try:
            # Convert to Convex format
            convex_protocol = self._convert_to_convex(protocol)
            
            # Get auth headers
            headers = await get_auth_headers()
            
            # Make API request to Convex
            async with aiohttp.ClientSession() as session:
                url = f"{self.convex_api_url}api/protocols"
                async with session.post(url, json=convex_protocol, headers=headers) as response:
                    if response.status != 201:
                        error_text = await response.text()
                        raise ValidationError(f"Failed to create protocol: {error_text}")
                    
                    result = await response.json()
                    return result.get("id")
        
        except aiohttp.ClientError as e:
            logger.error(f"HTTP error while creating protocol: {str(e)}")
            raise ValidationError(f"Failed to create protocol: {str(e)}")
        
        except Exception as e:
            logger.error(f"Error while creating protocol: {str(e)}")
            raise ValidationError(f"Failed to create protocol: {str(e)}")
    
    async def update_protocol(self, protocol_id: str, protocol: Protocol) -> str:
        """
        Update a protocol in Convex API from experimental data format.
        
        Args:
            protocol_id: ID of the protocol in Convex
            protocol: Protocol in experimental data format
            
        Returns:
            ID of the updated protocol in Convex
            
        Raises:
            ValidationError: If the update fails
        """
        try:
            # Convert to Convex format
            convex_protocol = self._convert_to_convex(protocol)
            
            # Get auth headers
            headers = await get_auth_headers()
            
            # Make API request to Convex
            async with aiohttp.ClientSession() as session:
                url = f"{self.convex_api_url}api/protocols/{protocol_id}"
                async with session.patch(url, json=convex_protocol, headers=headers) as response:
                    if response.status != 200:
                        error_text = await response.text()
                        raise ValidationError(f"Failed to update protocol: {error_text}")
                    
                    result = await response.json()
                    return result.get("id")
        
        except aiohttp.ClientError as e:
            logger.error(f"HTTP error while updating protocol: {str(e)}")
            raise ValidationError(f"Failed to update protocol: {str(e)}")
        
        except Exception as e:
            logger.error(f"Error while updating protocol: {str(e)}")
            raise ValidationError(f"Failed to update protocol: {str(e)}")
    
    async def delete_protocol(self, protocol_id: str) -> bool:
        """
        Delete a protocol in Convex API.
        
        Args:
            protocol_id: ID of the protocol in Convex
            
        Returns:
            True if deletion was successful
            
        Raises:
            ValidationError: If the deletion fails
        """
        try:
            # Get auth headers
            headers = await get_auth_headers()
            
            # Make API request to Convex
            async with aiohttp.ClientSession() as session:
                url = f"{self.convex_api_url}api/protocols/{protocol_id}"
                async with session.delete(url, headers=headers) as response:
                    if response.status != 200:
                        error_text = await response.text()
                        raise ValidationError(f"Failed to delete protocol: {error_text}")
                    
                    return True
        
        except aiohttp.ClientError as e:
            logger.error(f"HTTP error while deleting protocol: {str(e)}")
            raise ValidationError(f"Failed to delete protocol: {str(e)}")
        
        except Exception as e:
            logger.error(f"Error while deleting protocol: {str(e)}")
            raise ValidationError(f"Failed to delete protocol: {str(e)}")
    
    async def create_protocol_version(self, protocol_id: str, changes: Dict[str, Any] = None) -> str:
        """
        Create a new version of a protocol in Convex API.
        
        Args:
            protocol_id: ID of the protocol in Convex
            changes: Changes to apply to the new version
            
        Returns:
            ID of the new protocol version in Convex
            
        Raises:
            ValidationError: If the creation fails
        """
        try:
            # Get auth headers
            headers = await get_auth_headers()
            
            # Prepare request data
            request_data = {
                "protocol": {
                    "parentId": protocol_id
                }
            }
            
            # Add changes if provided
            if changes:
                for key, value in changes.items():
                    request_data["protocol"][key] = value
            
            # Make API request to Convex
            async with aiohttp.ClientSession() as session:
                url = f"{self.convex_api_url}api/protocols/{protocol_id}/versions"
                async with session.post(url, json=request_data, headers=headers) as response:
                    if response.status != 201:
                        error_text = await response.text()
                        raise ValidationError(f"Failed to create protocol version: {error_text}")
                    
                    result = await response.json()
                    return result.get("id")
        
        except aiohttp.ClientError as e:
            logger.error(f"HTTP error while creating protocol version: {str(e)}")
            raise ValidationError(f"Failed to create protocol version: {str(e)}")
        
        except Exception as e:
            logger.error(f"Error while creating protocol version: {str(e)}")
            raise ValidationError(f"Failed to create protocol version: {str(e)}")
    
    async def link_protocol_to_experiment(self, protocol_id: str, experiment_id: str) -> Dict[str, Any]:
        """
        Link a protocol to an experiment for experimental data tracking.
        
        Args:
            protocol_id: ID of the protocol in Convex
            experiment_id: ID of the experiment
            
        Returns:
            Details of the link created
            
        Raises:
            ValidationError: If the linking fails
        """
        try:
            # Get auth headers
            headers = await get_auth_headers()
            
            # Prepare request data
            request_data = {
                "protocol_id": protocol_id,
                "experiment_id": experiment_id,
                "linked_at": int(datetime.now().timestamp() * 1000),
            }
            
            # Make API request to Convex
            async with aiohttp.ClientSession() as session:
                url = f"{self.convex_api_url}api/protocol-experiment-links"
                async with session.post(url, json=request_data, headers=headers) as response:
                    if response.status != 201:
                        error_text = await response.text()
                        raise ValidationError(f"Failed to link protocol to experiment: {error_text}")
                    
                    return await response.json()
        
        except aiohttp.ClientError as e:
            logger.error(f"HTTP error while linking protocol to experiment: {str(e)}")
            raise ValidationError(f"Failed to link protocol to experiment: {str(e)}")
        
        except Exception as e:
            logger.error(f"Error while linking protocol to experiment: {str(e)}")
            raise ValidationError(f"Failed to link protocol to experiment: {str(e)}")
    
    async def get_experiments_for_protocol(self, protocol_id: str) -> List[Dict[str, Any]]:
        """
        Get experiments linked to a protocol.
        
        Args:
            protocol_id: ID of the protocol in Convex
            
        Returns:
            List of experiment IDs linked to the protocol
            
        Raises:
            ValidationError: If the retrieval fails
        """
        try:
            # Get auth headers
            headers = await get_auth_headers()
            
            # Make API request to Convex
            async with aiohttp.ClientSession() as session:
                url = f"{self.convex_api_url}api/protocol-experiment-links"
                params = {"protocol_id": protocol_id}
                async with session.get(url, params=params, headers=headers) as response:
                    if response.status != 200:
                        error_text = await response.text()
                        raise ValidationError(f"Failed to get experiments for protocol: {error_text}")
                    
                    result = await response.json()
                    return result.get("data", [])
        
        except aiohttp.ClientError as e:
            logger.error(f"HTTP error while getting experiments for protocol: {str(e)}")
            raise ValidationError(f"Failed to get experiments for protocol: {str(e)}")
        
        except Exception as e:
            logger.error(f"Error while getting experiments for protocol: {str(e)}")
            raise ValidationError(f"Failed to get experiments for protocol: {str(e)}")
    
    async def get_protocols_for_experiment(self, experiment_id: str) -> List[Dict[str, Any]]:
        """
        Get protocols linked to an experiment.
        
        Args:
            experiment_id: ID of the experiment
            
        Returns:
            List of protocol IDs linked to the experiment
            
        Raises:
            ValidationError: If the retrieval fails
        """
        try:
            # Get auth headers
            headers = await get_auth_headers()
            
            # Make API request to Convex
            async with aiohttp.ClientSession() as session:
                url = f"{self.convex_api_url}api/protocol-experiment-links"
                params = {"experiment_id": experiment_id}
                async with session.get(url, params=params, headers=headers) as response:
                    if response.status != 200:
                        error_text = await response.text()
                        raise ValidationError(f"Failed to get protocols for experiment: {error_text}")
                    
                    result = await response.json()
                    return result.get("data", [])
        
        except aiohttp.ClientError as e:
            logger.error(f"HTTP error while getting protocols for experiment: {str(e)}")
            raise ValidationError(f"Failed to get protocols for experiment: {str(e)}")
        
        except Exception as e:
            logger.error(f"Error while getting protocols for experiment: {str(e)}")
            raise ValidationError(f"Failed to get protocols for experiment: {str(e)}")
    
    def _convert_from_convex(self, convex_protocol: Dict[str, Any]) -> Protocol:
        """
        Convert a protocol from Convex format to experimental data format.
        
        Args:
            convex_protocol: Protocol in Convex format
            
        Returns:
            Protocol in experimental data format
        """
        # Validate with Pydantic model
        validated = ConvexProtocolDTO(**convex_protocol)
        
        # Map template_type from category
        template_type = "standard"
        if validated.category:
            template_type = validated.category
        
        # Map date fields
        created_at = None
        if validated.created_at:
            created_at = datetime.fromtimestamp(validated.created_at / 1000)
        
        updated_at = None
        if validated.updated_at:
            updated_at = datetime.fromtimestamp(validated.updated_at / 1000)
        
        # Map steps
        steps = []
        for step_data in validated.steps:
            # Map step fields
            step = ProtocolStep(
                id=step_data.get("id", str(uuid.uuid4())),
                name=step_data.get("name", ""),
                description=step_data.get("description", ""),
                step_type=step_data.get("type", "manual"),
                parameters=step_data.get("parameters", {}),
            )
            
            # Map duration
            duration_value = step_data.get("duration")
            duration_unit = step_data.get("durationUnit", "min")
            if duration_value is not None:
                step.duration = {
                    "value": duration_value,
                    "unit": duration_unit
                }
            
            # Map temperature
            temperature_value = step_data.get("temperature")
            temperature_unit = step_data.get("temperatureUnit", "°C")
            if temperature_value is not None:
                step.temperature = {
                    "value": temperature_value,
                    "unit": temperature_unit,
                    "tolerance": step_data.get("temperatureTolerance", 2)
                }
            
            steps.append(step)
        
        # Create Protocol instance
        protocol = Protocol(
            id=validated.id,
            name=validated.name,
            description=validated.description or "",
            steps=steps,
            version=validated.version or 1,
            parent_id=validated.parent_id,
            template_type=template_type,
            author=validated.created_by,
            parameters=validated.parameters or {},
            created_at=created_at,
            updated_at=updated_at,
        )
        
        return protocol
    
    def _convert_to_convex(self, protocol: Protocol) -> Dict[str, Any]:
        """
        Convert a protocol from experimental data format to Convex format.
        
        Args:
            protocol: Protocol in experimental data format
            
        Returns:
            Protocol in Convex format
        """
        # Convert steps
        steps = []
        for step in protocol.steps:
            step_data = {
                "id": step.id,
                "name": step.name,
                "description": step.description,
                "type": step.step_type,
                "parameters": step.parameters,
            }
            
            # Add duration if present
            if step.duration:
                step_data["duration"] = step.duration.get("value")
                step_data["durationUnit"] = step.duration.get("unit", "min")
            
            # Add temperature if present
            if step.temperature:
                step_data["temperature"] = step.temperature.get("value")
                step_data["temperatureUnit"] = step.temperature.get("unit", "°C")
                step_data["temperatureTolerance"] = step.temperature.get("tolerance", 2)
            
            steps.append(step_data)
        
        # Create Convex protocol object
        convex_protocol = {
            "name": protocol.name,
            "description": protocol.description,
            "steps": steps,
            "version": protocol.version,
            "parentId": protocol.parent_id,
            "category": protocol.template_type,
            "parameters": protocol.parameters,
            "isTemplate": protocol.template_type != "standard"
        }
        
        # Add ID if present
        if protocol.id:
            convex_protocol["id"] = protocol.id
        
        return convex_protocol