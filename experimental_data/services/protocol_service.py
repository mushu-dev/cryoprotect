#!/usr/bin/env python3
"""
Protocol Service for CryoProtect Enhanced Experimental Data System.

This module provides services for protocol management, including creating,
retrieving, updating, and versioning experimental protocols.
"""

from typing import Dict, List, Any, Optional, Union, Tuple
import uuid
from datetime import datetime
import json
import difflib

from ..models import (
    Protocol,
    ProtocolStep,
    ValidationError,
    Provenance
)

class ProtocolService:
    """Service for managing experimental protocols."""
    
    def __init__(self, db_adapter):
        """
        Initialize protocol service.
        
        Args:
            db_adapter: Database adapter for persistence
        """
        self.db_adapter = db_adapter
    
    async def create_protocol(self, protocol_data: Dict[str, Any]) -> Protocol:
        """
        Create a new protocol.
        
        Args:
            protocol_data: Data for the new protocol
            
        Returns:
            Created protocol
            
        Raises:
            ValidationError: If protocol data is invalid
        """
        # Create protocol model
        protocol = Protocol.from_dict(protocol_data)
        
        # Validate protocol
        protocol.validate()
        
        # Persist to database
        protocol_dict = protocol.to_dict()
        protocol_id = await self.db_adapter.create('protocols', protocol_dict)
        protocol.id = protocol_id
        
        return protocol
    
    async def get_protocol(self, protocol_id: str) -> Optional[Protocol]:
        """
        Get a protocol by ID.
        
        Args:
            protocol_id: ID of the protocol
            
        Returns:
            Protocol if found, None otherwise
        """
        protocol_data = await self.db_adapter.get('protocols', protocol_id)
        
        if not protocol_data:
            return None
        
        return Protocol.from_dict(protocol_data)
    
    async def list_protocols(
        self, 
        filters: Optional[Dict[str, Any]] = None,
        page: int = 1,
        page_size: int = 20,
        sort_by: str = 'created_at',
        sort_order: str = 'desc'
    ) -> Tuple[List[Protocol], int]:
        """
        List protocols with pagination and filtering.
        
        Args:
            filters: Filters to apply
            page: Page number (1-indexed)
            page_size: Page size
            sort_by: Field to sort by
            sort_order: Sort order ('asc' or 'desc')
            
        Returns:
            Tuple of (protocols, total_count)
        """
        # Apply pagination
        offset = (page - 1) * page_size
        limit = page_size
        
        # Fetch protocols
        protocols_data, total_count = await self.db_adapter.list(
            'protocols',
            filters=filters,
            offset=offset,
            limit=limit,
            sort_by=sort_by,
            sort_order=sort_order
        )
        
        # Convert to models
        protocols = [Protocol.from_dict(data) for data in protocols_data]
        
        return protocols, total_count
    
    async def update_protocol(self, protocol_id: str, data: Dict[str, Any]) -> Optional[Protocol]:
        """
        Update a protocol.
        
        Args:
            protocol_id: ID of the protocol to update
            data: Fields to update
            
        Returns:
            Updated protocol if found, None otherwise
            
        Raises:
            ValidationError: If updated data is invalid
        """
        # Get protocol
        protocol = await self.get_protocol(protocol_id)
        
        if not protocol:
            return None
        
        # Update fields
        protocol.update(data)
        
        # Persist to database
        protocol_dict = protocol.to_dict()
        await self.db_adapter.update('protocols', protocol_id, protocol_dict)
        
        return protocol
    
    async def delete_protocol(self, protocol_id: str) -> bool:
        """
        Delete a protocol.
        
        Args:
            protocol_id: ID of the protocol to delete
            
        Returns:
            True if deleted, False if not found
        """
        # Check if protocol exists
        protocol = await self.get_protocol(protocol_id)
        
        if not protocol:
            return False
        
        # Check if protocol is used by any experiments
        experiments_using_protocol = await self.db_adapter.count(
            'experiments',
            filters={'protocol_id': protocol_id}
        )
        
        if experiments_using_protocol > 0:
            raise ValidationError(
                f"Cannot delete protocol: it is used by {experiments_using_protocol} experiments"
            )
        
        # Delete protocol
        result = await self.db_adapter.delete('protocols', protocol_id)
        
        return result
    
    async def create_protocol_version(self, protocol_id: str, changes: Dict[str, Any] = None) -> Protocol:
        """
        Create a new version of a protocol.
        
        Args:
            protocol_id: ID of the protocol to version
            changes: Changes to apply to the new version
            
        Returns:
            New protocol version
            
        Raises:
            ValidationError: If protocol not found or changes are invalid
        """
        # Get protocol
        protocol = await self.get_protocol(protocol_id)
        
        if not protocol:
            raise ValidationError(f"Protocol not found: {protocol_id}")
        
        # Create new version
        new_protocol = protocol.create_version()
        
        # Apply changes if provided
        if changes:
            new_protocol.update(changes)
        
        # Persist to database
        new_protocol_dict = new_protocol.to_dict()
        new_protocol_id = await self.db_adapter.create('protocols', new_protocol_dict)
        new_protocol.id = new_protocol_id
        
        return new_protocol
    
    async def get_protocol_version_history(self, protocol_id: str) -> List[Dict[str, Any]]:
        """
        Get version history for a protocol.
        
        Args:
            protocol_id: ID of the protocol
            
        Returns:
            List of version history entries
        """
        # Get protocol
        protocol = await self.get_protocol(protocol_id)
        
        if not protocol:
            return []
        
        # Initialize history with this protocol
        history = [{
            'id': protocol.id,
            'version': protocol.version,
            'created_at': protocol.created_at.isoformat() if protocol.created_at else None,
            'created_by': protocol.created_by,
            'name': protocol.name
        }]
        
        # Find ancestors
        current_parent_id = protocol.parent_id
        while current_parent_id:
            parent = await self.get_protocol(current_parent_id)
            if not parent:
                break
            
            history.append({
                'id': parent.id,
                'version': parent.version,
                'created_at': parent.created_at.isoformat() if parent.created_at else None,
                'created_by': parent.created_by,
                'name': parent.name
            })
            
            current_parent_id = parent.parent_id
        
        # Reverse to get chronological order
        history.reverse()
        
        return history
    
    async def compare_protocols(
        self, 
        protocol_id1: str, 
        protocol_id2: str
    ) -> Dict[str, Any]:
        """
        Compare two protocols and identify differences.
        
        Args:
            protocol_id1: ID of the first protocol
            protocol_id2: ID of the second protocol
            
        Returns:
            Comparison results
            
        Raises:
            ValidationError: If either protocol not found
        """
        # Get protocols
        protocol1 = await self.get_protocol(protocol_id1)
        protocol2 = await self.get_protocol(protocol_id2)
        
        if not protocol1:
            raise ValidationError(f"Protocol not found: {protocol_id1}")
        
        if not protocol2:
            raise ValidationError(f"Protocol not found: {protocol_id2}")
        
        # Compare basic metadata
        metadata_diff = {
            'name': protocol1.name != protocol2.name,
            'description': protocol1.description != protocol2.description,
            'template_type': protocol1.template_type != protocol2.template_type,
            'target_application': protocol1.target_application != protocol2.target_application,
            'parameters': protocol1.parameters != protocol2.parameters
        }
        
        # Compare steps
        steps1 = {step.id: step.to_dict() for step in protocol1.steps}
        steps2 = {step.id: step.to_dict() for step in protocol2.steps}
        
        # Find common, added, and removed steps
        step_ids1 = set(steps1.keys())
        step_ids2 = set(steps2.keys())
        
        common_step_ids = step_ids1.intersection(step_ids2)
        added_step_ids = step_ids2 - step_ids1
        removed_step_ids = step_ids1 - step_ids2
        
        # Compare common steps
        step_differences = {}
        for step_id in common_step_ids:
            step1 = steps1[step_id]
            step2 = steps2[step_id]
            
            if step1 != step2:
                # Find field differences
                field_diffs = {}
                for key in set(step1.keys()).union(step2.keys()):
                    if key not in step1:
                        field_diffs[key] = {'status': 'added', 'value': step2[key]}
                    elif key not in step2:
                        field_diffs[key] = {'status': 'removed', 'value': step1[key]}
                    elif step1[key] != step2[key]:
                        field_diffs[key] = {
                            'status': 'modified',
                            'old_value': step1[key],
                            'new_value': step2[key]
                        }
                
                step_differences[step_id] = field_diffs
        
        # Compile comparison results
        comparison = {
            'protocol1': {
                'id': protocol1.id,
                'version': protocol1.version,
                'name': protocol1.name,
                'created_at': protocol1.created_at.isoformat() if protocol1.created_at else None
            },
            'protocol2': {
                'id': protocol2.id,
                'version': protocol2.version,
                'name': protocol2.name,
                'created_at': protocol2.created_at.isoformat() if protocol2.created_at else None
            },
            'metadata_diff': metadata_diff,
            'steps': {
                'common': len(common_step_ids),
                'added': len(added_step_ids),
                'removed': len(removed_step_ids),
                'modified': len(step_differences)
            },
            'added_steps': [steps2[step_id] for step_id in added_step_ids],
            'removed_steps': [steps1[step_id] for step_id in removed_step_ids],
            'modified_steps': step_differences
        }
        
        return comparison
    
    async def export_protocol(self, protocol_id: str, format: str = 'json') -> str:
        """
        Export a protocol to a portable format.
        
        Args:
            protocol_id: ID of the protocol
            format: Export format ('json', 'yaml', or 'human')
            
        Returns:
            Protocol in the requested format
            
        Raises:
            ValidationError: If protocol not found or format invalid
        """
        # Get protocol
        protocol = await self.get_protocol(protocol_id)
        
        if not protocol:
            raise ValidationError(f"Protocol not found: {protocol_id}")
        
        # Export in requested format
        if format == 'json':
            return json.dumps(protocol.to_dict(), indent=2)
        elif format == 'yaml':
            import yaml
            return yaml.dump(protocol.to_dict())
        elif format == 'human':
            return self._format_protocol_for_humans(protocol)
        else:
            raise ValidationError(f"Unsupported export format: {format}")
    
    def _format_protocol_for_humans(self, protocol: Protocol) -> str:
        """
        Format protocol for human readability.
        
        Args:
            protocol: Protocol to format
            
        Returns:
            Human-readable protocol text
        """
        lines = [
            f"Protocol: {protocol.name} (Version {protocol.version})",
            f"Description: {protocol.description}",
            ""
        ]
        
        if protocol.target_application:
            lines.append(f"Target Application: {protocol.target_application}")
        
        if protocol.author:
            lines.append(f"Author: {protocol.author}")
        
        lines.append("")
        lines.append("Steps:")
        lines.append("")
        
        for i, step in enumerate(protocol.steps, 1):
            lines.append(f"Step {i}: {step.name}")
            lines.append("-" * (len(f"Step {i}: {step.name}")))
            lines.append(f"Description: {step.description}")
            
            if step.duration:
                lines.append(f"Duration: {step.duration['value']} {step.duration['unit']}")
            
            if step.temperature:
                lines.append(f"Temperature: {step.temperature['value']} {step.temperature['unit']} ± {step.temperature.get('tolerance', 0)} {step.temperature['unit']}")
            
            if step.equipment:
                lines.append("Equipment:")
                for eq in step.equipment:
                    lines.append(f"  - {eq.get('name', 'Unknown equipment')}")
            
            if step.parameters:
                lines.append("Parameters:")
                for key, value in step.parameters.items():
                    lines.append(f"  - {key}: {value}")
            
            lines.append("")
        
        lines.append("Parameters:")
        for key, value in protocol.parameters.items():
            lines.append(f"  {key}: {value}")
        
        return "\n".join(lines)
    
    async def validate_protocol(self, protocol_id: str) -> Dict[str, Any]:
        """
        Validate a protocol against scientific and logical criteria.
        
        Args:
            protocol_id: ID of the protocol
            
        Returns:
            Validation results
            
        Raises:
            ValidationError: If protocol not found
        """
        # Get protocol
        protocol = await self.get_protocol(protocol_id)
        
        if not protocol:
            raise ValidationError(f"Protocol not found: {protocol_id}")
        
        # Perform validation
        validation_results = {
            'protocol_id': protocol_id,
            'protocol_name': protocol.name,
            'validation_timestamp': datetime.now().isoformat(),
            'is_valid': True,
            'issues': []
        }
        
        # Check step prerequisites
        step_ids = {step.id for step in protocol.steps}
        for step in protocol.steps:
            for prereq_id in step.prerequisites:
                if prereq_id not in step_ids:
                    validation_results['is_valid'] = False
                    validation_results['issues'].append({
                        'severity': 'error',
                        'type': 'missing_prerequisite',
                        'message': f"Step {step.id} references missing prerequisite {prereq_id}",
                        'step_id': step.id
                    })
        
        # Check for empty steps
        for step in protocol.steps:
            if not step.name:
                validation_results['is_valid'] = False
                validation_results['issues'].append({
                    'severity': 'error',
                    'type': 'empty_step_name',
                    'message': "Step name cannot be empty",
                    'step_id': step.id
                })
            
            if not step.description:
                validation_results['issues'].append({
                    'severity': 'warning',
                    'type': 'empty_step_description',
                    'message': "Step description is empty",
                    'step_id': step.id
                })
        
        # Check for duplicate step names
        step_names = [step.name for step in protocol.steps]
        for name in set(step_names):
            if step_names.count(name) > 1:
                validation_results['issues'].append({
                    'severity': 'warning',
                    'type': 'duplicate_step_name',
                    'message': f"Duplicate step name: {name}",
                    'value': name
                })
        
        # Check for extreme temperatures
        for step in protocol.steps:
            if step.temperature and 'value' in step.temperature:
                temp_value = step.temperature['value']
                temp_unit = step.temperature.get('unit', '°C')
                
                if temp_unit == '°C':
                    if temp_value < -200:
                        validation_results['issues'].append({
                            'severity': 'warning',
                            'type': 'extreme_temperature',
                            'message': f"Temperature is extremely low: {temp_value}°C",
                            'step_id': step.id
                        })
                    elif temp_value > 100:
                        validation_results['issues'].append({
                            'severity': 'warning',
                            'type': 'extreme_temperature',
                            'message': f"Temperature is extremely high: {temp_value}°C",
                            'step_id': step.id
                        })
                elif temp_unit == 'K':
                    if temp_value < 70:
                        validation_results['issues'].append({
                            'severity': 'warning',
                            'type': 'extreme_temperature',
                            'message': f"Temperature is extremely low: {temp_value}K",
                            'step_id': step.id
                        })
                    elif temp_value > 373:
                        validation_results['issues'].append({
                            'severity': 'warning',
                            'type': 'extreme_temperature',
                            'message': f"Temperature is extremely high: {temp_value}K",
                            'step_id': step.id
                        })
        
        # Check for unreasonable durations
        for step in protocol.steps:
            if step.duration and 'value' in step.duration:
                duration_value = step.duration['value']
                duration_unit = step.duration.get('unit', 'min')
                
                if duration_unit == 's' and duration_value > 86400:
                    validation_results['issues'].append({
                        'severity': 'warning',
                        'type': 'long_duration',
                        'message': f"Step duration is very long: {duration_value} seconds",
                        'step_id': step.id
                    })
                elif duration_unit == 'min' and duration_value > 1440:
                    validation_results['issues'].append({
                        'severity': 'warning',
                        'type': 'long_duration',
                        'message': f"Step duration is very long: {duration_value} minutes",
                        'step_id': step.id
                    })
                elif duration_unit == 'h' and duration_value > 72:
                    validation_results['issues'].append({
                        'severity': 'warning',
                        'type': 'long_duration',
                        'message': f"Step duration is very long: {duration_value} hours",
                        'step_id': step.id
                    })
        
        # Check for missing equipment
        for step in protocol.steps:
            if step.step_type in ['automated', 'measurement'] and (not step.equipment or len(step.equipment) == 0):
                validation_results['issues'].append({
                    'severity': 'warning',
                    'type': 'missing_equipment',
                    'message': f"Step type '{step.step_type}' should specify required equipment",
                    'step_id': step.id
                })
        
        return validation_results
    
    async def import_protocol(self, protocol_data: Dict[str, Any]) -> Protocol:
        """
        Import a protocol from external data.
        
        Args:
            protocol_data: Protocol data to import
            
        Returns:
            Imported protocol
            
        Raises:
            ValidationError: If protocol data is invalid
        """
        # Create protocol from data
        protocol = Protocol.from_dict(protocol_data)
        
        # Add import provenance
        protocol.add_provenance_step(
            step_name="import",
            parameters={"origin": "external_import", "timestamp": datetime.now().isoformat()}
        )
        
        # Validate protocol
        protocol.validate()
        
        # Generate new ID to avoid collisions
        protocol.id = str(uuid.uuid4())
        
        # Persist to database
        protocol_dict = protocol.to_dict()
        protocol_id = await self.db_adapter.create('protocols', protocol_dict)
        protocol.id = protocol_id
        
        return protocol