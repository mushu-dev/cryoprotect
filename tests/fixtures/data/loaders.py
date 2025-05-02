"""
Test data loaders for CryoProtect v2.

This module provides functions for loading standard test data sets
from JSON files.
"""

import os
import json
from typing import Dict, List, Any, Optional

def get_data_dir() -> str:
    """
    Get the path to the data directory.
    
    Returns:
        Path to the data directory
    """
    return os.path.dirname(os.path.abspath(__file__))

def load_data_file(filename: str) -> List[Dict[str, Any]]:
    """
    Load data from a JSON file in the data directory.
    
    Args:
        filename: Name of the JSON file
        
    Returns:
        List of dictionaries with data
    """
    file_path = os.path.join(get_data_dir(), filename)
    
    with open(file_path, 'r') as f:
        return json.load(f)

def load_molecules() -> List[Dict[str, Any]]:
    """
    Load standard test molecules.
    
    Returns:
        List of molecule dictionaries
    """
    return load_data_file('molecules.json')

def load_property_types() -> List[Dict[str, Any]]:
    """
    Load standard test property types.
    
    Returns:
        List of property type dictionaries
    """
    return load_data_file('property_types.json')

def load_mixtures() -> List[Dict[str, Any]]:
    """
    Load standard test mixtures.
    
    Returns:
        List of mixture dictionaries
    """
    return load_data_file('mixtures.json')

def load_experiments() -> List[Dict[str, Any]]:
    """
    Load standard test experiments.
    
    Returns:
        List of experiment dictionaries
    """
    return load_data_file('experiments.json')

def load_calculation_methods() -> List[Dict[str, Any]]:
    """
    Load standard test calculation methods.
    
    Returns:
        List of calculation method dictionaries
    """
    return load_data_file('calculation_methods.json')

def load_teams() -> List[Dict[str, Any]]:
    """
    Load standard test teams.
    
    Returns:
        List of team dictionaries
    """
    return load_data_file('teams.json')

def load_projects() -> List[Dict[str, Any]]:
    """
    Load standard test projects.
    
    Returns:
        List of project dictionaries
    """
    return load_data_file('projects.json')

def load_related_data(entity_type: str, entity_id: str) -> Dict[str, Any]:
    """
    Load related data for a specific entity.
    
    Args:
        entity_type: Type of entity (e.g., 'molecule', 'mixture')
        entity_id: ID of the entity
        
    Returns:
        Dictionary with related data
    """
    if entity_type == 'molecule':
        # Load molecular properties for the molecule
        properties = []
        try:
            all_properties = load_data_file('molecular_properties.json')
            properties = [p for p in all_properties if p.get('molecule_id') == entity_id]
        except FileNotFoundError:
            pass
        
        # Load predictions for the molecule
        predictions = []
        try:
            all_predictions = load_data_file('predictions.json')
            predictions = [p for p in all_predictions if p.get('molecule_id') == entity_id]
        except FileNotFoundError:
            pass
        
        return {
            'properties': properties,
            'predictions': predictions
        }
    
    elif entity_type == 'mixture':
        # Load mixture components
        components = []
        try:
            all_components = load_data_file('mixture_components.json')
            components = [c for c in all_components if c.get('mixture_id') == entity_id]
            
            # Add molecule data to components
            molecules = load_molecules()
            molecules_by_id = {m['id']: m for m in molecules}
            
            for component in components:
                molecule_id = component.get('molecule_id')
                if molecule_id in molecules_by_id:
                    component['molecule'] = molecules_by_id[molecule_id]
        except FileNotFoundError:
            pass
        
        # Load experiments for the mixture
        experiments = []
        try:
            all_experiments = load_data_file('experiments.json')
            experiments = [e for e in all_experiments if e.get('mixture_id') == entity_id]
            
            # Add experiment properties to experiments
            try:
                all_properties = load_data_file('experiment_properties.json')
                for experiment in experiments:
                    experiment_id = experiment.get('id')
                    properties = [p for p in all_properties if p.get('experiment_id') == experiment_id]
                    experiment['properties'] = properties
            except FileNotFoundError:
                pass
        except FileNotFoundError:
            pass
        
        # Load predictions for the mixture
        predictions = []
        try:
            all_predictions = load_data_file('predictions.json')
            predictions = [p for p in all_predictions if p.get('mixture_id') == entity_id]
        except FileNotFoundError:
            pass
        
        return {
            'components': components,
            'experiments': experiments,
            'predictions': predictions
        }
    
    elif entity_type == 'experiment':
        # Load experiment properties
        properties = []
        try:
            all_properties = load_data_file('experiment_properties.json')
            properties = [p for p in all_properties if p.get('experiment_id') == entity_id]
        except FileNotFoundError:
            pass
        
        return {
            'properties': properties
        }
    
    elif entity_type == 'project':
        # Load team data
        team = None
        try:
            project = None
            all_projects = load_data_file('projects.json')
            for p in all_projects:
                if p.get('id') == entity_id:
                    project = p
                    break
            
            if project and 'team_id' in project:
                all_teams = load_data_file('teams.json')
                for t in all_teams:
                    if t.get('id') == project['team_id']:
                        team = t
                        break
        except FileNotFoundError:
            pass
        
        return {
            'team': team
        }
    
    # Default: return empty dict for unknown entity types
    return {}