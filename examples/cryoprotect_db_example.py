#!/usr/bin/env python3
"""
CryoProtect Analyzer - Database Usage Example (Python)

This example demonstrates how to interact with the CryoProtect database
using the Supabase Python client.

Prerequisites:
- Python 3.6+ installed
- Supabase project with the CryoProtect schema applied
- supabase-py package installed (pip install supabase)
"""

import os
import json
from datetime import date
from supabase import create_client, Client

# Replace with your Supabase URL and anon key
supabase_url = "https://your-project-ref.supabase.co"
supabase_key = "your-anon-key"
supabase: Client = create_client(supabase_url, supabase_key)


def import_molecule_from_pubchem(cid: int) -> str:
    """
    Example 1: Import a molecule from PubChem
    
    Args:
        cid: PubChem Compound ID
        
    Returns:
        UUID of the imported molecule
    """
    print(f"Importing molecule with CID {cid} from PubChem...")
    
    # Get the current user ID
    user_id = supabase.auth.current_user.id if supabase.auth.current_user else None
    
    # Call the database function to import the molecule
    response = supabase.rpc(
        "import_molecule_from_pubchem",
        {"p_cid": cid, "p_user_id": user_id}
    ).execute()
    
    if response.error:
        print(f"Error importing molecule: {response.error}")
        raise Exception(response.error)
    
    molecule_id = response.data
    print(f"Successfully imported molecule with ID: {molecule_id}")
    return molecule_id


def add_molecular_properties(molecule_id: str, properties: dict) -> list:
    """
    Example 2: Add molecular properties
    
    Args:
        molecule_id: UUID of the molecule
        properties: Dictionary with property name-value pairs
        
    Returns:
        List of inserted property records
    """
    print(f"Adding properties for molecule {molecule_id}...")
    
    # Get property types
    response = supabase.table("property_types").select("id, name, data_type").execute()
    if response.error:
        print(f"Error fetching property types: {response.error}")
        raise Exception(response.error)
    
    property_types = response.data
    
    # Get the current user ID
    user_id = supabase.auth.current_user.id if supabase.auth.current_user else None
    
    # Prepare property inserts
    property_inserts = []
    
    for property_name, value in properties.items():
        property_type = next((pt for pt in property_types if pt["name"] == property_name), None)
        if not property_type:
            print(f"Property type '{property_name}' not found, skipping")
            continue
        
        property_insert = {
            "molecule_id": molecule_id,
            "property_type_id": property_type["id"],
            "created_by": user_id
        }
        
        # Set the appropriate value field based on data type
        if property_type["data_type"] == "numeric":
            property_insert["numeric_value"] = float(value) if value is not None else None
        elif property_type["data_type"] == "text":
            property_insert["text_value"] = str(value) if value is not None else None
        elif property_type["data_type"] == "boolean":
            property_insert["boolean_value"] = bool(value) if value is not None else None
        else:
            print(f"Unknown data type '{property_type['data_type']}' for property '{property_name}', skipping")
            continue
        
        property_inserts.append(property_insert)
    
    # Insert properties
    response = supabase.table("molecular_properties").upsert(
        property_inserts, 
        on_conflict="molecule_id,property_type_id"
    ).execute()
    
    if response.error:
        print(f"Error adding properties: {response.error}")
        raise Exception(response.error)
    
    print(f"Successfully added {len(property_inserts)} properties")
    return response.data


def create_mixture(name: str, description: str, components: list) -> str:
    """
    Example 3: Create a mixture
    
    Args:
        name: Name of the mixture
        description: Description of the mixture
        components: List of dicts with moleculeId, concentration, unit
        
    Returns:
        UUID of the created mixture
    """
    print(f"Creating mixture '{name}'...")
    
    # Get the current user ID
    user_id = supabase.auth.current_user.id if supabase.auth.current_user else None
    
    # Create mixture
    response = supabase.table("mixtures").insert({
        "name": name,
        "description": description,
        "created_by": user_id
    }).execute()
    
    if response.error:
        print(f"Error creating mixture: {response.error}")
        raise Exception(response.error)
    
    mixture_id = response.data[0]["id"]
    
    # Add components
    component_inserts = [{
        "mixture_id": mixture_id,
        "molecule_id": comp["moleculeId"],
        "concentration": comp["concentration"],
        "concentration_unit": comp["unit"],
        "created_by": user_id
    } for comp in components]
    
    response = supabase.table("mixture_components").insert(component_inserts).execute()
    
    if response.error:
        print(f"Error adding mixture components: {response.error}")
        raise Exception(response.error)
    
    print(f"Successfully created mixture with ID: {mixture_id}")
    return mixture_id


def add_prediction(mixture_id: str, property_name: str, value, confidence: float, method_name: str):
    """
    Example 4: Add a prediction for a mixture
    
    Args:
        mixture_id: UUID of the mixture
        property_name: Name of the property being predicted
        value: Predicted value
        confidence: Confidence level (0-1)
        method_name: Name of the calculation method
    
    Returns:
        The inserted prediction record
    """
    print(f"Adding prediction for mixture {mixture_id}, property '{property_name}'...")
    
    # Get property type ID
    response = supabase.table("property_types").select("id, data_type").eq("name", property_name).execute()
    if response.error or not response.data:
        print(f"Error fetching property type '{property_name}': {response.error}")
        raise Exception(f"Property type '{property_name}' not found")
    
    property_type = response.data[0]
    
    # Get calculation method ID
    response = supabase.table("calculation_methods").select("id").eq("name", method_name).execute()
    if response.error or not response.data:
        print(f"Error fetching calculation method '{method_name}': {response.error}")
        raise Exception(f"Calculation method '{method_name}' not found")
    
    calculation_method = response.data[0]
    
    # Get the current user ID
    user_id = supabase.auth.current_user.id if supabase.auth.current_user else None
    
    # Prepare prediction insert
    prediction = {
        "mixture_id": mixture_id,
        "property_type_id": property_type["id"],
        "calculation_method_id": calculation_method["id"],
        "confidence": confidence,
        "created_by": user_id
    }
    
    # Set the appropriate value field based on data type
    if property_type["data_type"] == "numeric":
        prediction["numeric_value"] = float(value) if value is not None else None
    elif property_type["data_type"] == "text":
        prediction["text_value"] = str(value) if value is not None else None
    elif property_type["data_type"] == "boolean":
        prediction["boolean_value"] = bool(value) if value is not None else None
    else:
        print(f"Unknown data type '{property_type['data_type']}'")
        raise Exception(f"Unknown data type '{property_type['data_type']}'")
    
    # Insert prediction
    response = supabase.table("predictions").upsert(
        prediction,
        on_conflict="mixture_id,property_type_id,calculation_method_id"
    ).execute()
    
    if response.error:
        print(f"Error adding prediction: {response.error}")
        raise Exception(response.error)
    
    print("Successfully added prediction")
    return response.data


def record_experiment(mixture_id: str, property_name: str, value, conditions: str, experiment_date: str):
    """
    Example 5: Record experimental results
    
    Args:
        mixture_id: UUID of the mixture
        property_name: Name of the property being measured
        value: Measured value
        conditions: Experimental conditions
        experiment_date: Date of experiment (YYYY-MM-DD)
    
    Returns:
        The inserted experiment record
    """
    print(f"Recording experiment for mixture {mixture_id}, property '{property_name}'...")
    
    # Get property type ID
    response = supabase.table("property_types").select("id, data_type").eq("name", property_name).execute()
    if response.error or not response.data:
        print(f"Error fetching property type '{property_name}': {response.error}")
        raise Exception(f"Property type '{property_name}' not found")
    
    property_type = response.data[0]
    
    # Get the current user ID
    user_id = supabase.auth.current_user.id if supabase.auth.current_user else None
    
    # Prepare experiment insert
    experiment = {
        "mixture_id": mixture_id,
        "property_type_id": property_type["id"],
        "experimental_conditions": conditions,
        "date_performed": experiment_date,
        "created_by": user_id
    }
    
    # Set the appropriate value field based on data type
    if property_type["data_type"] == "numeric":
        experiment["numeric_value"] = float(value) if value is not None else None
    elif property_type["data_type"] == "text":
        experiment["text_value"] = str(value) if value is not None else None
    elif property_type["data_type"] == "boolean":
        experiment["boolean_value"] = bool(value) if value is not None else None
    else:
        print(f"Unknown data type '{property_type['data_type']}'")
        raise Exception(f"Unknown data type '{property_type['data_type']}'")
    
    # Insert experiment
    response = supabase.table("experiments").insert(experiment).execute()
    
    if response.error:
        print(f"Error recording experiment: {response.error}")
        raise Exception(response.error)
    
    print("Successfully recorded experiment")
    return response.data


def compare_prediction_with_experiment(mixture_id: str, property_name: str):
    """
    Example 6: Compare prediction with experiment
    
    Args:
        mixture_id: UUID of the mixture
        property_name: Name of the property
    
    Returns:
        Comparison result
    """
    print(f"Comparing prediction with experiment for mixture {mixture_id}, property '{property_name}'...")
    
    # Get property type ID
    response = supabase.table("property_types").select("id").eq("name", property_name).execute()
    if response.error or not response.data:
        print(f"Error fetching property type '{property_name}': {response.error}")
        raise Exception(f"Property type '{property_name}' not found")
    
    property_type = response.data[0]
    
    # Call the database function to compare
    response = supabase.rpc(
        "compare_prediction_with_experiment",
        {
            "p_mixture_id": mixture_id,
            "p_property_type_id": property_type["id"]
        }
    ).execute()
    
    if response.error:
        print(f"Error comparing prediction with experiment: {response.error}")
        raise Exception(response.error)
    
    print(f"Comparison result: {json.dumps(response.data, indent=2)}")
    return response.data


def get_molecule_with_properties(molecule_id: str):
    """
    Example 7: Get molecule with all its properties
    
    Args:
        molecule_id: UUID of the molecule
    
    Returns:
        Molecule with properties
    """
    print(f"Getting molecule {molecule_id} with properties...")
    
    response = supabase.table("molecule_with_properties").select("*").eq("id", molecule_id).execute()
    
    if response.error:
        print(f"Error getting molecule with properties: {response.error}")
        raise Exception(response.error)
    
    if not response.data:
        print(f"Molecule {molecule_id} not found")
        return None
    
    print(f"Molecule with properties: {json.dumps(response.data[0], indent=2)}")
    return response.data[0]


def get_mixture_with_components(mixture_id: str):
    """
    Example 8: Get mixture with all its components
    
    Args:
        mixture_id: UUID of the mixture
    
    Returns:
        Mixture with components
    """
    print(f"Getting mixture {mixture_id} with components...")
    
    response = supabase.table("mixture_with_components").select("*").eq("id", mixture_id).execute()
    
    if response.error:
        print(f"Error getting mixture with components: {response.error}")
        raise Exception(response.error)
    
    if not response.data:
        print(f"Mixture {mixture_id} not found")
        return None
    
    print(f"Mixture with components: {json.dumps(response.data[0], indent=2)}")
    return response.data[0]


def run_example():
    """Example usage"""
    try:
        # Sign in (replace with your auth method)
        response = supabase.auth.sign_in_with_password({
            "email": "user@example.com",
            "password": "password"
        })
        
        if response.error:
            print(f"Authentication error: {response.error}")
            return
        
        # Example 1: Import a molecule from PubChem (Glycerol, CID: 753)
        molecule_id = import_molecule_from_pubchem(753)
        
        # Example 2: Add molecular properties
        add_molecular_properties(molecule_id, {
            "Molecular Weight": 92.09,
            "LogP": -1.76,
            "TPSA": 60.69,
            "H-Bond Donors": 3,
            "H-Bond Acceptors": 3,
            "Toxicity": "Low toxicity",
            "Stability": "Stable under normal conditions",
            "Environmental Safety": "Biodegradable",
            "Total Score": 180
        })
        
        # Example 3: Create a mixture
        mixture_id = create_mixture(
            "Glycerol-Water Solution",
            "A 30% glycerol solution in water",
            [
                {"moleculeId": molecule_id, "concentration": 30, "unit": "%"},
                # Water would be another component, but we'd need to import it first
            ]
        )
        
        # Example 4: Add a prediction
        add_prediction(
            mixture_id,
            "Freezing Point",
            -15.3,
            0.9,
            "CryoProtect Scoring"
        )
        
        # Example 5: Record experimental results
        record_experiment(
            mixture_id,
            "Freezing Point",
            -14.8,
            "Standard pressure, cooling rate 1°C/min",
            "2025-04-14"
        )
        
        # Example 6: Compare prediction with experiment
        comparison = compare_prediction_with_experiment(
            mixture_id,
            "Freezing Point"
        )
        
        # Example 7: Get molecule with properties
        molecule_with_props = get_molecule_with_properties(molecule_id)
        
        # Example 8: Get mixture with components
        mixture_with_comps = get_mixture_with_components(mixture_id)
        
        print("Example completed successfully!")
    except Exception as e:
        print(f"Example failed: {e}")


if __name__ == "__main__":
    run_example()