#!/usr/bin/env python3
"""
Enhanced property verification that can perform comprehensive analysis 
of our mock_rdkit_formula calculations compared to reference values.

This script:
1. Tests against known reference values from literature
2. Can use RDKit if available, but works without it
3. Provides detailed accuracy analysis by property type
4. Generates recommendations for improvement
"""

import os
import sys
import json
import logging
import statistics
from collections import defaultdict

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Import our mock RDKit formula module
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import mock_rdkit_formula

# Try to import rdkit
try:
    from rdkit import Chem, __version__ as rdkit_version
    from rdkit.Chem import Descriptors, Lipinski, MolSurf, rdMolDescriptors
    RDKIT_AVAILABLE = True
    logger.info(f"Using RDKit version {rdkit_version} for calculations")
except ImportError:
    RDKIT_AVAILABLE = False
    logger.info("RDKit not available, calculations will use only reference values and mock implementation")

# Known reference values for common molecules
# These values are standardized reference values from literature
REFERENCE_MOLECULES = [
    {
        "name": "Glycerol",
        "smiles": "C(C(CO)O)O",
        "properties": {
            "molecular_formula": "C3H8O3",
            "molecular_weight": 92.09,
            "logp": -1.76,
            "tpsa": 60.69,
            "h_donors": 3,
            "h_acceptors": 3,
            "rotatable_bonds": 2,
            "ring_count": 0,
            "aromatic_ring_count": 0,
            "heavy_atom_count": 6
        }
    },
    {
        "name": "DMSO",
        "smiles": "CS(=O)C",
        "properties": {
            "molecular_formula": "C2H6OS",
            "molecular_weight": 78.13,
            "logp": -1.35,
            "tpsa": 36.28,
            "h_donors": 0,
            "h_acceptors": 1,
            "rotatable_bonds": 1,
            "ring_count": 0,
            "aromatic_ring_count": 0,
            "heavy_atom_count": 4
        }
    },
    {
        "name": "Ethylene Glycol",
        "smiles": "OCCO",
        "properties": {
            "molecular_formula": "C2H6O2",
            "molecular_weight": 62.07,
            "logp": -1.20,
            "tpsa": 40.46,
            "h_donors": 2,
            "h_acceptors": 2,
            "rotatable_bonds": 1,
            "ring_count": 0,
            "aromatic_ring_count": 0,
            "heavy_atom_count": 4
        }
    },
    {
        "name": "Sucrose",
        "smiles": "C(C1C(C(C(C(O1)OC2C(C(C(C(O2)CO)O)O)O)O)CO)O)O",
        "properties": {
            "molecular_formula": "C12H22O11",
            "molecular_weight": 342.30,
            "logp": -3.76,
            "tpsa": 189.53,
            "h_donors": 8,
            "h_acceptors": 11,
            "rotatable_bonds": 5,
            "ring_count": 2,
            "aromatic_ring_count": 0,
            "heavy_atom_count": 23
        }
    },
    {
        "name": "Trehalose",
        "smiles": "C(C1C(C(C(C(O1)OC2C(C(C(C(O2)O)O)O)O)O)O)O)O",
        "properties": {
            "molecular_formula": "C12H22O11",
            "molecular_weight": 342.30,
            "logp": -4.20,
            "tpsa": 189.53,
            "h_donors": 8,
            "h_acceptors": 11,
            "rotatable_bonds": 3,
            "ring_count": 2,
            "aromatic_ring_count": 0,
            "heavy_atom_count": 23
        }
    },
    {
        "name": "Mannitol",
        "smiles": "OCC(O)C(O)C(O)C(O)CO",
        "properties": {
            "molecular_formula": "C6H14O6",
            "molecular_weight": 182.17,
            "logp": -3.10,
            "tpsa": 121.38,
            "h_donors": 6,
            "h_acceptors": 6,
            "rotatable_bonds": 5,
            "ring_count": 0,
            "aromatic_ring_count": 0,
            "heavy_atom_count": 12
        }
    },
    {
        "name": "Proline",
        "smiles": "OC(=O)C1CCCN1",
        "properties": {
            "molecular_formula": "C5H9NO2",
            "molecular_weight": 115.13,
            "logp": -0.06,
            "tpsa": 49.33,
            "h_donors": 1,
            "h_acceptors": 3,
            "rotatable_bonds": 1,
            "ring_count": 1,
            "aromatic_ring_count": 0,
            "heavy_atom_count": 8
        }
    },
    {
        "name": "Benzene",
        "smiles": "c1ccccc1",
        "properties": {
            "molecular_formula": "C6H6",
            "molecular_weight": 78.11,
            "logp": 2.13,
            "tpsa": 0.00,
            "h_donors": 0,
            "h_acceptors": 0,
            "rotatable_bonds": 0,
            "ring_count": 1,
            "aromatic_ring_count": 1,
            "heavy_atom_count": 6
        }
    },
    {
        "name": "Aspirin",
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "properties": {
            "molecular_formula": "C9H8O4",
            "molecular_weight": 180.16,
            "logp": 1.43,
            "tpsa": 63.60,
            "h_donors": 1,
            "h_acceptors": 4,
            "rotatable_bonds": 3,
            "ring_count": 1,
            "aromatic_ring_count": 1,
            "heavy_atom_count": 13
        }
    },
    {
        "name": "Caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "properties": {
            "molecular_formula": "C8H10N4O2",
            "molecular_weight": 194.19,
            "logp": -0.07,
            "tpsa": 58.44,
            "h_donors": 0,
            "h_acceptors": 6,
            "rotatable_bonds": 0,
            "ring_count": 2,
            "aromatic_ring_count": 1,
            "heavy_atom_count": 14
        }
    }
]

def calculate_mock_properties(smiles):
    """Calculate properties using our mock_rdkit_formula functions."""
    try:
        properties = {
            "molecular_formula": mock_rdkit_formula.calculate_molecular_formula(smiles),
            "molecular_weight": mock_rdkit_formula.calculate_molecular_weight(smiles),
            "logp": mock_rdkit_formula.calculate_logp(smiles),
            "tpsa": mock_rdkit_formula.calculate_tpsa(smiles),
            "h_donors": mock_rdkit_formula.calculate_h_donors(smiles),
            "h_acceptors": mock_rdkit_formula.calculate_h_acceptors(smiles),
            "rotatable_bonds": mock_rdkit_formula.calculate_rotatable_bonds(smiles),
            "ring_count": mock_rdkit_formula.calculate_ring_count(smiles),
            "aromatic_ring_count": mock_rdkit_formula.calculate_aromatic_ring_count(smiles),
            "heavy_atom_count": mock_rdkit_formula.calculate_heavy_atom_count(smiles)
        }
        return properties
    except Exception as e:
        logger.error(f"Error calculating mock properties: {e}")
        return None

def calculate_rdkit_properties(smiles):
    """Calculate properties using RDKit (if available)."""
    if not RDKIT_AVAILABLE:
        return None
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            properties = {
                "molecular_formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
                "molecular_weight": round(Descriptors.MolWt(mol), 2),
                "logp": round(Descriptors.MolLogP(mol), 2),
                "tpsa": round(MolSurf.TPSA(mol), 2),
                "h_donors": Lipinski.NumHDonors(mol),
                "h_acceptors": Lipinski.NumHAcceptors(mol),
                "rotatable_bonds": Descriptors.NumRotatableBonds(mol),
                "ring_count": Chem.rdMolDescriptors.CalcNumRings(mol),
                "aromatic_ring_count": Chem.rdMolDescriptors.CalcNumAromaticRings(mol),
                "heavy_atom_count": mol.GetNumHeavyAtoms()
            }
            return properties
        else:
            logger.warning(f"Failed to parse SMILES: {smiles}")
            return None
    except Exception as e:
        logger.error(f"Error calculating RDKit properties: {e}")
        return None

def compare_properties(reference, calculated):
    """Compare reference values with calculated values and compute accuracy."""
    results = {}
    
    for prop in reference:
        if reference[prop] is not None and calculated[prop] is not None:
            ref_val = reference[prop]
            calc_val = calculated[prop]
            
            # For numeric values, calculate percent error
            if isinstance(ref_val, (int, float)) and isinstance(calc_val, (int, float)):
                if ref_val == 0:
                    # Avoid division by zero
                    accuracy = 100.0 if calc_val == 0 else 0.0
                else:
                    percent_error = abs((calc_val - ref_val) / ref_val) * 100
                    accuracy = max(0, 100 - min(percent_error, 100))
            
            # For string values (like formula), exact match only
            elif isinstance(ref_val, str) and isinstance(calc_val, str):
                accuracy = 100.0 if ref_val == calc_val else 0.0
            
            # For other types, just check equality
            else:
                accuracy = 100.0 if ref_val == calc_val else 0.0
            
            results[prop] = {
                "reference": ref_val,
                "calculated": calc_val,
                "accuracy": accuracy
            }
    
    # Calculate overall accuracy
    if results:
        overall_accuracy = sum(item["accuracy"] for item in results.values()) / len(results)
    else:
        overall_accuracy = 0.0
    
    return {
        "properties": results,
        "overall_accuracy": overall_accuracy
    }

def analyze_accuracy_by_property(all_results):
    """Analyze accuracy statistics for each property type."""
    property_stats = defaultdict(list)
    
    # Collect accuracy values by property type
    for result in all_results:
        if "accuracy" in result and "properties" in result["accuracy"]:
            for prop, details in result["accuracy"]["properties"].items():
                property_stats[prop].append(details["accuracy"])
    
    # Calculate statistics for each property
    stats = {}
    for prop, accuracies in property_stats.items():
        if accuracies:
            stats[prop] = {
                "mean": statistics.mean(accuracies),
                "median": statistics.median(accuracies),
                "min": min(accuracies),
                "max": max(accuracies),
                "std_dev": statistics.stdev(accuracies) if len(accuracies) > 1 else 0,
                "sample_size": len(accuracies)
            }
    
    return stats

def generate_recommendations(accuracy_stats):
    """Generate recommendations for improving the mock implementation."""
    recommendations = []
    
    # Identify low-accuracy properties
    low_accuracy_props = [
        (prop, stats["mean"]) 
        for prop, stats in accuracy_stats.items() 
        if stats["mean"] < 70
    ]
    
    if low_accuracy_props:
        low_accuracy_props.sort(key=lambda x: x[1])  # Sort by accuracy (lowest first)
        
        for prop, accuracy in low_accuracy_props:
            if prop == "molecular_formula":
                recommendations.append({
                    "property": prop,
                    "accuracy": accuracy,
                    "recommendation": (
                        "Improve molecular formula calculation by implementing better handling "
                        "of implicit hydrogens and explicit valence states. Consider using lookup "
                        "tables for common fragments."
                    )
                })
            elif prop == "logp":
                recommendations.append({
                    "property": prop,
                    "accuracy": accuracy,
                    "recommendation": (
                        "Improve LogP calculation by implementing fragment-based contribution method "
                        "and better handling of ring systems. Consider implementing Wildman-Crippen "
                        "or similar algorithm."
                    )
                })
            elif prop == "tpsa":
                recommendations.append({
                    "property": prop,
                    "accuracy": accuracy,
                    "recommendation": (
                        "Enhance TPSA calculation by using more accurate contribution values for "
                        "functional groups. Consider implementing the PSA contribution method by Ertl et al."
                    )
                })
            elif prop == "rotatable_bonds":
                recommendations.append({
                    "property": prop,
                    "accuracy": accuracy,
                    "recommendation": (
                        "Improve rotatable bond counting by better identifying single bonds that are "
                        "not in rings and are not terminal."
                    )
                })
            else:
                recommendations.append({
                    "property": prop,
                    "accuracy": accuracy,
                    "recommendation": f"Consider revising the calculation algorithm for {prop}."
                })
    
    # General recommendations
    if RDKIT_AVAILABLE:
        recommendations.append({
            "property": "general",
            "accuracy": None,
            "recommendation": (
                "When RDKit is available, use it for the most accurate calculations. "
                "The mock implementation should be used only as a fallback."
            )
        })
    else:
        recommendations.append({
            "property": "general",
            "accuracy": None,
            "recommendation": (
                "Continue working to resolve RDKit container issues, as RDKit will provide "
                "more accurate calculations than any mock implementation."
            )
        })
    
    # If there are high-accuracy properties, note them
    high_accuracy_props = [
        prop for prop, stats in accuracy_stats.items() if stats["mean"] >= 90
    ]
    
    if high_accuracy_props:
        recommendations.append({
            "property": "high_accuracy",
            "accuracy": None,
            "recommendation": (
                f"The current implementation performs well for these properties: "
                f"{', '.join(high_accuracy_props)}. No immediate improvements needed."
            )
        })
    
    return recommendations

def main():
    """Main function to verify property calculations."""
    import argparse
    parser = argparse.ArgumentParser(description="Enhanced property calculation verification")
    parser.add_argument("--output", type=str, default="property_verification_report.json", 
                        help="Output file for detailed report")
    args = parser.parse_args()
    
    # Process reference molecules
    all_results = []
    
    print("\nEnhanced Property Calculation Verification")
    print("=========================================")
    
    for molecule in REFERENCE_MOLECULES:
        name = molecule["name"]
        smiles = molecule["smiles"]
        reference_props = molecule["properties"]
        
        # Calculate with our mock implementation
        mock_props = calculate_mock_properties(smiles)
        
        # Calculate with RDKit if available
        rdkit_props = calculate_rdkit_properties(smiles) if RDKIT_AVAILABLE else None
        
        # Compare with reference values
        if mock_props:
            accuracy = compare_properties(reference_props, mock_props)
            
            result = {
                "name": name,
                "smiles": smiles,
                "reference": reference_props,
                "mock": mock_props,
                "accuracy": accuracy
            }
            
            if rdkit_props:
                rdkit_accuracy = compare_properties(reference_props, rdkit_props)
                result["rdkit"] = rdkit_props
                result["rdkit_accuracy"] = rdkit_accuracy
            
            all_results.append(result)
            
            # Print summary for this molecule
            print(f"\nMolecule: {name}")
            print(f"SMILES: {smiles}")
            print(f"Mock Implementation Accuracy: {accuracy['overall_accuracy']:.2f}%")
            
            if rdkit_props:
                print(f"RDKit Accuracy: {rdkit_accuracy['overall_accuracy']:.2f}%")
    
    # Analyze accuracy by property type
    accuracy_stats = analyze_accuracy_by_property(all_results)
    
    # Generate recommendations
    recommendations = generate_recommendations(accuracy_stats)
    
    # Print summary statistics
    print("\nAccuracy by Property Type")
    print("========================")
    for prop, stats in sorted(accuracy_stats.items(), key=lambda x: x[1]["mean"], reverse=True):
        print(f"{prop}: {stats['mean']:.2f}% (min: {stats['min']:.2f}%, max: {stats['max']:.2f}%)")
    
    # Print overall accuracy
    overall_accuracies = [r["accuracy"]["overall_accuracy"] for r in all_results]
    avg_accuracy = statistics.mean(overall_accuracies) if overall_accuracies else 0
    
    print(f"\nOverall Average Accuracy: {avg_accuracy:.2f}%")
    
    # Print recommendations
    print("\nRecommendations for Improvement")
    print("==============================")
    for rec in recommendations:
        if rec["accuracy"] is not None:
            print(f"{rec['property']} ({rec['accuracy']:.2f}% accuracy): {rec['recommendation']}")
        else:
            print(f"{rec['property']}: {rec['recommendation']}")
    
    # Save detailed report
    report = {
        "summary": {
            "sample_size": len(all_results),
            "overall_accuracy": avg_accuracy,
            "property_accuracy": accuracy_stats
        },
        "recommendations": recommendations,
        "detailed_results": all_results
    }
    
    with open(args.output, 'w') as f:
        json.dump(report, f, indent=2)
    
    print(f"\nDetailed report saved to {args.output}")

if __name__ == "__main__":
    main()