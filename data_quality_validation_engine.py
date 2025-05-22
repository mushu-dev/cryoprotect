#!/usr/bin/env python3
"""
Data Quality Validation Engine for CryoProtect
Real-time data validation and automated quality scoring system

This module provides comprehensive data validation, quality scoring,
and automated monitoring for the molecular database.
"""

import json
import requests
import time
import logging
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any
from dataclasses import dataclass, asdict

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@dataclass
class QualityMetrics:
    """Data quality metrics for a molecule"""
    molecule_id: str
    overall_score: float
    completeness_score: float
    accuracy_score: float
    consistency_score: float
    validation_errors: List[str]
    validation_warnings: List[str]
    last_validated: str
    
class DataQualityValidator:
    """
    Comprehensive data quality validation engine
    """
    
    def __init__(self, convex_url: str = "https://upbeat-parrot-866.convex.cloud"):
        self.convex_url = convex_url
        self.validation_rules = self._initialize_validation_rules()
        self.session = requests.Session()
        self.session.headers.update({
            'Content-Type': 'application/json',
            'User-Agent': 'CryoProtect-DataQuality/1.0'
        })
        self.session.timeout = 30
    
    def call_convex_function(self, function_path: str, args: Dict[str, Any], is_mutation: bool = False) -> Any:
        """Call a Convex function via HTTP API (using the working method)."""
        endpoint = "mutation" if is_mutation else "query"
        url = f"{self.convex_url}/api/{endpoint}"
        
        payload = {
            "path": function_path,
            "args": args,
            "format": "json"
        }
        
        try:
            response = self.session.post(url, json=payload, timeout=30)
            
            if response.status_code == 200:
                result = response.json()
                
                # Handle Convex response format
                if isinstance(result, dict) and 'status' in result:
                    if result['status'] == 'success':
                        return result.get('value')
                    else:
                        raise Exception(f"Convex error: {result.get('error', 'Unknown error')}")
                else:
                    return result
            else:
                raise Exception(f"HTTP {response.status_code}: {response.text}")
                
        except requests.RequestException as e:
            raise Exception(f"Request failed: {e}")
        
    def _initialize_validation_rules(self) -> Dict[str, Any]:
        """Initialize comprehensive validation rules for molecular data"""
        return {
            "required_fields": [
                "name", "smiles", "molecularWeight", "logP", "tpsa",
                "hbondDonors", "hbondAcceptors", "rotatablebonds"
            ],
            "critical_fields": [
                "exactMass", "complexity", "heavyAtomCount",
                "morganFingerprint", "rdkitFingerprint"
            ],
            "molecular_weight_range": (50, 2000),  # Da
            "logp_range": (-5, 10),
            "tpsa_range": (0, 300),  # Ų
            "hbond_range": (0, 20),
            "rotatable_bonds_range": (0, 50),
            "smiles_patterns": {
                "valid_chars": "CNOPSFClBrI[]()=+-#@",
                "min_length": 3,
                "max_length": 500
            },
            "scoring_weights": {
                "completeness": 0.40,
                "accuracy": 0.35,
                "consistency": 0.25
            }
        }
    
    def validate_molecule(self, molecule: Dict[str, Any]) -> QualityMetrics:
        """
        Validate a single molecule and return quality metrics
        """
        errors = []
        warnings = []
        
        # 1. Completeness validation
        completeness_score = self._validate_completeness(molecule, errors, warnings)
        
        # 2. Accuracy validation
        accuracy_score = self._validate_accuracy(molecule, errors, warnings)
        
        # 3. Consistency validation
        consistency_score = self._validate_consistency(molecule, errors, warnings)
        
        # Calculate overall score
        weights = self.validation_rules["scoring_weights"]
        overall_score = (
            completeness_score * weights["completeness"] +
            accuracy_score * weights["accuracy"] +
            consistency_score * weights["consistency"]
        )
        
        return QualityMetrics(
            molecule_id=molecule.get("_id", "unknown"),
            overall_score=round(overall_score, 2),
            completeness_score=round(completeness_score, 2),
            accuracy_score=round(accuracy_score, 2),
            consistency_score=round(consistency_score, 2),
            validation_errors=errors,
            validation_warnings=warnings,
            last_validated=datetime.now().isoformat()
        )
    
    def _validate_completeness(self, molecule: Dict[str, Any], errors: List[str], warnings: List[str]) -> float:
        """Validate data completeness"""
        required_fields = self.validation_rules["required_fields"]
        critical_fields = self.validation_rules["critical_fields"]
        
        missing_required = []
        missing_critical = []
        
        for field in required_fields:
            if not molecule.get(field):
                missing_required.append(field)
                
        for field in critical_fields:
            if not molecule.get(field):
                missing_critical.append(field)
        
        if missing_required:
            errors.append(f"Missing required fields: {', '.join(missing_required)}")
            
        if missing_critical:
            warnings.append(f"Missing critical fields: {', '.join(missing_critical)}")
        
        total_fields = len(required_fields) + len(critical_fields)
        missing_total = len(missing_required) + len(missing_critical)
        
        # Weight required fields more heavily than critical fields
        required_weight = 0.7
        critical_weight = 0.3
        
        required_completeness = (len(required_fields) - len(missing_required)) / len(required_fields)
        critical_completeness = (len(critical_fields) - len(missing_critical)) / len(critical_fields)
        
        return (required_completeness * required_weight + critical_completeness * critical_weight) * 100
    
    def _validate_accuracy(self, molecule: Dict[str, Any], errors: List[str], warnings: List[str]) -> float:
        """Validate data accuracy"""
        accuracy_score = 100.0
        
        # Validate molecular weight
        mw = molecule.get("molecularWeight")
        if mw:
            mw_range = self.validation_rules["molecular_weight_range"]
            if not (mw_range[0] <= mw <= mw_range[1]):
                errors.append(f"Molecular weight {mw} outside valid range {mw_range}")
                accuracy_score -= 20
        
        # Validate LogP
        logp = molecule.get("logP")
        if logp:
            logp_range = self.validation_rules["logp_range"]
            if not (logp_range[0] <= logp <= logp_range[1]):
                warnings.append(f"LogP {logp} outside typical range {logp_range}")
                accuracy_score -= 10
        
        # Validate TPSA
        tpsa = molecule.get("tpsa")
        if tpsa:
            tpsa_range = self.validation_rules["tpsa_range"]
            if not (tpsa_range[0] <= tpsa <= tpsa_range[1]):
                warnings.append(f"TPSA {tpsa} outside typical range {tpsa_range}")
                accuracy_score -= 10
        
        # Validate H-bond donors/acceptors
        for field in ["hbondDonors", "hbondAcceptors"]:
            value = molecule.get(field)
            if value:
                hbond_range = self.validation_rules["hbond_range"]
                if not (hbond_range[0] <= value <= hbond_range[1]):
                    warnings.append(f"{field} {value} outside typical range {hbond_range}")
                    accuracy_score -= 5
        
        # Validate SMILES format
        smiles = molecule.get("smiles")
        if smiles:
            if not self._validate_smiles_format(smiles):
                errors.append(f"Invalid SMILES format: {smiles}")
                accuracy_score -= 25
        
        return max(0, accuracy_score)
    
    def _validate_consistency(self, molecule: Dict[str, Any], errors: List[str], warnings: List[str]) -> float:
        """Validate data consistency"""
        consistency_score = 100.0
        
        # Check if exactMass and molecularWeight are consistent
        exact_mass = molecule.get("exactMass")
        molecular_weight = molecule.get("molecularWeight")
        
        if exact_mass and molecular_weight:
            # Allow 5% difference between exact mass and molecular weight
            diff_percent = abs(exact_mass - molecular_weight) / molecular_weight * 100
            if diff_percent > 5:
                warnings.append(f"Large difference between exactMass ({exact_mass}) and molecularWeight ({molecular_weight})")
                consistency_score -= 15
        
        # Check if heavy atom count is reasonable for molecular weight
        heavy_atoms = molecule.get("heavyAtomCount")
        if heavy_atoms and molecular_weight:
            # Typical range: 7-15 Da per heavy atom
            expected_mw_min = heavy_atoms * 7
            expected_mw_max = heavy_atoms * 15
            if not (expected_mw_min <= molecular_weight <= expected_mw_max):
                warnings.append(f"Heavy atom count ({heavy_atoms}) inconsistent with molecular weight ({molecular_weight})")
                consistency_score -= 10
        
        # Check if fingerprints are present and have reasonable length
        morgan_fp = molecule.get("morganFingerprint")
        rdkit_fp = molecule.get("rdkitFingerprint")
        
        if morgan_fp and len(morgan_fp) < 512:  # Morgan fingerprints should be 512+ bits
            warnings.append(f"Morgan fingerprint length ({len(morgan_fp)}) shorter than expected")
            consistency_score -= 5
            
        if rdkit_fp and len(rdkit_fp) < 1024:  # RDKit fingerprints should be 1024+ bits
            warnings.append(f"RDKit fingerprint length ({len(rdkit_fp)}) shorter than expected")
            consistency_score -= 5
        
        return max(0, consistency_score)
    
    def _validate_smiles_format(self, smiles: str) -> bool:
        """Basic SMILES format validation"""
        rules = self.validation_rules["smiles_patterns"]
        
        # Check length
        if not (rules["min_length"] <= len(smiles) <= rules["max_length"]):
            return False
        
        # Check for balanced brackets
        bracket_count = smiles.count('[') - smiles.count(']')
        paren_count = smiles.count('(') - smiles.count(')')
        
        if bracket_count != 0 or paren_count != 0:
            return False
        
        return True
    
    def validate_all_molecules(self) -> Dict[str, Any]:
        """
        Validate all molecules in the database and return comprehensive report
        """
        logger.info("Starting comprehensive database validation...")
        
        try:
            # Get all molecules first
            molecules = self.call_convex_function("molecules:listMolecules", {"limit": 1000})
            
            # Then get properties for each molecule
            molecules_with_data = []
            for molecule in molecules:
                try:
                    mol_with_props = self.call_convex_function(
                        "molecules:getMoleculeWithProperties", 
                        {"id": molecule["_id"]}
                    )
                    
                    if mol_with_props:
                        # Flatten the structure for easier validation
                        flat_molecule = {
                            "_id": mol_with_props["_id"],
                            "name": mol_with_props["name"],
                            "pubchemCid": mol_with_props.get("pubchemCid"),
                            "smiles": mol_with_props.get("canonicalSmiles"),
                            "inchiKey": mol_with_props.get("inchiKey"),
                            "formula": mol_with_props.get("formula"),
                            "status": mol_with_props.get("status"),
                        }
                        
                        # Add molecular properties if they exist
                        if mol_with_props.get("properties"):
                            props = mol_with_props["properties"]
                            flat_molecule.update({
                                "molecularWeight": props.get("molecularWeight"),
                                "exactMass": props.get("exactMass"),
                                "logP": props.get("logP"),
                                "tpsa": props.get("tpsa"),
                                "hbondDonors": props.get("hbondDonors"),
                                "hbondAcceptors": props.get("hbondAcceptors"),
                                "rotatableBonds": props.get("rotatableBonds"),
                                "aromaticRings": props.get("aromaticRings"),
                                "aliphaticRings": props.get("aliphaticRings"),
                                "complexity": props.get("complexity"),
                                "heavyAtomCount": props.get("heavyAtomCount"),
                                "morganFingerprint": props.get("morganFingerprint"),
                                "rdkitFingerprint": props.get("rdkitFingerprint"),
                            })
                        
                        # Add cryoprotectant scores if they exist
                        if mol_with_props.get("scores"):
                            scores = mol_with_props["scores"]
                            flat_molecule.update({
                                "glassTempScore": scores.get("glassTempScore"),
                                "viscosityScore": scores.get("viscosityScore"),
                                "permeabilityScore": scores.get("permeabilityScore"),
                                "toxicityScore": scores.get("toxicityScore"),
                                "overallScore": scores.get("overallScore"),
                                "category": scores.get("category"),
                                "confidence": scores.get("confidence"),
                            })
                        
                        molecules_with_data.append(flat_molecule)
                        
                except Exception as e:
                    logger.warning(f"Error getting properties for {molecule.get('name', 'unknown')}: {e}")
                    # Add molecule without properties for completeness validation
                    molecules_with_data.append({
                        "_id": molecule["_id"],
                        "name": molecule["name"],
                        "pubchemCid": molecule.get("pubchemCid"),
                        "smiles": molecule.get("canonicalSmiles"),
                        "inchiKey": molecule.get("inchiKey"),
                        "formula": molecule.get("formula"),
                        "status": molecule.get("status"),
                    })
            
            molecules = molecules_with_data
            
            if not molecules:
                logger.warning("No molecules found in database")
                return {"error": "No molecules found", "validation_results": []}
            
            validation_results = []
            total_molecules = len(molecules)
            
            logger.info(f"Validating {total_molecules} molecules...")
            
            for i, molecule in enumerate(molecules, 1):
                logger.info(f"Validating molecule {i}/{total_molecules}: {molecule.get('name', 'unknown')}")
                
                metrics = self.validate_molecule(molecule)
                validation_results.append(asdict(metrics))
                
                # Progress indicator
                if i % 5 == 0:
                    logger.info(f"Progress: {i}/{total_molecules} molecules validated")
            
            # Generate summary statistics
            summary = self._generate_validation_summary(validation_results)
            
            report = {
                "validation_timestamp": datetime.now().isoformat(),
                "total_molecules": total_molecules,
                "summary": summary,
                "validation_results": validation_results
            }
            
            # Save report
            report_filename = f"data_quality_validation_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            with open(report_filename, 'w') as f:
                json.dump(report, f, indent=2)
            
            logger.info(f"Validation complete! Report saved to {report_filename}")
            return report
            
        except Exception as e:
            logger.error(f"Validation failed: {str(e)}")
            return {"error": str(e), "validation_results": []}
    
    def _generate_validation_summary(self, validation_results: List[Dict]) -> Dict[str, Any]:
        """Generate summary statistics from validation results"""
        if not validation_results:
            return {}
        
        scores = [result["overall_score"] for result in validation_results]
        completeness_scores = [result["completeness_score"] for result in validation_results]
        accuracy_scores = [result["accuracy_score"] for result in validation_results]
        consistency_scores = [result["consistency_score"] for result in validation_results]
        
        total_errors = sum(len(result["validation_errors"]) for result in validation_results)
        total_warnings = sum(len(result["validation_warnings"]) for result in validation_results)
        
        # Quality tiers
        excellent = len([s for s in scores if s >= 90])
        good = len([s for s in scores if 80 <= s < 90])
        fair = len([s for s in scores if 70 <= s < 80])
        poor = len([s for s in scores if s < 70])
        
        return {
            "overall_quality": {
                "average_score": round(sum(scores) / len(scores), 2),
                "median_score": round(sorted(scores)[len(scores) // 2], 2),
                "min_score": min(scores),
                "max_score": max(scores)
            },
            "score_breakdown": {
                "completeness_avg": round(sum(completeness_scores) / len(completeness_scores), 2),
                "accuracy_avg": round(sum(accuracy_scores) / len(accuracy_scores), 2),
                "consistency_avg": round(sum(consistency_scores) / len(consistency_scores), 2)
            },
            "quality_distribution": {
                "excellent_90_plus": excellent,
                "good_80_89": good,
                "fair_70_79": fair,
                "poor_below_70": poor
            },
            "issues": {
                "total_errors": total_errors,
                "total_warnings": total_warnings,
                "molecules_with_errors": len([r for r in validation_results if r["validation_errors"]]),
                "molecules_with_warnings": len([r for r in validation_results if r["validation_warnings"]])
            }
        }
    
    def get_quality_recommendations(self, validation_results: List[Dict]) -> List[str]:
        """Generate specific recommendations based on validation results"""
        recommendations = []
        
        # Analyze common issues
        all_errors = []
        all_warnings = []
        
        for result in validation_results:
            all_errors.extend(result["validation_errors"])
            all_warnings.extend(result["validation_warnings"])
        
        # Count error types
        missing_required = len([e for e in all_errors if "Missing required fields" in e])
        invalid_smiles = len([e for e in all_errors if "Invalid SMILES" in e])
        invalid_ranges = len([e for e in all_errors if "outside valid range" in e])
        
        if missing_required > 0:
            recommendations.append(f"🚨 HIGH PRIORITY: {missing_required} molecules missing required fields - implement data completion pipeline")
        
        if invalid_smiles > 0:
            recommendations.append(f"🔍 CRITICAL: {invalid_smiles} molecules have invalid SMILES - requires manual curation")
        
        if invalid_ranges > 0:
            recommendations.append(f"⚠️ MEDIUM: {invalid_ranges} molecules have values outside valid ranges - review calculations")
        
        # Performance recommendations
        avg_score = sum(r["overall_score"] for r in validation_results) / len(validation_results)
        
        if avg_score < 80:
            recommendations.append("📊 Implement automated data completion pipeline to improve overall quality scores")
        
        if avg_score >= 90:
            recommendations.append("✅ Database quality is excellent - focus on maintaining standards for new data")
        
        return recommendations

def main():
    """Main function to run comprehensive validation"""
    print("🔍 CryoProtect Data Quality Validation Engine")
    print("=" * 50)
    
    validator = DataQualityValidator()
    
    # Run validation
    print("Running comprehensive database validation...")
    report = validator.validate_all_molecules()
    
    if "error" in report:
        print(f"❌ Validation failed: {report['error']}")
        return
    
    # Display summary
    summary = report["summary"]
    print(f"\n📊 VALIDATION SUMMARY")
    print(f"Total molecules validated: {report['total_molecules']}")
    print(f"Average quality score: {summary['overall_quality']['average_score']}/100")
    print(f"Quality distribution:")
    print(f"  • Excellent (90+): {summary['quality_distribution']['excellent_90_plus']}")
    print(f"  • Good (80-89): {summary['quality_distribution']['good_80_89']}")
    print(f"  • Fair (70-79): {summary['quality_distribution']['fair_70_79']}")
    print(f"  • Poor (<70): {summary['quality_distribution']['poor_below_70']}")
    
    print(f"\n🚨 ISSUES DETECTED")
    print(f"Total errors: {summary['issues']['total_errors']}")
    print(f"Total warnings: {summary['issues']['total_warnings']}")
    print(f"Molecules with errors: {summary['issues']['molecules_with_errors']}")
    print(f"Molecules with warnings: {summary['issues']['molecules_with_warnings']}")
    
    # Generate recommendations
    recommendations = validator.get_quality_recommendations(report['validation_results'])
    print(f"\n💡 RECOMMENDATIONS")
    for i, rec in enumerate(recommendations, 1):
        print(f"{i}. {rec}")
    
    print(f"\n✅ Validation complete! Detailed report available in the generated JSON file.")

if __name__ == "__main__":
    main()