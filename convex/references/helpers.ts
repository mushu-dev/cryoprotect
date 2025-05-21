/**
 * Helper functions for molecular cross references and synonyms
 */

import { Id } from "../_generated/dataModel";
import { QueryCtx } from "../_generated/server";
import { 
  MoleculeCrossReference, 
  MoleculeSynonym,
  CrossReferenceWithMolecule,
  SynonymWithMolecule
} from "./types";

/**
 * Format a cross reference URL based on database and identifier
 */
export function formatCrossReferenceUrl(
  databaseName: string,
  identifier: string,
  existingUrl?: string
): string {
  // If a URL is already provided and valid, use it
  if (existingUrl && isValidUrl(existingUrl)) {
    return existingUrl;
  }
  
  // Common database URL patterns
  const urlPatterns: Record<string, string> = {
    "PubChem": "https://pubchem.ncbi.nlm.nih.gov/compound/{id}",
    "ChEMBL": "https://www.ebi.ac.uk/chembl/compound_report_card/{id}/",
    "DrugBank": "https://go.drugbank.com/drugs/{id}",
    "CAS": "https://commonchemistry.cas.org/detail?cas_rn={id}",
    "ZINC": "https://zinc.docking.org/substances/{id}/",
    "ChEBI": "https://www.ebi.ac.uk/chebi/searchId.do?chebiId={id}",
    "KEGG": "https://www.genome.jp/dbget-bin/www_bget?cpd:{id}"
  };
  
  // Normalize database name for lookup
  const normalizedName = databaseName.trim();
  
  // If we have a pattern for this database, use it
  if (normalizedName in urlPatterns) {
    return urlPatterns[normalizedName].replace("{id}", identifier);
  }
  
  // Otherwise, return an empty string
  return "";
}

/**
 * Expand a cross reference with its molecule data
 */
export async function expandCrossReferenceWithMolecule(
  ctx: QueryCtx,
  crossRef: MoleculeCrossReference
): Promise<CrossReferenceWithMolecule> {
  const molecule = await ctx.db.get(crossRef.moleculeId);
  
  if (!molecule) {
    throw new Error(`Molecule with ID ${crossRef.moleculeId} not found`);
  }
  
  return {
    ...crossRef,
    molecule: {
      _id: molecule._id,
      name: molecule.name,
      formula: molecule.formula
    }
  };
}

/**
 * Expand a synonym with its molecule data
 */
export async function expandSynonymWithMolecule(
  ctx: QueryCtx,
  synonym: MoleculeSynonym
): Promise<SynonymWithMolecule> {
  const molecule = await ctx.db.get(synonym.moleculeId);
  
  if (!molecule) {
    throw new Error(`Molecule with ID ${synonym.moleculeId} not found`);
  }
  
  return {
    ...synonym,
    molecule: {
      _id: molecule._id,
      name: molecule.name,
      formula: molecule.formula
    }
  };
}

/**
 * Group synonyms by type
 */
export function groupSynonymsByType(
  synonyms: MoleculeSynonym[]
): Record<string, MoleculeSynonym[]> {
  const grouped: Record<string, MoleculeSynonym[]> = {};
  
  for (const synonym of synonyms) {
    const type = synonym.type || "other";
    
    if (!grouped[type]) {
      grouped[type] = [];
    }
    
    grouped[type].push(synonym);
  }
  
  return grouped;
}

/**
 * Format a synonym based on its type
 */
export function formatSynonym(synonym: MoleculeSynonym): string {
  if (synonym.type) {
    return `${synonym.name} (${synonym.type})`;
  }
  
  return synonym.name;
}

/**
 * Helper to validate URL format
 */
function isValidUrl(url: string): boolean {
  try {
    // Simple URL validation
    new URL(url);
    return true;
  } catch (error) {
    return false;
  }
}