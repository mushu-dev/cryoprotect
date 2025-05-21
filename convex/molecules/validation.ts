/**
 * Validation functions for molecule data
 */

import { CreateMoleculeInput, UpdateMoleculeInput } from "./types";

/**
 * Validate a molecule name
 * @param name - The molecule name to validate
 * @returns True if the name is valid
 */
export const isValidMoleculeName = (name: string): boolean => {
  return name.trim().length > 0;
};

/**
 * Validate a PubChem CID
 * @param cid - The PubChem CID to validate
 * @returns True if the CID is valid
 */
export const isValidPubchemCid = (cid: string): boolean => {
  return /^\d+$/.test(cid);
};

/**
 * Validate an InChI Key
 * @param inchiKey - The InChI Key to validate
 * @returns True if the InChI Key is valid
 */
export const isValidInchiKey = (inchiKey: string): boolean => {
  return /^[A-Z]{14}-[A-Z]{10}-[A-Z]$/.test(inchiKey);
};

/**
 * Validate a SMILES string
 * @param smiles - The SMILES string to validate
 * @returns True if the SMILES string is valid (basic check)
 */
export const isValidSmiles = (smiles: string): boolean => {
  // Basic SMILES validation - could be enhanced with more sophisticated checks
  return smiles.trim().length > 0 && /^[A-Za-z0-9@\[\]\(\)\\\/.=#$:*+\-]+$/.test(smiles);
};

/**
 * Validate a chemical formula
 * @param formula - The chemical formula to validate
 * @returns True if the formula is valid
 */
export const isValidFormula = (formula: string): boolean => {
  // Basic formula validation - checks for letter-number patterns like C6H12O6
  return /^([A-Z][a-z]?\d*)+$/.test(formula);
};

/**
 * Validate a molecule creation input
 * @param input - The molecule creation input to validate
 * @returns An object with validation result and optional error message
 */
export const validateCreateMoleculeInput = (
  input: CreateMoleculeInput
): { valid: boolean; error?: string } => {
  if (!isValidMoleculeName(input.name)) {
    return { valid: false, error: "Invalid molecule name" };
  }

  if (input.pubchemCid && !isValidPubchemCid(input.pubchemCid)) {
    return { valid: false, error: "Invalid PubChem CID format" };
  }

  if (input.inchiKey && !isValidInchiKey(input.inchiKey)) {
    return { valid: false, error: "Invalid InChI Key format" };
  }

  if (input.canonicalSmiles && !isValidSmiles(input.canonicalSmiles)) {
    return { valid: false, error: "Invalid SMILES format" };
  }

  if (input.formula && !isValidFormula(input.formula)) {
    return { valid: false, error: "Invalid chemical formula format" };
  }

  return { valid: true };
};

/**
 * Validate a molecule update input
 * @param input - The molecule update input to validate
 * @returns An object with validation result and optional error message
 */
export const validateUpdateMoleculeInput = (
  input: UpdateMoleculeInput
): { valid: boolean; error?: string } => {
  if (input.name && !isValidMoleculeName(input.name)) {
    return { valid: false, error: "Invalid molecule name" };
  }

  if (input.pubchemCid && !isValidPubchemCid(input.pubchemCid)) {
    return { valid: false, error: "Invalid PubChem CID format" };
  }

  if (input.inchiKey && !isValidInchiKey(input.inchiKey)) {
    return { valid: false, error: "Invalid InChI Key format" };
  }

  if (input.canonicalSmiles && !isValidSmiles(input.canonicalSmiles)) {
    return { valid: false, error: "Invalid SMILES format" };
  }

  if (input.formula && !isValidFormula(input.formula)) {
    return { valid: false, error: "Invalid chemical formula format" };
  }
  
  if (input.status && !["active", "deprecated", "consolidated"].includes(input.status)) {
    return { valid: false, error: "Invalid status value" };
  }

  return { valid: true };
};