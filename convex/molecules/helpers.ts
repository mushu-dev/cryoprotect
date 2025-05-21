import { DatabaseReader } from "./_generated/server";
import { Id } from "./_generated/dataModel";

/**
 * Helper functions for working with molecules
 */

/**
 * Check if a molecule exists by ID
 */
export async function doesMoleculeExist(
  db: DatabaseReader,
  id: Id<"molecules">
): Promise<boolean> {
  const molecule = await db.get(id);
  return molecule !== null;
}

/**
 * Check if a molecule with the given name already exists
 */
export async function doesNameExist(
  db: DatabaseReader, 
  name: string
): Promise<boolean> {
  const molecules = await db
    .query("molecules")
    .filter(q => q.eq(q.field("name"), name))
    .take(1);
  
  return molecules.length > 0;
}

/**
 * Find a molecule by its PubChem CID
 */
export async function findByPubchemCid(
  db: DatabaseReader,
  pubchemCid: string
): Promise<Id<"molecules"> | null> {
  const molecules = await db
    .query("molecules")
    .withIndex("by_pubchemCid", q => q.eq("pubchemCid", pubchemCid))
    .take(1);
  
  return molecules.length > 0 ? molecules[0]._id : null;
}

/**
 * Find a molecule by its InChI Key
 */
export async function findByInchiKey(
  db: DatabaseReader,
  inchiKey: string
): Promise<Id<"molecules"> | null> {
  const molecules = await db
    .query("molecules")
    .withIndex("by_inchiKey", q => q.eq("inchiKey", inchiKey))
    .take(1);
  
  return molecules.length > 0 ? molecules[0]._id : null;
}

/**
 * Get active molecules count
 */
export async function getActiveMoleculesCount(
  db: DatabaseReader
): Promise<number> {
  const molecules = await db
    .query("molecules")
    .filter(q => q.eq(q.field("status"), "active"))
    .collect();
  
  return molecules.length;
}

/**
 * Format molecule data for display
 */
export function formatMoleculeForDisplay(molecule: any) {
  return {
    id: molecule._id,
    name: molecule.name,
    formula: molecule.formula || "Not available",
    pubchemCid: molecule.pubchemCid || "Not available",
    smiles: molecule.canonicalSmiles || "Not available",
    inchiKey: molecule.inchiKey || "Not available",
    status: molecule.status,
    lastUpdated: new Date(molecule.updatedAt).toLocaleString()
  };
}
