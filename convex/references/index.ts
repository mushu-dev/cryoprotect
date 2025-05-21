/**
 * Export all reference-related functionality
 */

// Export from crossReferences.ts
export { 
  createCrossReference,
  getCrossReference,
  updateCrossReference,
  deleteCrossReference,
  listCrossReferences,
  getMoleculeCrossReferences,
  batchAddCrossReferences
} from "./crossReferences";

// Export from synonyms.ts
export {
  createSynonym,
  getSynonym,
  updateSynonym,
  deleteSynonym,
  listSynonyms,
  getMoleculeSynonyms,
  searchBySynonym,
  batchAddSynonyms
} from "./synonyms";

// Export types
export type * from "./types";

// Export helpers
export {
  formatCrossReferenceUrl,
  expandCrossReferenceWithMolecule,
  expandSynonymWithMolecule,
  groupSynonymsByType,
  formatSynonym
} from "./helpers";