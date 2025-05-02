"""
Search terms for ChEMBL data import.

This module defines search terms for finding cryoprotectant-related compounds
in the ChEMBL database.
"""

CRYOPROTECTANT_CATEGORIES = [
    {
        "category": "polyols",
        "terms": ["glycerol", "sorbitol", "mannitol", "xylitol", "erythritol"]
    },
    {
        "category": "sugars",
        "terms": ["trehalose", "sucrose", "glucose", "lactose", "maltose", "dextran"]
    },
    {
        "category": "alcohols",
        "terms": ["ethanol", "methanol", "propanol", "butanol", "ethylene glycol", "propylene glycol"]
    },
    {
        "category": "amides",
        "terms": ["formamide", "acetamide", "urea", "hydroxyethyl starch"]
    },
    {
        "category": "sulfoxides",
        "terms": ["dimethyl sulfoxide", "DMSO"]
    },
    {
        "category": "amino_acids",
        "terms": ["glycine", "alanine", "proline", "glutamine", "lysine", "arginine"]
    },
]

def get_all_search_terms():
    """Get flattened list of all search terms."""
    terms = []
    for category in CRYOPROTECTANT_CATEGORIES:
        terms.extend(category["terms"])
    return terms

def get_terms_by_category(category):
    """Get search terms for a specific category."""
    for cat in CRYOPROTECTANT_CATEGORIES:
        if cat["category"] == category:
            return cat["terms"]
    return []