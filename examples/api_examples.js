/**
 * CryoProtect API Examples - JavaScript
 * 
 * This file contains examples of how to use the CryoProtect API with JavaScript.
 * It demonstrates common operations such as authentication, working with molecules,
 * mixtures, predictions, RDKit integration, scoring, and export functionality.
 */

// Base URL for the API
const BASE_URL = "http://localhost:5000";

// Authentication token (would be obtained through login)
let TOKEN = null;

// Helper function to get authentication token
function getToken() {
  return TOKEN;
}

// Authentication Examples
async function signUp(email, password) {
  try {
    const response = await fetch(`${BASE_URL}/auth/register`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ email, password })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error signing up');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error signing up:', error.message);
    throw error;
  }
}

async function signIn(email, password) {
  try {
    const response = await fetch(`${BASE_URL}/auth/login`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ email, password })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error signing in');
    }
    
    const data = await response.json();
    TOKEN = data.access_token;
    
    return data;
  } catch (error) {
    console.error('Error signing in:', error.message);
    throw error;
  }
}

async function signOut() {
  try {
    const response = await fetch(`${BASE_URL}/auth/logout`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      }
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error signing out');
    }
    
    TOKEN = null;
    return { success: true };
  } catch (error) {
    console.error('Error signing out:', error.message);
    throw error;
  }
}

// Molecule Examples
async function getMolecules() {
  try {
    const response = await fetch(`${BASE_URL}/api/v1/molecules`, {
      method: 'GET',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      }
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error fetching molecules');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error fetching molecules:', error.message);
    throw error;
  }
}

async function createMolecule(moleculeData) {
  try {
    const response = await fetch(`${BASE_URL}/api/v1/molecules`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify(moleculeData)
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error creating molecule');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error creating molecule:', error.message);
    throw error;
  }
}

// RDKit Integration Examples
async function calculateProperties(moleculeData, inputFormat = 'smiles') {
  try {
    const response = await fetch(`${BASE_URL}/api/v1/rdkit/properties`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        molecule_data: moleculeData,
        input_format: inputFormat
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error calculating properties');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error calculating molecular properties:', error.message);
    throw error;
  }
}

async function generateVisualization(moleculeData, inputFormat = 'smiles', width = 400, height = 300, highlightAtoms = []) {
  try {
    const response = await fetch(`${BASE_URL}/api/v1/rdkit/visualization`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        molecule_data: moleculeData,
        input_format: inputFormat,
        width: width,
        height: height,
        highlight_atoms: highlightAtoms
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error generating visualization');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error generating molecular visualization:', error.message);
    throw error;
  }
}

async function performSubstructureSearch(queryMolData, targetMolData, queryFormat = 'smarts', targetFormat = 'smiles') {
  try {
    const response = await fetch(`${BASE_URL}/api/v1/rdkit/substructure`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        query_mol_data: queryMolData,
        target_mol_data: targetMolData,
        query_format: queryFormat,
        target_format: targetFormat
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error performing substructure search');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error performing substructure search:', error.message);
    throw error;
  }
}

// Scoring Examples
async function scoreMolecule(moleculeData, inputFormat = 'smiles', storeResult = false) {
  try {
    const response = await fetch(`${BASE_URL}/api/v1/scoring/molecules`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        molecule_data: moleculeData,
        input_format: inputFormat,
        store_result: storeResult
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error scoring molecule');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error scoring molecule:', error.message);
    throw error;
  }
}

async function scoreMixtureById(mixtureId, storeResult = true) {
  try {
    const response = await fetch(`${BASE_URL}/api/v1/mixtures/${mixtureId}/score`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        store_result: storeResult
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error scoring mixture');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error scoring mixture:', error.message);
    throw error;
  }
}

async function batchScore(entityType, ids, storeResults = true) {
  try {
    const response = await fetch(`${BASE_URL}/api/v1/scoring/batch`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        entity_type: entityType,
        ids: ids,
        store_results: storeResults
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error performing batch scoring');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Error performing batch scoring:', error.message);
    throw error;
  }
}

// Export and Sharing Examples
async function exportData(dataType, format, id = null, includeRelated = false) {
  try {
    const response = await fetch(`${BASE_URL}/api/v1/export`, {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json',
        'Authorization': `Bearer ${getToken()}`
      },
      body: JSON.stringify({
        data_type: dataType,
        format: format,
        id: id,
        include_related: includeRelated
      })
    });
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Error exporting data');
    }
    
    // Handle file download
    const blob = await response.blob();
    const filename = response.headers.get('Content-Disposition')
      ? response.headers.get('Content-Disposition').split('filename=')[1].replace(/"/g, '')
      : `${dataType}_export.${format}`;
    
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.style.display = 'none';
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    window.URL.revokeObjectURL(url);
    
    return { success: true, filename };
  } catch (error) {
    console.error('Error exporting data:', error.message);
    throw error;
  }
}

// Rate Limiting Examples
async function handleRateLimitedRequest(url, options) {
  try {
    const response = await fetch(url, options);
    
    // Check for rate limit headers
    const rateLimitLimit = response.headers.get('X-RateLimit-Limit');
    const rateLimitRemaining = response.headers.get('X-RateLimit-Remaining');
    const rateLimitReset = response.headers.get('X-RateLimit-Reset');
    
    if (rateLimitRemaining && parseInt(rateLimitRemaining) < 5) {
      console.warn(`Rate limit warning: ${rateLimitRemaining}/${rateLimitLimit} requests remaining. Resets at ${new Date(parseInt(rateLimitReset) * 1000).toLocaleTimeString()}`);
    }
    
    if (response.status === 429) {
      // Rate limit exceeded
      const retryAfter = response.headers.get('Retry-After') || 60;
      console.warn(`Rate limit exceeded. Retrying after ${retryAfter} seconds.`);
      
      // Wait for the specified time and retry
      await new Promise(resolve => setTimeout(resolve, parseInt(retryAfter) * 1000));
      return handleRateLimitedRequest(url, options);
    }
    
    if (!response.ok) {
      const errorData = await response.json();
      throw new Error(errorData.message || 'Request failed');
    }
    
    return await response.json();
  } catch (error) {
    console.error('Request error:', error.message);
    throw error;
  }
}

// Example usage
async function runExamples() {
  console.log("CryoProtect API Examples - JavaScript");
  console.log("=====================================");
  
  // Sign in
  console.log("\nSigning in...");
  try {
    // Replace with your credentials
    const authData = await signIn("user@example.com", "password");
    console.log(`Signed in successfully. Token: ${TOKEN.substring(0, 10)}...`);
  } catch (error) {
    console.error(`Failed to sign in: ${error.message}`);
    return;
  }
  
  // Calculate properties for ethanol
  console.log("\nCalculating properties for ethanol...");
  try {
    const properties = await calculateProperties("CCO");
    console.log(`LogP: ${properties.logp}`);
    console.log(`TPSA: ${properties.tpsa}`);
    console.log(`H-Bond Donors: ${properties.hydrogen_bonding.donors}`);
    console.log(`H-Bond Acceptors: ${properties.hydrogen_bonding.acceptors}`);
  } catch (error) {
    console.error(`Failed to calculate properties: ${error.message}`);
  }
  
  // Score ethanol
  console.log("\nScoring ethanol...");
  try {
    const score = await scoreMolecule("CCO");
    console.log(`Overall Score: ${score.overall_score}`);
    console.log(`Component Scores:`, score.component_scores);
  } catch (error) {
    console.error(`Failed to score molecule: ${error.message}`);
  }
  
  // Sign out
  console.log("\nSigning out...");
  try {
    await signOut();
    console.log("Signed out successfully.");
  } catch (error) {
    console.error(`Failed to sign out: ${error.message}`);
  }
}

// Run examples if this file is executed directly
if (typeof window !== 'undefined') {
  // Browser environment
  document.addEventListener('DOMContentLoaded', () => {
    const runButton = document.getElementById('run-examples');
    if (runButton) {
      runButton.addEventListener('click', runExamples);
    }
  });
} else if (typeof require !== 'undefined') {
  // Node.js environment
  if (require.main === module) {
    runExamples().catch(console.error);
  }
}

// Export functions for use in other modules
if (typeof module !== 'undefined' && typeof module.exports !== 'undefined') {
  module.exports = {
    signUp,
    signIn,
    signOut,
    getMolecules,
    createMolecule,
    calculateProperties,
    generateVisualization,
    performSubstructureSearch,
    scoreMolecule,
    scoreMixtureById,
    batchScore,
    exportData,
    handleRateLimitedRequest
  };
}
