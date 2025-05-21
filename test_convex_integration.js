// Simple script to test Convex integration
// Run with: node test_convex_integration.js

const { ConvexHttpClient } = require("convex/browser");

// Get Convex URL from environment or use default
const convexUrl = process.env.CONVEX_URL || "https://hallowed-malamute-424.convex.cloud";

async function testConvexIntegration() {
  console.log(`Testing Convex integration with URL: ${convexUrl}`);
  
  // Create Convex HTTP client
  const client = new ConvexHttpClient(convexUrl);
  
  try {
    // Query the molecules table directly since we don't have the generated API yet
    console.log("Fetching molecules directly...");
    const molecules = await client.query("api/query", { 
      table: "molecules",
      limit: 10
    });
    
    console.log(`Successfully fetched ${molecules.length} molecules:`);
    
    // For each molecule, display basic information
    for (let i = 0; i < molecules.length; i++) {
      const molecule = molecules[i];
      console.log(`${i + 1}. ${molecule.name} (${molecule.formula || 'No formula'})`);
      console.log(`   ID: ${molecule._id}`);
      console.log(`   SMILES: ${molecule.canonicalSmiles || 'N/A'}`);
      
      // Fetch properties for this molecule
      try {
        const properties = await client.query("api/query", {
          table: "molecularProperties",
          filters: { moleculeId: molecule._id }
        });
        
        if (properties && properties.length > 0) {
          console.log(`   Properties (${properties.length}):`);
          
          // For each property, try to get the property type
          for (const property of properties) {
            try {
              // Get property type for better display
              const propertyType = await client.query("api/query", {
                table: "propertyTypes",
                filters: { _id: property.propertyTypeId }
              });
              
              const typeName = propertyType[0]?.displayName || 'Unknown Property';
              console.log(`   - ${typeName}: ${property.value} ${property.units || ''}`);
            } catch (propError) {
              console.log(`   - Property ${property.propertyTypeId}: ${property.value}`);
            }
          }
        }
      } catch (propError) {
        console.log(`   Error fetching properties: ${propError.message}`);
      }
      
      console.log(""); // Empty line between molecules
    }
    
    // Test successful
    console.log("\nCONVEX INTEGRATION TEST: SUCCESS");
    return true;
  } catch (error) {
    console.error("Error testing Convex integration:");
    console.error(error);
    console.log("\nCONVEX INTEGRATION TEST: FAILED");
    return false;
  }
}

// Run the test
testConvexIntegration()
  .then(success => {
    process.exit(success ? 0 : 1);
  })
  .catch(err => {
    console.error("Unexpected error:", err);
    process.exit(1);
  });