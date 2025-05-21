# CryoProtect Enhanced Schema Usage Guide

This guide explains how to use the enhanced Convex schema we've implemented for the CryoProtect application. The enhanced schema provides comprehensive support for scientific cryopreservation research with specialized fields and relationships.

## Table of Contents

1. [Overview](#overview)
2. [Enhanced Schema Features](#enhanced-schema-features)
3. [Setting Up Your Environment](#setting-up-your-environment)
4. [Populating the Database](#populating-the-database)
5. [API Usage Examples](#api-usage-examples)
6. [Migration from Supabase](#migration-from-supabase)
7. [Scientific Query Examples](#scientific-query-examples)
8. [Best Practices](#best-practices)

## Overview

The enhanced Convex schema has been designed to support comprehensive cryopreservation research with specialized data structures for:

- Lab verification and experimental reproducibility
- Detailed viability and recovery metrics
- Cryopreservation-specific protocol parameters
- Comprehensive biospecimen management
- Cryoprotectant effectiveness tracking

This document provides guidance on how to leverage these enhancements in your application.

## Enhanced Schema Features

### Lab Verification System

The lab verification system tracks the verification status of experimental results, ensuring reproducibility and quality control.

**Key Features:**
- Verification status tracking (pending, verified, rejected)
- Equipment and methodology documentation
- Evidence and protocol step verification
- Quality scoring system

### Enhanced Viability Measurements

Specialized fields for tracking cell viability, recovery, and cryopreservation-specific metrics.

**Key Features:**
- Method-specific viability measurements
- Recovery rate and functionality scoring
- Cell integrity metrics
- Cryopreservation-specific observations (ice formation, morphology, etc.)

### Cryopreservation Protocol Specifics

Detailed protocol fields capturing the unique parameters critical to cryopreservation research.

**Key Features:**
- Cooling and warming rate specification
- Cryoprotectant addition methodology
- Pre-freezing and post-thawing treatments
- Critical step identification and quality controls

### Biospecimen Management

Comprehensive tracking of biospecimen characteristics and individual specimen samples.

**Key Features:**
- Detailed tissue and cell type information
- Donor and preparation metadata
- Storage and handling history
- Specialized physical properties (osmolality, water content, etc.)

### Cryoprotectant Effectiveness Metrics

Specialized metrics for tracking cryoprotectant action and effectiveness.

**Key Features:**
- Overall effectiveness scoring
- Mechanism-specific metrics (membrane permeability, glass transition temps, etc.)
- Experimental context documentation
- Literature reference integration

## Setting Up Your Environment

### Prerequisites

1. Node.js 14+ and npm
2. Convex CLI (`npm install -g convex`)
3. Access to your Convex deployment

### Configuration

1. Set up your Convex URL:

```bash
export CONVEX_URL="https://your-project-id.convex.cloud"
# or for Windows
set CONVEX_URL=https://your-project-id.convex.cloud
```

2. Initialize your Convex project if not already done:

```bash
npx convex init
```

## Populating the Database

You can populate your database with sample data for development and testing using our sample data generator.

### Sample Data Population

```bash
# Navigate to your project directory
cd /path/to/cryoprotect

# Run the sample data population script (generates 10 records per type)
node convex/scripts/populateSampleData.js 10

# Or force without confirmation
node convex/scripts/populateSampleData.js 20 --force
```

## API Usage Examples

### Lab Verification Operations

```typescript
// Get verification status for an experiment
const verification = await client.query(api.labVerifications.getByExperimentId, {
  experimentId: "abc123"
});

// Create a new verification
const verificationId = await client.mutation(api.labVerifications.create, {
  experimentId: "abc123", 
  verificationStatus: "pending",
  verifierId: "user123",
  equipmentUsed: "Olympus BX51 Microscope",
  comments: "Initial verification pending protocol review"
});

// Update verification status
await client.mutation(api.labVerifications.update, {
  id: verificationId,
  verificationStatus: "verified",
  qualityScore: 8
});

// Get verification statistics
const stats = await client.query(api.labVerifications.getStats);
console.log(`Verification success rate: ${stats.verifiedPercentage}%`);
```

### Working with Cryopreservation Protocols

```typescript
// Create a new slow-freezing protocol
const protocolId = await client.mutation(api.protocols.create, {
  name: "Slow-freezing protocol for human oocytes",
  description: "Standard slow-freezing protocol optimized for human oocytes",
  protocolType: "slow_freezing",
  coolingRate: -0.3,
  coolingRateUnit: "°C/min",
  warmingRate: 10,
  warmingRateUnit: "°C/min",
  holdTemperature: -80,
  cryoprotectantAdditionMethod: "step-wise",
  // ... additional fields
  steps: [
    {
      id: "step1",
      name: "Equilibration",
      description: "Equilibrate cells with cryoprotectant",
      duration: 10,
      durationUnit: "minutes",
      temperature: 4,
      temperatureUnit: "°C",
      substancesAdded: ["DMSO", "Sucrose"],
      criticalStep: true
    },
    // ... more steps
  ],
  createdAt: Date.now(),
  updatedAt: Date.now(),
  public: true
});

// Get protocols by type
const protocols = await client.query(api.protocols.getByType, {
  type: "vitrification"
});
```

### Recording Experimental Results

```typescript
// Create an experiment result with enhanced viability data
const resultId = await client.mutation(api.enhancedExperimentResults.create, {
  experimentId: "exp123",
  tissueTypeId: "tissue456",
  concentration: 10,
  concentrationUnit: "% v/v",
  
  // Enhanced viability measurements
  viabilityPercentage: 78.5,
  viabilityMethod: "trypan_blue",
  recoveryRate: 65.3,
  recoveryRateUnit: "%",
  functionalityScore: 7,
  functionalAssay: "ATP production",
  
  // Cryopreservation observations
  iceFormationObserved: false,
  postThawMorphologyScore: 8,
  stressResponseMarkers: ["HSP70", "ROS"],
  timeToRecovery: 24,
  
  createdAt: Date.now(),
  updatedAt: Date.now()
});

// Query results by viability range
const highViabilityResults = await client.query(api.enhancedExperimentResults.getByViabilityRange, {
  min: 75,
  max: 100
});
```

## Migration from Supabase

If you're migrating from a Supabase implementation, you can use our migration utilities.

```bash
# Run the migration script
node convex/scripts/populateSampleData.js --migrate

# Set environment variables for migration
export SUPABASE_URL="https://your-project.supabase.co"
export SUPABASE_KEY="your-anon-key"
```

You can also use the programmatic API for more controlled migration:

```javascript
const { ConvexClient } = require("convex/browser");
const { runFullMigration } = require("./convex/utils/migrationUtility");

const client = new ConvexClient(process.env.CONVEX_URL);

// Run the migration
const migrationContext = await runFullMigration(
  client,
  process.env.SUPABASE_URL,
  process.env.SUPABASE_KEY,
  { logLevel: 'info', batchSize: 50 }
);

console.log(migrationContext.stats);
```

## Scientific Query Examples

Here are some examples of scientifically relevant queries enabled by our enhanced schema:

### Compare Viability Across Protocols

```typescript
// Get viability results grouped by protocol
const viabilityByProtocol = await client.query(api.analytics.compareViabilityByProtocol);

// Result format:
// [
//   {
//     protocolName: "Slow-freezing Protocol A",
//     protocolType: "slow_freezing",
//     avgViability: 72.3,
//     sampleSize: 15
//   },
//   {
//     protocolName: "Vitrification Protocol B",
//     protocolType: "vitrification",
//     avgViability: 85.7,
//     sampleSize: 12
//   }
// ]
```

### Analyze Cryoprotectant Effectiveness

```typescript
// Get ranked cryoprotectants by tissue type
const rankedCryoprotectants = await client.query(api.cryoprotectantEffectiveness.getRankedByTissueType, {
  tissueTypeId: "tissue123"
});

// Result format:
// [
//   {
//     moleculeName: "DMSO",
//     effectivenessScore: 87.2,
//     optimalConcentration: "10% v/v",
//     keyMechanisms: ["membrane permeability", "ice inhibition"]
//   },
//   // ... more ranked cryoprotectants
// ]
```

### Track Protocol Verification Success

```typescript
// Get verification statistics by protocol type
const verificationStats = await client.query(api.labVerifications.getStatsByProtocolType);

// Result format:
// {
//   "slow_freezing": {
//     totalCount: 34,
//     verifiedCount: 28,
//     verifiedPercentage: 82.4
//   },
//   "vitrification": {
//     totalCount: 42,
//     verifiedCount: 38,
//     verifiedPercentage: 90.5
//   }
// }
```

## Best Practices

### Data Entry and Validation

1. **Standardized Units**: Always specify units when entering numerical values
2. **Method Documentation**: Record the method used for measurements like viability
3. **Comprehensive Protocol Steps**: Create detailed protocol steps with all parameters
4. **Evidence Links**: Provide URLs to evidence files for verifications
5. **Uncertainty Tracking**: Include uncertainty measurements where applicable

### Scientific Workflows

1. **Protocol-First Approach**: Define protocols before creating experiments
2. **Verification Workflow**: Implement a peer verification process for critical results
3. **Biospecimen Tracking**: Maintain detailed biospecimen records for reproducibility
4. **Literature Connection**: Link protocols to published literature when available
5. **Complete Methodology**: Document all methodology details in protocol steps

### Schema Evolution Strategy

As your research needs evolve, follow these principles for schema evolution:

1. **Backward Compatibility**: Maintain compatibility with existing data
2. **Field Documentation**: Document the scientific meaning of all fields
3. **Controlled Vocabulary**: Use consistent terminology for enum-like fields
4. **Versioning**: Version your protocols and schemas appropriately
5. **Migration Planning**: Plan data migrations carefully when making breaking changes

## Conclusion

The enhanced Convex schema provides comprehensive support for cryopreservation research with specialized data structures. By following this guide, you can leverage these enhancements to build scientifically rigorous applications that support advanced research workflows.

For more details, refer to:
- [Implementation Plan](/convex/ENHANCED_SCHEMA_IMPLEMENTATION_PLAN.md)
- [Schema Gaps Analysis](/convex/SCHEMA_GAPS_ANALYSIS.md)
- [Schema Enhancements Implementation](/convex/SCHEMA_ENHANCEMENTS_IMPLEMENTATION.md)