/**
 * Sample Data Generator for Enhanced Schema
 * 
 * This module provides utilities for generating realistic sample data
 * for the enhanced schema tables. This is used for testing, development,
 * and creating demonstration datasets.
 */

import { Id } from "../_generated/dataModel";
import { faker } from "@faker-js/faker";
import { nanoid } from "nanoid";

// Define reusable scientific constants for realistic data
const COMMON_CRYOPROTECTANTS = [
  "Glycerol", "DMSO", "Ethylene glycol", "Propylene glycol", 
  "Trehalose", "Sucrose", "Mannitol", "Polyvinylpyrrolidone"
];

const TISSUE_TYPES = [
  "Hepatocytes", "Oocytes", "Embryos", "Sperm", 
  "Stem cells", "Fibroblasts", "Erythrocytes", "Neurons"
];

const SPECIES = [
  "Human", "Mouse", "Rat", "Bovine", 
  "Porcine", "Rabbit", "Zebrafish", "Drosophila"
];

const VERIFICATION_STATUSES = ["pending", "verified", "rejected"];
const PROTOCOL_TYPES = ["slow_freezing", "vitrification", "controlled_rate_freezing", "ultra-rapid_freezing"];
const VIABILITY_METHODS = ["trypan_blue", "fluorescent_dye", "metabolic_assay", "flow_cytometry", "impedance_based"];

/**
 * Generates random scientific value within biologically plausible ranges
 */
function randomScientificValue(min: number, max: number, precision: number = 2): number {
  const value = faker.number.float({ min, max, precision });
  return Number(value.toFixed(precision));
}

/**
 * Generate sample lab verification data
 */
export function generateLabVerification(
  experimentId: Id<"enhancedExperiments">,
  verifierId: Id<"users">
) {
  const now = Date.now();
  const verificationDate = faker.date.recent({ days: 14 }).getTime();
  
  return {
    experimentId,
    verificationStatus: faker.helpers.arrayElement(VERIFICATION_STATUSES),
    verifierId,
    equipmentUsed: faker.helpers.arrayElement([
      "Olympus BX51 Microscope", 
      "Leica DMi8 Microscope", 
      "Zeiss LSM 880 Confocal", 
      "Thermo Fisher Flow Cytometer"
    ]),
    verificationDate,
    comments: faker.datatype.boolean() ? faker.lorem.paragraph(2) : undefined,
    evidenceUrls: faker.datatype.boolean() ? 
      Array.from({ length: faker.number.int({ min: 1, max: 3 }) }, () => 
        faker.image.url()
      ) : undefined,
    reviewedProtocolSteps: faker.datatype.boolean() ?
      Array.from({ length: faker.number.int({ min: 1, max: 5 }) }, () => 
        nanoid(8)
      ) : undefined,
    qualityScore: faker.datatype.boolean() ? 
      faker.number.int({ min: 1, max: 10 }) : undefined,
    createdAt: now,
    updatedAt: now
  };
}

/**
 * Generate sample enhanced experiment results
 */
export function generateEnhancedExperimentResult(
  experimentId: Id<"enhancedExperiments"> | string,
  tissueTypeId: Id<"tissueTypes"> | string,
  moleculeId?: Id<"molecules">,
  mixtureId?: Id<"mixtures">
) {
  const now = Date.now();
  
  // Generate realistic viability percentage - bimodal distribution for realistic outcomes
  let viabilityPercentage: number;
  if (faker.datatype.boolean({ probability: 0.7 })) {
    // Successful experiment (higher viability)
    viabilityPercentage = randomScientificValue(65, 95);
  } else {
    // Less successful experiment (lower viability)
    viabilityPercentage = randomScientificValue(20, 60);
  }
  
  return {
    experimentId,
    tissueTypeId,
    moleculeId,
    mixtureId,
    concentration: randomScientificValue(0.5, 15),
    concentrationUnit: faker.helpers.arrayElement(["% w/v", "% v/v", "M", "mM"]),
    
    // Enhanced viability measurements
    viabilityPercentage,
    viabilityMethod: faker.helpers.arrayElement(VIABILITY_METHODS),
    recoveryRate: randomScientificValue(viabilityPercentage * 0.8, viabilityPercentage * 1.1, 1),
    recoveryRateUnit: "%/hour",
    functionalityScore: randomScientificValue(1, 10, 1),
    functionalAssay: faker.helpers.arrayElement([
      "ATP production", "Protein synthesis", "Membrane integrity", 
      "Enzyme activity", "Oxygen consumption"
    ]),
    integrityMeasure: randomScientificValue(60, 95),
    integrityMethod: faker.helpers.arrayElement([
      "Permeability assay", "Electron microscopy", "Impedance measurement", 
      "Fluorescent dye exclusion"
    ]),
    
    // Cryopreservation-specific metrics
    iceFormationObserved: faker.datatype.boolean({ probability: 0.3 }),
    postThawMorphologyScore: faker.number.int({ min: 1, max: 10 }),
    stressResponseMarkers: faker.helpers.arrayElements(
      ["HSP70", "ROS", "Apoptosis", "DNA damage", "Membrane lipid peroxidation"],
      { min: 0, max: 3 }
    ),
    timeToRecovery: randomScientificValue(2, 36, 1),
    
    // Additional metadata
    uncertainty: {
      viabilityPercentage: {
        type: "standard_deviation",
        value: randomScientificValue(1, 5)
      }
    },
    resultDetails: {
      temperature: randomScientificValue(-196, 37),
      humidity: randomScientificValue(30, 95),
      plateId: `P${faker.string.numeric(4)}`
    },
    notes: faker.datatype.boolean({ probability: 0.3 }) ? 
      faker.lorem.paragraph() : undefined,
    createdBy: undefined,
    createdAt: now,
    updatedAt: now
  };
}

/**
 * Generate sample cryopreservation protocol
 */
export function generateCryoProtocol(creatorId?: Id<"users">) {
  const now = Date.now();
  const protocolType = faker.helpers.arrayElement(PROTOCOL_TYPES);
  
  // Generate realistic protocol parameters based on type
  let coolingRate, warmingRate, holdTemperature, cpaEquilibrationTime;
  
  if (protocolType === "slow_freezing") {
    coolingRate = randomScientificValue(-0.5, -3);
    warmingRate = randomScientificValue(10, 30);
    holdTemperature = randomScientificValue(-80, -20);
    cpaEquilibrationTime = randomScientificValue(10, 30);
  } else if (protocolType === "vitrification") {
    coolingRate = randomScientificValue(-1000, -10000);
    warmingRate = randomScientificValue(1000, 5000);
    holdTemperature = -196;
    cpaEquilibrationTime = randomScientificValue(5, 15);
  } else if (protocolType === "controlled_rate_freezing") {
    coolingRate = randomScientificValue(-0.3, -1);
    warmingRate = randomScientificValue(5, 15);
    holdTemperature = randomScientificValue(-60, -30);
    cpaEquilibrationTime = randomScientificValue(15, 45);
  } else {
    coolingRate = randomScientificValue(-5000, -20000);
    warmingRate = randomScientificValue(2000, 8000);
    holdTemperature = -196;
    cpaEquilibrationTime = randomScientificValue(3, 8);
  }
  
  // Generate steps based on protocol type
  const stepsCount = faker.number.int({ min: 3, max: 8 });
  const steps = Array.from({ length: stepsCount }, (_, i) => {
    const isCritical = i === 1 || i === stepsCount - 2; // Make a couple steps critical
    const step = {
      id: nanoid(8),
      name: faker.helpers.arrayElement([
        `Add cryoprotectant ${i+1}`,
        `Cool to ${randomScientificValue(-100, -20)}°C`,
        `Hold at ${randomScientificValue(-80, -30)}°C`,
        "Transfer to storage container",
        `Add ${faker.helpers.arrayElement(COMMON_CRYOPROTECTANTS)}`,
        "Seal container",
        "Plunge in liquid nitrogen",
        "Equilibrate with CPA"
      ]),
      description: faker.lorem.sentence(),
      parameters: {
        concentration: randomScientificValue(1, 15),
        exposureTime: randomScientificValue(1, 30)
      },
      duration: randomScientificValue(1, 60),
      durationUnit: "minutes",
      temperature: randomScientificValue(-196, 37),
      temperatureUnit: "°C",
      
      // Step-specific cryopreservation fields
      rampRate: randomScientificValue(-10, -0.1),
      rampRateUnit: "°C/min",
      pressureApplied: faker.datatype.boolean() ? randomScientificValue(1, 5) : undefined,
      pressureUnit: "atm",
      substancesAdded: faker.helpers.arrayElements(COMMON_CRYOPROTECTANTS, { min: 0, max: 2 }),
      equipmentRequired: faker.helpers.arrayElements([
        "Controlled-rate freezer", "LN2 container", "Cryo vials", 
        "Transfer pipettes", "Temperature logger", "Microscope"
      ], { min: 1, max: 3 }),
      criticalStep: isCritical,
      qualityControlChecks: isCritical ? 
        ["Visual inspection", "Temperature verification", "Timing verification"] : undefined
    };
    
    return step;
  });
  
  return {
    name: `${faker.helpers.arrayElement(["Standard", "Modified", "Optimized", "Experimental"])} ${protocolType} for ${faker.helpers.arrayElement(TISSUE_TYPES)}`,
    description: faker.lorem.paragraph(),
    
    // Cryopreservation-specific fields
    protocolType,
    coolingRate,
    coolingRateUnit: "°C/min",
    warmingRate,
    warmingRateUnit: "°C/min",
    holdTemperature,
    holdDuration: randomScientificValue(10, 120),
    holdDurationUnit: "minutes",
    cryoprotectantAdditionMethod: faker.helpers.arrayElement([
      "step-wise", "single-step", "gradient", "perfusion"
    ]),
    preFreezingTreatment: faker.datatype.boolean() ? 
      faker.helpers.arrayElement([
        "Sucrose loading", "Cell washing", "Medium replacement", 
        "Protein supplementation", "None"
      ]) : undefined,
    postThawingTreatment: faker.datatype.boolean() ?
      faker.helpers.arrayElement([
        "Stepwise dilution", "Centrifugation washing", "Membrane stabilization", 
        "Culture recovery", "None"
      ]) : undefined,
    storageTemperature: faker.helpers.arrayElement([-196, -150, -80, -20]),
    storageContainerType: faker.helpers.arrayElement([
      "Cryovial", "Straw", "Cryobag", "CryoBox", "Dewar"
    ]),
    cpaEquilibrationTime,
    
    steps,
    
    version: faker.number.int({ min: 1, max: 5 }),
    parentId: undefined,
    category: faker.helpers.arrayElement([
      "Research", "Clinical", "Standard", "Experimental", 
      faker.helpers.arrayElement(TISSUE_TYPES)
    ]),
    isTemplate: faker.datatype.boolean({ probability: 0.3 }),
    parameters: {
      targetCellType: faker.helpers.arrayElement(TISSUE_TYPES),
      expectedViability: randomScientificValue(70, 90),
      difficultyLevel: faker.helpers.arrayElement(["Beginner", "Intermediate", "Advanced", "Expert"])
    },
    
    validationStatus: faker.helpers.arrayElement(["unvalidated", "validated", "failed_validation"]),
    publications: faker.datatype.boolean() ? 
      Array.from({ length: faker.number.int({ min: 1, max: 3 }) }, () => 
        `https://doi.org/10.1016/j.cryobiol.${faker.string.numeric(4)}.${faker.string.numeric(2)}`
      ) : undefined,
    successRate: randomScientificValue(50, 95),
    
    createdBy: creatorId,
    createdAt: now,
    updatedAt: now,
    public: faker.datatype.boolean({ probability: 0.7 })
  };
}

/**
 * Generate sample tissue type with enhanced fields
 */
export function generateTissueType(creatorId?: Id<"users">) {
  const now = Date.now();
  const species = faker.helpers.arrayElement(SPECIES);
  const cellType = faker.helpers.arrayElement(TISSUE_TYPES);
  
  return {
    name: `${species} ${cellType}`,
    description: faker.lorem.sentence(),
    species,
    taxonomyId: species === "Human" ? 9606 : faker.number.int({ min: 8000, max: 12000 }),
    
    // Enhanced biospecimen fields
    cellType,
    tissueOrigin: faker.helpers.arrayElement([
      "Liver", "Brain", "Blood", "Bone marrow", "Kidney", 
      "Reproductive tract", "Skin", "Muscle"
    ]),
    cellDensity: randomScientificValue(1, 10),
    cellDensityUnit: faker.helpers.arrayElement(["million/mL", "thousand/μL"]),
    passageNumber: faker.datatype.boolean() ? faker.number.int({ min: 1, max: 20 }) : undefined,
    cellDiameter: randomScientificValue(5, 50),
    cellDiameterUnit: "μm",
    waterContent: randomScientificValue(60, 90),
    lipidContent: randomScientificValue(1, 20),
    osmolality: randomScientificValue(280, 340),
    osmolalityUnit: "mOsm/kg",
    preparationMethod: faker.helpers.arrayElement([
      "Enzymatic digestion", "Mechanical dissociation", 
      "Density centrifugation", "Flow sorting", "Outgrowth"
    ]),
    storageConditions: {
      temperature: -196,
      medium: faker.helpers.arrayElement([
        "DMEM + 10% FBS", "PBS", "Cell freezing medium", 
        "Customized medium", "Commercial cryomedium"
      ]),
      container: faker.helpers.arrayElement([
        "Cryovial", "Straw", "Bag", "Plate"
      ])
    },
    cryopreservationHistory: faker.datatype.boolean({ probability: 0.2 }),
    
    properties: {
      growth_rate: randomScientificValue(0.2, 1.5),
      membrane_permeability: randomScientificValue(1, 10),
      freezing_point: randomScientificValue(-0.5, -3)
    },
    createdBy: creatorId,
    createdAt: now,
    updatedAt: now
  };
}

/**
 * Generate sample biospecimen data
 */
export function generateBiospecimen(tissueTypeId: Id<"tissueTypes">) {
  const now = Date.now();
  const collectionDate = faker.date.past({ years: 2 }).getTime();
  const processingDate = new Date(collectionDate + (1000 * 60 * 60 * faker.number.int({ min: 1, max: 48 }))).getTime();
  
  return {
    tissueTypeId,
    identifier: `B${faker.string.alphanumeric(8).toUpperCase()}`,
    donorId: `D${faker.string.numeric(6)}`,
    donorAge: faker.number.int({ min: 18, max: 85 }),
    donorSex: faker.helpers.arrayElement(["M", "F"]),
    donorSpecies: faker.helpers.arrayElement(SPECIES),
    collectionDate,
    processingMethod: faker.helpers.arrayElement([
      "Fresh isolation", "Biopsy", "Surgical specimen", 
      "Cadaveric collection", "Donation"
    ]),
    processingDate,
    qualityScore: faker.number.int({ min: 1, max: 10 }),
    viabilityAtCollection: randomScientificValue(80, 98),
    storageLocation: `${faker.helpers.arrayElement(["A", "B", "C", "D"])}-${faker.number.int({ min: 1, max: 10 })}-${faker.number.int({ min: 1, max: 100 })}`,
    freezeThawCycles: faker.number.int({ min: 0, max: 3 }),
    clinicalDiagnosis: faker.datatype.boolean() ? 
      faker.helpers.arrayElement([
        "Healthy", "Type 2 Diabetes", "Hypertension", 
        "Cancer", "None", "Unknown"
      ]) : undefined,
    consentInformation: "Full research consent obtained",
    metadata: {
      collection_center: faker.company.name(),
      transport_time: faker.number.int({ min: 1, max: 24 }),
      special_handling: faker.datatype.boolean()
    },
    createdAt: now,
    updatedAt: now
  };
}

/**
 * Generate sample cryoprotectant effectiveness data
 */
export function generateCryoprotectantEffectiveness(
  moleculeId?: Id<"molecules">,
  mixtureId?: Id<"mixtures">,
  tissueTypeId?: Id<"tissueTypes">,
  experimentId?: Id<"enhancedExperiments">
) {
  if (!moleculeId && !mixtureId) {
    throw new Error("Either moleculeId or mixtureId must be provided");
  }
  
  const now = Date.now();
  
  // Generate overall effectiveness score - weighted toward higher scores for known cryoprotectants
  let effectivenessScore: number;
  if (faker.datatype.boolean({ probability: 0.7 })) {
    // Known effective agent
    effectivenessScore = randomScientificValue(65, 95);
  } else {
    // Less effective agent
    effectivenessScore = randomScientificValue(20, 60);
  }
  
  return {
    moleculeId,
    mixtureId,
    tissueTypeId,
    
    // Overall effectiveness metrics
    effectivenessScore,
    effectiveConcentration: randomScientificValue(1, 15),
    effectiveConcentrationUnit: faker.helpers.arrayElement(["% w/v", "% v/v", "M", "mM"]),
    
    // Mechanism-specific metrics
    membranePermeability: randomScientificValue(0.5, 10),
    membranePermeabilityUnit: "μm/s",
    glassTransitionTemp: randomScientificValue(-120, -80),
    toxicityThreshold: randomScientificValue(5, 25),
    toxicityThresholdUnit: "% w/v",
    osmoticTolerance: randomScientificValue(20, 50),
    iceInhibitionCapacity: randomScientificValue(0.1, 0.9),
    
    // Experimental context
    experimentalConditions: {
      minTemperature: randomScientificValue(-196, -20),
      exposureTime: randomScientificValue(1, 60),
      cellType: faker.helpers.arrayElement(TISSUE_TYPES)
    },
    measuredAt: randomScientificValue(-20, 25),
    measuredAtUnit: "°C",
    freezingRate: randomScientificValue(-10, -0.1),
    freezingRateUnit: "°C/min",
    
    // Data provenance
    source: undefined,
    experimentId,
    methodologyDescription: faker.lorem.paragraph(),
    confidenceLevel: randomScientificValue(0.6, 0.99, 2),
    literatureReference: faker.datatype.boolean() ? 
      `https://doi.org/10.1016/j.cryobiol.${faker.string.numeric(4)}.${faker.string.numeric(2)}` : undefined,
    
    createdAt: now,
    updatedAt: now
  };
}

/**
 * Generate a complete set of sample data
 */
export async function generateCompleteSampleDataset(db: any, count: number = 10) {
  console.log(`Generating sample dataset with ${count} records per type...`);
  
  // Generate users
  const userIds = [];
  for (let i = 0; i < 5; i++) {
    const userId = await db.insert("users", {
      email: faker.internet.email(),
      name: faker.person.fullName(),
      role: faker.helpers.arrayElement(["admin", "scientist", "viewer"]),
      lastLogin: faker.date.recent().getTime(),
      createdAt: Date.now(),
      updatedAt: Date.now()
    });
    userIds.push(userId);
  }
  
  // Generate tissue types
  const tissueTypeIds = [];
  for (let i = 0; i < count; i++) {
    const tissueTypeId = await db.insert(
      "tissueTypes", 
      generateTissueType(faker.helpers.arrayElement(userIds))
    );
    tissueTypeIds.push(tissueTypeId);
  }
  
  // Generate biospecimens
  const biospecimenIds = [];
  for (let i = 0; i < count; i++) {
    const biospecimenId = await db.insert(
      "biospecimens", 
      generateBiospecimen(faker.helpers.arrayElement(tissueTypeIds))
    );
    biospecimenIds.push(biospecimenId);
  }
  
  // Generate protocols
  const protocolIds = [];
  for (let i = 0; i < count; i++) {
    const protocolId = await db.insert(
      "protocols", 
      generateCryoProtocol(faker.helpers.arrayElement(userIds))
    );
    protocolIds.push(protocolId);
  }
  
  // Generate experiments (simplified for this example)
  const experimentIds = [];
  for (let i = 0; i < count; i++) {
    const experimentId = await db.insert("enhancedExperiments", {
      title: `Experiment ${i+1}`,
      description: faker.lorem.sentence(),
      experimentTypeId: faker.helpers.arrayElement([
        "cryopreservation", "viability_test", "recovery_assessment"
      ]),
      protocolId: faker.helpers.arrayElement(protocolIds),
      datePerformed: faker.date.recent({ days: 90 }).getTime(),
      status: faker.helpers.arrayElement(["planned", "in-progress", "completed", "failed"]),
      createdBy: faker.helpers.arrayElement(userIds),
      createdAt: Date.now(),
      updatedAt: Date.now(),
      public: true,
      version: 1
    });
    experimentIds.push(experimentId);
  }
  
  // Generate lab verifications
  for (let i = 0; i < Math.floor(count * 0.7); i++) {
    await db.insert(
      "labVerifications", 
      generateLabVerification(
        faker.helpers.arrayElement(experimentIds),
        faker.helpers.arrayElement(userIds)
      )
    );
  }
  
  // Generate experiment results
  for (let i = 0; i < count * 2; i++) {
    await db.insert(
      "enhancedExperimentResults", 
      generateEnhancedExperimentResult(
        faker.helpers.arrayElement(experimentIds),
        faker.helpers.arrayElement(tissueTypeIds)
      )
    );
  }
  
  // Generate cryoprotectant effectiveness (simplified without real molecules)
  for (let i = 0; i < count; i++) {
    // Note: In a real implementation, we would use actual molecule IDs
    await db.insert(
      "cryoprotectantEffectiveness",
      generateCryoprotectantEffectiveness(
        undefined,  // moleculeId 
        undefined,  // mixtureId
        faker.helpers.arrayElement(tissueTypeIds),
        faker.helpers.arrayElement(experimentIds)
      )
    );
  }
  
  console.log("Sample dataset generation complete!");
  return {
    userIds,
    tissueTypeIds,
    biospecimenIds,
    protocolIds,
    experimentIds
  };
}