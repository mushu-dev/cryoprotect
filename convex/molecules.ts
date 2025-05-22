import { v } from "convex/values";
import { mutation, query } from "./_generated/server";

// Create a new molecule
export const createMolecule = mutation({
  args: {
    name: v.string(),
    pubchemCid: v.optional(v.string()),
    canonicalSmiles: v.optional(v.string()),
    inchiKey: v.optional(v.string()),
    formula: v.optional(v.string()),
    status: v.string(),
  },
  handler: async (ctx, args) => {
    const moleculeId = await ctx.db.insert("molecules", {
      name: args.name,
      pubchemCid: args.pubchemCid,
      canonicalSmiles: args.canonicalSmiles,
      inchiKey: args.inchiKey,
      formula: args.formula,
      status: args.status,
      createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString(),
    });
    
    return moleculeId;
  },
});

// List molecules
export const listMolecules = query({
  args: { limit: v.optional(v.number()) },
  handler: async (ctx, args) => {
    const limit = args.limit || 10;
    return await ctx.db.query("molecules").take(limit);
  },
});

// Get molecule by ID
export const getMolecule = query({
  args: { id: v.id("molecules") },
  handler: async (ctx, args) => {
    return await ctx.db.get(args.id);
  },
});

// Get molecule with all properties
export const getMoleculeWithProperties = query({
  args: { id: v.id("molecules") },
  handler: async (ctx, args) => {
    const molecule = await ctx.db.get(args.id);
    if (!molecule) return null;
    
    // Get molecular properties (get the latest one)
    const properties = await ctx.db
      .query("molecularProperties")
      .withIndex("by_molecule", (q) => q.eq("moleculeId", args.id))
      .order("desc")
      .first();
    
    // Get cryoprotectant scores (get the latest one)
    const scores = await ctx.db
      .query("cryoprotectantScores")
      .withIndex("by_molecule", (q) => q.eq("moleculeId", args.id))
      .order("desc")
      .first();
    
    // Get experimental data
    const experimentalData = await ctx.db
      .query("experimentalData")
      .withIndex("by_molecule", (q) => q.eq("moleculeId", args.id))
      .collect();
    
    return {
      ...molecule,
      properties,
      scores,
      experimentalData,
    };
  },
});

// List molecules with scores for ranking
export const listMoleculesWithScores = query({
  args: { 
    limit: v.optional(v.number()),
    category: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    const limit = args.limit || 10;
    
    // Get molecules
    const molecules = await ctx.db.query("molecules").take(1000); // Get all, then filter
    
    // Get scores for each molecule
    const moleculesWithScores = await Promise.all(
      molecules.map(async (molecule) => {
        const scores = await ctx.db
          .query("cryoprotectantScores")
          .withIndex("by_molecule", (q) => q.eq("moleculeId", molecule._id))
          .order("desc")
          .first();
        
        return {
          ...molecule,
          scores,
        };
      })
    );
    
    // Filter by category if specified
    let filteredMolecules = moleculesWithScores;
    if (args.category) {
      filteredMolecules = moleculesWithScores.filter(
        (mol) => mol.scores?.category === args.category
      );
    }
    
    // Sort by overall score (descending)
    filteredMolecules.sort((a, b) => {
      const scoreA = a.scores?.overallScore || 0;
      const scoreB = b.scores?.overallScore || 0;
      return scoreB - scoreA;
    });
    
    return filteredMolecules.slice(0, limit);
  },
});

// Create molecular properties
export const createMolecularProperties = mutation({
  args: {
    moleculeId: v.id("molecules"),
    molecularWeight: v.optional(v.number()),
    exactMass: v.optional(v.number()),
    logP: v.optional(v.number()),
    tpsa: v.optional(v.number()),
    hbondDonors: v.optional(v.number()),
    hbondAcceptors: v.optional(v.number()),
    rotatableBonds: v.optional(v.number()),
    aromaticRings: v.optional(v.number()),
    aliphaticRings: v.optional(v.number()),
    complexity: v.optional(v.number()),
    heavyAtomCount: v.optional(v.number()),
    morganFingerprint: v.optional(v.string()),
    rdkitFingerprint: v.optional(v.string()),
    calculationVersion: v.optional(v.string()),
    status: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    const propertiesId = await ctx.db.insert("molecularProperties", {
      moleculeId: args.moleculeId,
      molecularWeight: args.molecularWeight,
      exactMass: args.exactMass,
      logP: args.logP,
      tpsa: args.tpsa,
      hbondDonors: args.hbondDonors,
      hbondAcceptors: args.hbondAcceptors,
      rotatableBonds: args.rotatableBonds,
      aromaticRings: args.aromaticRings,
      aliphaticRings: args.aliphaticRings,
      complexity: args.complexity,
      heavyAtomCount: args.heavyAtomCount,
      morganFingerprint: args.morganFingerprint,
      rdkitFingerprint: args.rdkitFingerprint,
      calculatedAt: new Date().toISOString(),
      calculationVersion: args.calculationVersion || "v1.0",
      status: args.status || "calculated",
    });
    
    return propertiesId;
  },
});

// Create cryoprotectant scores
export const createCryoprotectantScores = mutation({
  args: {
    moleculeId: v.id("molecules"),
    glassTempScore: v.optional(v.number()),
    viscosityScore: v.optional(v.number()),
    permeabilityScore: v.optional(v.number()),
    toxicityScore: v.optional(v.number()),
    overallScore: v.number(),
    category: v.string(),
    scoringAlgorithmVersion: v.optional(v.string()),
    confidence: v.optional(v.number()),
  },
  handler: async (ctx, args) => {
    const scoresId = await ctx.db.insert("cryoprotectantScores", {
      moleculeId: args.moleculeId,
      glassTempScore: args.glassTempScore,
      viscosityScore: args.viscosityScore,
      permeabilityScore: args.permeabilityScore,
      toxicityScore: args.toxicityScore,
      overallScore: args.overallScore,
      category: args.category,
      scoringAlgorithmVersion: args.scoringAlgorithmVersion || "v1.0",
      calculatedAt: new Date().toISOString(),
      confidence: args.confidence,
    });
    
    return scoresId;
  },
});

// Create calculation job
export const createCalculationJob = mutation({
  args: {
    moleculeId: v.id("molecules"),
    jobType: v.string(),
    parameters: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    const jobId = await ctx.db.insert("calculationJobs", {
      moleculeId: args.moleculeId,
      jobType: args.jobType,
      status: "pending",
      parameters: args.parameters,
      createdAt: new Date().toISOString(),
    });
    
    return jobId;
  },
});

// Update calculation job status
export const updateCalculationJob = mutation({
  args: {
    jobId: v.id("calculationJobs"),
    status: v.string(),
    result: v.optional(v.string()),
    errorMessage: v.optional(v.string()),
  },
  handler: async (ctx, args) => {
    const job = await ctx.db.get(args.jobId);
    if (!job) throw new Error("Job not found");
    
    const updates: any = {
      status: args.status,
    };
    
    if (args.status === "running" && !job.startedAt) {
      updates.startedAt = new Date().toISOString();
    }
    
    if (args.status === "completed" || args.status === "failed") {
      updates.completedAt = new Date().toISOString();
    }
    
    if (args.result) {
      updates.result = args.result;
    }
    
    if (args.errorMessage) {
      updates.errorMessage = args.errorMessage;
    }
    
    await ctx.db.patch(args.jobId, updates);
    
    return args.jobId;
  },
});

// Get pending calculation jobs
export const getPendingJobs = query({
  args: { jobType: v.optional(v.string()) },
  handler: async (ctx, args) => {
    let query = ctx.db.query("calculationJobs").withIndex("by_status", (q) => q.eq("status", "pending"));
    
    if (args.jobType) {
      // Filter by job type after querying
      const allPendingJobs = await query.collect();
      return allPendingJobs.filter(job => job.jobType === args.jobType);
    }
    
    return await query.collect();
  },
});

// List all molecular properties for analysis
export const listAllMolecularProperties = query({
  args: { limit: v.optional(v.number()) },
  handler: async (ctx, args) => {
    const limit = args.limit || 1000; // Default to large limit for analysis
    return await ctx.db.query("molecularProperties").take(limit);
  },
});

// List all molecules with their properties for validation
export const listMoleculesWithProperties = query({
  args: { limit: v.optional(v.number()) },
  handler: async (ctx, args) => {
    const limit = args.limit || 1000; // Default to large limit for validation
    
    // Get all molecules
    const molecules = await ctx.db.query("molecules").take(limit);
    
    // Get properties and scores for each molecule
    const moleculesWithData = await Promise.all(
      molecules.map(async (molecule) => {
        // Get molecular properties (get the latest one)
        const properties = await ctx.db
          .query("molecularProperties")
          .withIndex("by_molecule", (q) => q.eq("moleculeId", molecule._id))
          .order("desc")
          .first();
        
        // Get cryoprotectant scores (get the latest one)
        const scores = await ctx.db
          .query("cryoprotectantScores")
          .withIndex("by_molecule", (q) => q.eq("moleculeId", molecule._id))
          .order("desc")
          .first();
        
        // Flatten the structure for easier validation
        return {
          _id: molecule._id,
          name: molecule.name,
          pubchemCid: molecule.pubchemCid,
          smiles: molecule.canonicalSmiles,
          inchiKey: molecule.inchiKey,
          formula: molecule.formula,
          status: molecule.status,
          createdAt: molecule.createdAt,
          updatedAt: molecule.updatedAt,
          // Molecular properties (flattened)
          molecularWeight: properties?.molecularWeight,
          exactMass: properties?.exactMass,
          logP: properties?.logP,
          tpsa: properties?.tpsa,
          hbondDonors: properties?.hbondDonors,
          hbondAcceptors: properties?.hbondAcceptors,
          rotatableBonds: properties?.rotatableBonds,
          aromaticRings: properties?.aromaticRings,
          aliphaticRings: properties?.aliphaticRings,
          complexity: properties?.complexity,
          heavyAtomCount: properties?.heavyAtomCount,
          morganFingerprint: properties?.morganFingerprint,
          rdkitFingerprint: properties?.rdkitFingerprint,
          // Cryoprotectant scores (flattened)
          glassTempScore: scores?.glassTempScore,
          viscosityScore: scores?.viscosityScore,
          permeabilityScore: scores?.permeabilityScore,
          toxicityScore: scores?.toxicityScore,
          overallScore: scores?.overallScore,
          category: scores?.category,
          confidence: scores?.confidence,
        };
      })
    );
    
    return moleculesWithData;
  },
});