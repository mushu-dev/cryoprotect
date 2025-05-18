# Organism-Based Cryoprotectant Recommendation System

This document outlines the implementation plan for a new feature that allows researchers to suggest biological organisms or use cases for cryoprotectants and receive appropriate recommendations. This system will significantly enhance the platform's scientific value and research utility.

## Architecture Overview

![Architecture Diagram]

The system consists of five main components:

1. **Convex Database** - Stores organism classification and effectiveness data
2. **Flask Backend (Heroku)** - Hosts the API for organism recommendations
3. **RDKit Service (fly.io)** - Provides molecular property calculations and similarity analysis
4. **React Frontend (Netlify)** - User interface for organism selection and recommendations
5. **Recommendation Engine** - Core algorithm that connects organisms to appropriate cryoprotectants

### Data Flow

1. User selects or defines an organism through the frontend interface
2. Request is sent to the Flask backend API
3. Recommendation engine queries the Convex database for:
   - Direct evidence of effectiveness for this organism
   - Similar organisms with known effectiveness data
   - Molecular properties relevant to the organism type
4. RDKit service is consulted for property-based compatibility analysis
5. Results are compiled with scientific reasoning and citations
6. Recommendations are returned to the frontend for display

## Database Schema (Convex)

```typescript
// schema.ts
export default defineSchema({
  // Organism classification table
  organisms: defineTable({
    name: v.string(),
    category: v.string(), // cell_line, primary_cells, tissue, embryo, etc.
    taxonomy: v.object({
      domain: v.string(),
      kingdom: v.string(),
      phylum: v.string(),
      class: v.string(),
      order: v.string(),
      family: v.string(),
      genus: v.string(),
      species: v.string(),
    }),
    sensitivity: v.object({
      osmotic: v.number(), // 1-10 scale
      toxic: v.number(),   // 1-10 scale
      thermal: v.number()  // 1-10 scale
    }),
    membraneProperties: v.object({
      permeability: v.number(),
      fluidityTemp: v.number(),
      cholesterolContent: v.optional(v.number())
    }),
    metadata: v.object({
      addedBy: v.string(),
      dateAdded: v.number(),
      lastUpdated: v.number(),
      references: v.array(v.string())
    })
  }).index("by_category", ["category"]),
  
  // Effectiveness data linking organisms to mixtures
  effectiveness: defineTable({
    mixtureId: v.string(),
    organismId: v.id("organisms"),
    effectivenessScore: v.number(), // 0-1 scale
    metrics: v.object({
      viability: v.number(), // 0-1 scale
      functionalRecovery: v.number(), // 0-1 scale
      structuralIntegrity: v.number(), // 0-1 scale
      longTermStability: v.number() // 0-1 scale
    }),
    protocol: v.object({
      concentrationPercent: v.number(),
      equilibrationTimeMinutes: v.number(),
      coolingRateC_per_min: v.number(),
      thawingMethod: v.string(),
      postThawDilution: v.string()
    }),
    evidence: v.object({
      quality: v.string(), // "direct", "inferred", "theoretical"
      source: v.string(), // "publication", "user_report", "algorithm"
      confidence: v.number(), // 0-1 scale
      citations: v.array(v.string())
    }),
    metadata: v.object({
      addedBy: v.string(),
      dateAdded: v.number(),
      lastUpdated: v.number(),
      notes: v.string()
    })
  }).index("by_organism", ["organismId"])
  .index("by_mixture_organism", ["mixtureId", "organismId"])
  .index("by_effectiveness", ["effectivenessScore"]),
  
  // User feedback on recommendations
  recommendation_feedback: defineTable({
    userId: v.string(),
    recommendationId: v.string(),
    organismId: v.id("organisms"),
    mixtureId: v.string(),
    rating: v.number(), // 1-5 scale
    actualEffectiveness: v.optional(v.number()), // 0-1 scale
    comments: v.optional(v.string()),
    timestamp: v.number()
  }).index("by_user", ["userId"])
  .index("by_recommendation", ["recommendationId"])
});
```

## Backend API Specification

### New Python Modules

#### 1. `organism_recommendation.py`

```python
from typing import List, Dict, Any, Optional
import requests
from database.db import get_db
from models import Mixture, Protocol

class OrganismClassifier:
    """Manages organism classification and taxonomy mapping"""
    
    def classify_organism(self, organism_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Classify an organism based on provided data
        
        Args:
            organism_data: Dict containing organism details
            
        Returns:
            Enriched organism data with classification
        """
        # Implementation details...
        pass
        
    def find_similar_organisms(self, organism_id: str, limit: int = 5) -> List[Dict[str, Any]]:
        """
        Find organisms similar to the given one based on taxonomy and properties
        
        Args:
            organism_id: ID of the target organism
            limit: Maximum number of similar organisms to return
            
        Returns:
            List of similar organisms with similarity scores
        """
        # Implementation details...
        pass


class OrganismRecommendationEngine:
    """Core recommendation engine for organism-based suggestions"""
    
    def __init__(self):
        self.db = get_db()
        self.classifier = OrganismClassifier()
        
    def recommend_for_organism(self, organism_id: str, count: int = 5, 
                             evidence_threshold: float = 0.5) -> List[Dict[str, Any]]:
        """
        Generate cryoprotectant recommendations for a specific organism.
        
        Args:
            organism_id: ID of the organism
            count: Number of recommendations to return
            evidence_threshold: Minimum evidence score (0-1)
            
        Returns:
            List of recommended mixtures with protocols and reasoning
        """
        # 1. Get organism data
        organism = self.db.query_organism(organism_id)
        
        # 2. Check for direct evidence
        direct_matches = self.db.query_effectiveness(
            organism_id=organism_id, 
            min_score=evidence_threshold
        )
        
        recommendations = []
        
        # 3. Use direct evidence if available
        if len(direct_matches) >= count:
            for match in direct_matches[:count]:
                mixture = self.db.get_mixture(match["mixtureId"])
                protocol = self.generate_protocol(organism, mixture)
                
                recommendations.append({
                    "mixture": mixture,
                    "protocol": protocol,
                    "evidence": match["evidence"],
                    "effectiveness": match["effectivenessScore"],
                    "reasoning": self.generate_reasoning(organism, mixture, match)
                })
            return recommendations
        
        # 4. Find similar organisms if direct evidence is insufficient
        similar_organisms = self.classifier.find_similar_organisms(organism_id)
        
        # 5. Get effectiveness data for similar organisms
        for similar in similar_organisms:
            similar_matches = self.db.query_effectiveness(
                organism_id=similar["id"],
                min_score=evidence_threshold * 0.8  # Lower threshold for similar organisms
            )
            
            for match in similar_matches:
                # Adjust effectiveness based on similarity
                adjusted_effectiveness = match["effectivenessScore"] * similar["similarity"]
                
                # Skip if below threshold after adjustment
                if adjusted_effectiveness < evidence_threshold * 0.7:
                    continue
                
                mixture = self.db.get_mixture(match["mixtureId"])
                protocol = self.generate_protocol(organism, mixture, 
                                               base_protocol=match["protocol"])
                
                # Create recommendation with inference from similar organism
                recommendations.append({
                    "mixture": mixture,
                    "protocol": protocol,
                    "evidence": {
                        "quality": "inferred",
                        "source": "similar_organism",
                        "confidence": similar["similarity"],
                        "original_organism": similar["name"],
                        "citations": match["evidence"]["citations"]
                    },
                    "effectiveness": adjusted_effectiveness,
                    "reasoning": self.generate_reasoning(organism, mixture, match, 
                                                      similar_organism=similar)
                })
                
                # Exit if we have enough recommendations
                if len(recommendations) >= count:
                    return sorted(recommendations, key=lambda x: x["effectiveness"], reverse=True)
        
        # 6. If still insufficient, use molecular property-based recommendations
        if len(recommendations) < count:
            property_recommendations = self.recommend_by_properties(organism, 
                                                                 count - len(recommendations))
            recommendations.extend(property_recommendations)
            
        return sorted(recommendations, key=lambda x: x["effectiveness"], reverse=True)
    
    def recommend_by_properties(self, organism: Dict[str, Any], count: int) -> List[Dict[str, Any]]:
        """
        Generate recommendations based on molecular properties suitable for the organism
        
        Args:
            organism: Organism data dictionary
            count: Number of recommendations to return
            
        Returns:
            List of property-based recommendations
        """
        # Call RDKit service for property compatibility analysis
        rdkit_endpoint = "https://your-rdkit-service.fly.io/analyze_compatibility"
        response = requests.post(rdkit_endpoint, json={
            "organism_category": organism["category"],
            "sensitivity": organism["sensitivity"],
            "membrane_properties": organism.get("membraneProperties", {})
        })
        
        if response.status_code != 200:
            return []
            
        compatible_mixtures = response.json()["compatible_mixtures"]
        
        recommendations = []
        for mixture_data in compatible_mixtures[:count]:
            mixture = self.db.get_mixture(mixture_data["mixture_id"])
            protocol = self.generate_protocol(organism, mixture)
            
            recommendations.append({
                "mixture": mixture,
                "protocol": protocol,
                "evidence": {
                    "quality": "theoretical",
                    "source": "property_analysis",
                    "confidence": mixture_data["compatibility_score"],
                    "citations": []
                },
                "effectiveness": mixture_data["predicted_effectiveness"],
                "reasoning": self.generate_property_reasoning(organism, mixture, mixture_data)
            })
            
        return recommendations
    
    def generate_protocol(self, organism: Dict[str, Any], mixture: Dict[str, Any], 
                        base_protocol: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """
        Generate optimized protocol for the organism and mixture
        
        Args:
            organism: Organism data
            mixture: Mixture data
            base_protocol: Optional starting protocol to adjust
            
        Returns:
            Optimized protocol parameters
        """
        # Start with base protocol or default values
        protocol = base_protocol or {
            "concentrationPercent": 10.0,
            "equilibrationTimeMinutes": 30,
            "coolingRateC_per_min": 1.0,
            "thawingMethod": "rapid",
            "postThawDilution": "stepwise"
        }
        
        # Adjust based on organism sensitivity
        if organism["sensitivity"]["osmotic"] > 7:
            # High osmotic sensitivity requires gradual concentration changes
            protocol["concentrationPercent"] = min(protocol["concentrationPercent"], 8.0)
            protocol["equilibrationTimeMinutes"] = max(protocol["equilibrationTimeMinutes"], 45)
            protocol["postThawDilution"] = "stepwise"
            
        if organism["sensitivity"]["toxic"] > 7:
            # High toxicity sensitivity requires lower concentration
            protocol["concentrationPercent"] = min(protocol["concentrationPercent"], 7.0)
            protocol["equilibrationTimeMinutes"] = min(protocol["equilibrationTimeMinutes"], 20)
            
        if organism["sensitivity"]["thermal"] > 7:
            # High thermal sensitivity requires slower cooling
            protocol["coolingRateC_per_min"] = min(protocol["coolingRateC_per_min"], 0.5)
            
        # Adjust based on organism category
        if organism["category"] == "cell_line":
            protocol["thawingMethod"] = "rapid"
        elif organism["category"] == "tissue":
            protocol["thawingMethod"] = "controlled"
            protocol["coolingRateC_per_min"] = min(protocol["coolingRateC_per_min"], 0.3)
        elif organism["category"] == "embryo":
            protocol["thawingMethod"] = "ultrarapid"
            
        return protocol
    
    def generate_reasoning(self, organism: Dict[str, Any], mixture: Dict[str, Any], 
                         effectiveness_data: Dict[str, Any], 
                         similar_organism: Optional[Dict[str, Any]] = None) -> str:
        """
        Generate scientific reasoning for the recommendation
        
        Args:
            organism: Target organism data
            mixture: Recommended mixture
            effectiveness_data: Effectiveness data
            similar_organism: Optional similar organism data for inferred recommendations
            
        Returns:
            Detailed scientific reasoning with citations
        """
        # Implementation details...
        pass
    
    def generate_property_reasoning(self, organism: Dict[str, Any], mixture: Dict[str, Any],
                                  compatibility_data: Dict[str, Any]) -> str:
        """
        Generate reasoning based on molecular property compatibility
        
        Args:
            organism: Target organism data
            mixture: Recommended mixture
            compatibility_data: Property compatibility data
            
        Returns:
            Property-based reasoning text
        """
        # Implementation details...
        pass
```

#### 2. `organism_resources.py`

```python
from flask import request
from flask_restful import Resource, fields, marshal_with
from api_decorators import token_required
from organism_recommendation import OrganismRecommendationEngine

# Marshalling fields for API responses
organism_fields = {
    'id': fields.String,
    'name': fields.String,
    'category': fields.String,
    'taxonomy': fields.Nested({
        'domain': fields.String,
        'kingdom': fields.String,
        'phylum': fields.String,
        'class': fields.String,
        'order': fields.String,
        'family': fields.String,
        'genus': fields.String,
        'species': fields.String
    }),
    'sensitivity': fields.Nested({
        'osmotic': fields.Float,
        'toxic': fields.Float,
        'thermal': fields.Float
    }),
    'membraneProperties': fields.Nested({
        'permeability': fields.Float,
        'fluidityTemp': fields.Float,
        'cholesterolContent': fields.Float
    })
}

protocol_fields = {
    'concentrationPercent': fields.Float,
    'equilibrationTimeMinutes': fields.Integer,
    'coolingRateC_per_min': fields.Float,
    'thawingMethod': fields.String,
    'postThawDilution': fields.String
}

evidence_fields = {
    'quality': fields.String,
    'source': fields.String,
    'confidence': fields.Float,
    'citations': fields.List(fields.String),
    'original_organism': fields.String(attribute=lambda x: x.get('original_organism', ''))
}

recommendation_fields = {
    'mixture': fields.Nested({
        'id': fields.String,
        'name': fields.String,
        'components': fields.List(fields.Nested({
            'name': fields.String,
            'concentration': fields.Float,
            'unit': fields.String
        }))
    }),
    'protocol': fields.Nested(protocol_fields),
    'evidence': fields.Nested(evidence_fields),
    'effectiveness': fields.Float,
    'reasoning': fields.String
}

class OrganismResource(Resource):
    """Resource for managing organism data"""
    
    @token_required
    @marshal_with(organism_fields)
    def get(self, organism_id):
        """
        Get organism details by ID
        """
        # Implementation details...
        pass
        
    @token_required
    @marshal_with(organism_fields)
    def put(self, organism_id):
        """
        Update organism details
        """
        # Implementation details...
        pass
        
    @token_required
    def delete(self, organism_id):
        """
        Delete an organism
        """
        # Implementation details...
        pass


class OrganismsResource(Resource):
    """Resource for managing multiple organisms"""
    
    @token_required
    @marshal_with({'organisms': fields.List(fields.Nested(organism_fields))})
    def get(self):
        """
        Get list of all organisms
        """
        # Implementation details...
        pass
        
    @token_required
    @marshal_with(organism_fields)
    def post(self):
        """
        Create a new organism
        """
        # Implementation details...
        pass


class OrganismRecommendationResource(Resource):
    """Resource for organism-based recommendations"""
    
    @token_required
    @marshal_with({'recommendations': fields.List(fields.Nested(recommendation_fields))})
    def get(self, organism_id):
        """
        Get cryoprotectant recommendations for an organism
        """
        count = request.args.get('count', default=5, type=int)
        evidence_threshold = request.args.get('evidence_threshold', default=0.5, type=float)
        
        engine = OrganismRecommendationEngine()
        recommendations = engine.recommend_for_organism(
            organism_id=organism_id,
            count=count,
            evidence_threshold=evidence_threshold
        )
        
        return {'recommendations': recommendations}


class EffectivenessResource(Resource):
    """Resource for managing effectiveness data"""
    
    @token_required
    def post(self):
        """
        Record effectiveness data for an organism-mixture combination
        """
        # Implementation details...
        pass
        
    @token_required
    def get(self, organism_id=None, mixture_id=None):
        """
        Get effectiveness data, optionally filtered by organism or mixture
        """
        # Implementation details...
        pass


class RecommendationFeedbackResource(Resource):
    """Resource for user feedback on recommendations"""
    
    @token_required
    def post(self):
        """
        Submit feedback on a recommendation
        """
        # Implementation details...
        pass


def register_organism_resources(api):
    """
    Register all organism-related resources with the API
    """
    api.add_resource(OrganismResource, '/api/organisms/<string:organism_id>')
    api.add_resource(OrganismsResource, '/api/organisms')
    api.add_resource(OrganismRecommendationResource, 
                   '/api/recommendations/organism/<string:organism_id>')
    api.add_resource(EffectivenessResource, '/api/effectiveness', 
                   '/api/effectiveness/organism/<string:organism_id>',
                   '/api/effectiveness/mixture/<string:mixture_id>')
    api.add_resource(RecommendationFeedbackResource, '/api/recommendations/feedback')
```

#### 3. Integration with `app.py`

```python
# Add to existing app.py
from api.organism_resources import register_organism_resources

# Register organism resources
register_organism_resources(api)
```

## Frontend Component Requirements

### React Components

#### 1. `OrganismSelector.tsx`

```tsx
import React, { useState, useEffect } from 'react';
import { useQuery, useMutation } from '@tanstack/react-query';
import { fetchOrganisms, createOrganism } from '../api/organisms';
import { Box, TextField, Autocomplete, Button, Dialog, 
         Typography, Slider, Grid, FormControl, InputLabel } from '@mui/material';

interface OrganismSelectorProps {
  onOrganismSelected: (organismId: string) => void;
}

interface Organism {
  id: string;
  name: string;
  category: string;
  // other properties...
}

export const OrganismSelector: React.FC<OrganismSelectorProps> = ({ onOrganismSelected }) => {
  const [selectedOrganism, setSelectedOrganism] = useState<Organism | null>(null);
  const [dialogOpen, setDialogOpen] = useState(false);
  const [newOrganism, setNewOrganism] = useState({
    name: '',
    category: '',
    taxonomy: {
      domain: '',
      kingdom: '',
      phylum: '',
      class: '',
      order: '',
      family: '',
      genus: '',
      species: '',
    },
    sensitivity: {
      osmotic: 5,
      toxic: 5,
      thermal: 5
    },
    membraneProperties: {
      permeability: 5,
      fluidityTemp: 20,
      cholesterolContent: 50
    }
  });

  // Fetch organisms
  const { data: organisms = [], isLoading } = useQuery({
    queryKey: ['organisms'],
    queryFn: fetchOrganisms
  });

  // Create organism mutation
  const createOrganismMutation = useMutation({
    mutationFn: createOrganism,
    onSuccess: (data) => {
      setDialogOpen(false);
      setSelectedOrganism(data);
      onOrganismSelected(data.id);
    }
  });

  // Handle organism selection
  const handleOrganismChange = (_: any, value: Organism | null) => {
    setSelectedOrganism(value);
    if (value) {
      onOrganismSelected(value.id);
    }
  };

  // Handle new organism creation
  const handleCreateOrganism = () => {
    createOrganismMutation.mutate(newOrganism);
  };

  return (
    <Box sx={{ mb: 4 }}>
      <Typography variant="h6" gutterBottom>
        Select Organism
      </Typography>
      
      <Grid container spacing={2}>
        <Grid item xs={9}>
          <Autocomplete
            options={organisms}
            loading={isLoading}
            getOptionLabel={(option) => `${option.name} (${option.category})`}
            value={selectedOrganism}
            onChange={handleOrganismChange}
            renderInput={(params) => (
              <TextField
                {...params}
                label="Organism"
                variant="outlined"
                placeholder="Search organisms..."
              />
            )}
          />
        </Grid>
        <Grid item xs={3}>
          <Button 
            variant="outlined" 
            fullWidth 
            onClick={() => setDialogOpen(true)}
            sx={{ height: '56px' }}
          >
            Add New
          </Button>
        </Grid>
      </Grid>

      {/* New Organism Dialog */}
      <Dialog open={dialogOpen} onClose={() => setDialogOpen(false)} maxWidth="md" fullWidth>
        <Box sx={{ p: 4 }}>
          <Typography variant="h6" gutterBottom>
            Add New Organism
          </Typography>
          
          {/* Basic Info */}
          <Grid container spacing={3} sx={{ mb: 4 }}>
            <Grid item xs={6}>
              <TextField
                fullWidth
                label="Organism Name"
                value={newOrganism.name}
                onChange={(e) => setNewOrganism({...newOrganism, name: e.target.value})}
              />
            </Grid>
            <Grid item xs={6}>
              <FormControl fullWidth>
                <InputLabel>Category</InputLabel>
                <Autocomplete
                  options={['cell_line', 'primary_cells', 'tissue', 'embryo', 'other']}
                  value={newOrganism.category}
                  onChange={(_, value) => setNewOrganism({
                    ...newOrganism, 
                    category: value || ''
                  })}
                  renderInput={(params) => <TextField {...params} label="Category" />}
                />
              </FormControl>
            </Grid>
          </Grid>
          
          {/* Taxonomy */}
          <Typography variant="subtitle1" gutterBottom>
            Taxonomy
          </Typography>
          <Grid container spacing={2} sx={{ mb: 4 }}>
            {/* Taxonomy fields */}
            {Object.keys(newOrganism.taxonomy).map((key) => (
              <Grid item xs={3} key={key}>
                <TextField
                  fullWidth
                  label={key.charAt(0).toUpperCase() + key.slice(1)}
                  value={newOrganism.taxonomy[key as keyof typeof newOrganism.taxonomy]}
                  onChange={(e) => setNewOrganism({
                    ...newOrganism,
                    taxonomy: {
                      ...newOrganism.taxonomy,
                      [key]: e.target.value
                    }
                  })}
                />
              </Grid>
            ))}
          </Grid>
          
          {/* Sensitivity */}
          <Typography variant="subtitle1" gutterBottom>
            Sensitivity Profile (1-10)
          </Typography>
          <Grid container spacing={2} sx={{ mb: 4 }}>
            <Grid item xs={4}>
              <Typography>Osmotic Sensitivity: {newOrganism.sensitivity.osmotic}</Typography>
              <Slider
                min={1}
                max={10}
                step={1}
                value={newOrganism.sensitivity.osmotic}
                onChange={(_, value) => setNewOrganism({
                  ...newOrganism,
                  sensitivity: {
                    ...newOrganism.sensitivity,
                    osmotic: value as number
                  }
                })}
              />
            </Grid>
            <Grid item xs={4}>
              <Typography>Toxic Sensitivity: {newOrganism.sensitivity.toxic}</Typography>
              <Slider
                min={1}
                max={10}
                step={1}
                value={newOrganism.sensitivity.toxic}
                onChange={(_, value) => setNewOrganism({
                  ...newOrganism,
                  sensitivity: {
                    ...newOrganism.sensitivity,
                    toxic: value as number
                  }
                })}
              />
            </Grid>
            <Grid item xs={4}>
              <Typography>Thermal Sensitivity: {newOrganism.sensitivity.thermal}</Typography>
              <Slider
                min={1}
                max={10}
                step={1}
                value={newOrganism.sensitivity.thermal}
                onChange={(_, value) => setNewOrganism({
                  ...newOrganism,
                  sensitivity: {
                    ...newOrganism.sensitivity,
                    thermal: value as number
                  }
                })}
              />
            </Grid>
          </Grid>
          
          {/* Action Buttons */}
          <Box sx={{ display: 'flex', justifyContent: 'flex-end', gap: 2, mt: 4 }}>
            <Button variant="outlined" onClick={() => setDialogOpen(false)}>
              Cancel
            </Button>
            <Button 
              variant="contained" 
              onClick={handleCreateOrganism}
              disabled={!newOrganism.name || !newOrganism.category}
            >
              Create Organism
            </Button>
          </Box>
        </Box>
      </Dialog>
    </Box>
  );
};
```

#### 2. `RecommendationResults.tsx`

```tsx
import React, { useState } from 'react';
import { useQuery } from '@tanstack/react-query';
import { getOrganismRecommendations } from '../api/recommendations';
import { Box, Card, CardContent, Typography, Divider, Chip,
         Accordion, AccordionSummary, AccordionDetails, 
         Table, TableBody, TableCell, TableHead, TableRow,
         CircularProgress, Rating, Button } from '@mui/material';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';

interface RecommendationResultsProps {
  organismId: string | null;
}

export const RecommendationResults: React.FC<RecommendationResultsProps> = ({ organismId }) => {
  const [selectedIndex, setSelectedIndex] = useState<number | null>(null);

  // Fetch recommendations when organism is selected
  const { data: recommendationData, isLoading } = useQuery({
    queryKey: ['recommendations', organismId],
    queryFn: () => organismId ? getOrganismRecommendations(organismId) : null,
    enabled: !!organismId
  });

  const recommendations = recommendationData?.recommendations || [];

  if (!organismId) {
    return (
      <Box sx={{ mt: 4, textAlign: 'center' }}>
        <Typography variant="h6" color="textSecondary">
          Select an organism to see recommendations
        </Typography>
      </Box>
    );
  }

  if (isLoading) {
    return (
      <Box sx={{ mt: 4, display: 'flex', justifyContent: 'center' }}>
        <CircularProgress />
      </Box>
    );
  }

  if (recommendations.length === 0) {
    return (
      <Box sx={{ mt: 4, textAlign: 'center' }}>
        <Typography variant="h6" color="textSecondary">
          No recommendations found for this organism
        </Typography>
      </Box>
    );
  }

  return (
    <Box sx={{ mt: 4 }}>
      <Typography variant="h5" gutterBottom>
        Recommended Cryoprotectants
      </Typography>
      
      {recommendations.map((recommendation, index) => (
        <Card 
          key={index} 
          elevation={2} 
          sx={{ 
            mb: 3, 
            border: selectedIndex === index ? '2px solid #2196f3' : 'none' 
          }}
          onClick={() => setSelectedIndex(index)}
        >
          <CardContent>
            <Box sx={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center' }}>
              <Typography variant="h6">
                {recommendation.mixture.name}
              </Typography>
              <Box sx={{ display: 'flex', alignItems: 'center', gap: 2 }}>
                <Chip 
                  label={`Effectiveness: ${Math.round(recommendation.effectiveness * 100)}%`}
                  color={recommendation.effectiveness > 0.7 ? "success" : 
                         recommendation.effectiveness > 0.4 ? "warning" : "error"}
                />
                <Chip 
                  label={`Evidence: ${recommendation.evidence.quality}`}
                  color={recommendation.evidence.quality === "direct" ? "success" : 
                         recommendation.evidence.quality === "inferred" ? "warning" : "info"}
                />
              </Box>
            </Box>
            
            <Divider sx={{ my: 2 }} />
            
            {/* Mixture Components */}
            <Typography variant="subtitle1" gutterBottom>
              Components:
            </Typography>
            <Table size="small" sx={{ mb: 2 }}>
              <TableHead>
                <TableRow>
                  <TableCell>Component</TableCell>
                  <TableCell align="right">Concentration</TableCell>
                  <TableCell align="right">Unit</TableCell>
                </TableRow>
              </TableHead>
              <TableBody>
                {recommendation.mixture.components.map((component, idx) => (
                  <TableRow key={idx}>
                    <TableCell>{component.name}</TableCell>
                    <TableCell align="right">{component.concentration}</TableCell>
                    <TableCell align="right">{component.unit}</TableCell>
                  </TableRow>
                ))}
              </TableBody>
            </Table>
            
            {/* Protocol */}
            <Typography variant="subtitle1" gutterBottom>
              Recommended Protocol:
            </Typography>
            <Box sx={{ display: 'flex', flexWrap: 'wrap', gap: 2, mb: 2 }}>
              <Chip label={`Concentration: ${recommendation.protocol.concentrationPercent}%`} />
              <Chip label={`Equilibration: ${recommendation.protocol.equilibrationTimeMinutes} min`} />
              <Chip label={`Cooling Rate: ${recommendation.protocol.coolingRateC_per_min}°C/min`} />
              <Chip label={`Thawing: ${recommendation.protocol.thawingMethod}`} />
              <Chip label={`Post-thaw: ${recommendation.protocol.postThawDilution}`} />
            </Box>
            
            {/* Reasoning */}
            <Accordion>
              <AccordionSummary expandIcon={<ExpandMoreIcon />}>
                <Typography>Scientific Reasoning</Typography>
              </AccordionSummary>
              <AccordionDetails>
                <Typography variant="body2" color="textSecondary" sx={{ whiteSpace: 'pre-line' }}>
                  {recommendation.reasoning}
                </Typography>
                
                {recommendation.evidence.citations.length > 0 && (
                  <>
                    <Typography variant="subtitle2" sx={{ mt: 2 }}>Citations:</Typography>
                    <ul>
                      {recommendation.evidence.citations.map((citation, citIdx) => (
                        <li key={citIdx}>
                          <Typography variant="caption">{citation}</Typography>
                        </li>
                      ))}
                    </ul>
                  </>
                )}
              </AccordionDetails>
            </Accordion>
            
            {/* Feedback */}
            <Box sx={{ mt: 2, display: 'flex', alignItems: 'center', gap: 2 }}>
              <Typography variant="body2">Was this recommendation helpful?</Typography>
              <Rating 
                name={`feedback-${index}`} 
                size="small"
                onChange={(_, value) => {
                  // Submit feedback logic
                }}
              />
              <Button variant="outlined" size="small">
                Report Results
              </Button>
            </Box>
          </CardContent>
        </Card>
      ))}
    </Box>
  );
};
```

#### 3. Main Page Component: `organism-recommendations.tsx`

```tsx
import React, { useState } from 'react';
import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { Box, Container, Typography, Paper } from '@mui/material';
import { OrganismSelector } from '../components/OrganismSelector';
import { RecommendationResults } from '../components/RecommendationResults';

const queryClient = new QueryClient();

export default function OrganismRecommendationsPage() {
  const [selectedOrganismId, setSelectedOrganismId] = useState<string | null>(null);

  return (
    <QueryClientProvider client={queryClient}>
      <Container maxWidth="lg">
        <Box sx={{ py: 4 }}>
          <Typography variant="h4" component="h1" gutterBottom>
            Organism-Based Cryoprotectant Recommendations
          </Typography>
          
          <Typography variant="body1" color="textSecondary" paragraph>
            Get tailored cryoprotectant recommendations based on specific organisms or cell types.
            Select an existing organism from our database or create a new entry with detailed 
            characteristics to receive customized recommendations.
          </Typography>
          
          <Paper elevation={2} sx={{ p: 3, mt: 3 }}>
            <OrganismSelector onOrganismSelected={setSelectedOrganismId} />
            <RecommendationResults organismId={selectedOrganismId} />
          </Paper>
        </Box>
      </Container>
    </QueryClientProvider>
  );
}
```

## Integration Plan

### Phase 1: Database Setup and Backend Implementation (2 weeks)

1. **Week 1: Database Schema Implementation**
   - Create Convex database schema
   - Implement initial data migration for existing organisms
   - Set up test data for development

2. **Week 2: Backend Implementation**
   - Develop core recommendation engine
   - Implement API endpoints
   - Write unit tests
   - Integration with RDKit service

### Phase 2: Frontend Development (2 weeks)

1. **Week 3: Component Development**
   - Implement organism selector component
   - Create recommendation results display
   - Develop protocol visualization

2. **Week 4: Integration and Refinement**
   - Connect frontend to backend API
   - Implement user feedback mechanisms
   - Improve UI/UX
   - Add responsive design

### Phase 3: Testing and Deployment (2 weeks)

1. **Week 5: Testing**
   - User acceptance testing
   - Performance optimization
   - Bug fixes and refinements

2. **Week 6: Deployment and Documentation**
   - Deploy to production environment
   - Write user documentation
   - Create scientific documentation

### Dependencies and Prerequisites

- Backend-to-frontend connection stability (flagged as current issue)
- RDKit service availability on fly.io
- Convex database configured for new schema
- Scientific literature database for citations

## Conclusion

The Organism-Based Cryoprotectant Recommendation System will provide significant scientific value to researchers by:

1. **Personalizing recommendations** for specific organisms and biological materials
2. **Leveraging scientific evidence** to generate credible suggestions
3. **Explaining the reasoning** behind each recommendation
4. **Optimizing protocols** based on organism sensitivity profiles
5. **Learning from user feedback** to continuously improve recommendations

This feature represents an important advancement in making cryoprotection research more accessible and effective, allowing researchers to quickly identify appropriate cryoprotectants for their specific biological materials.