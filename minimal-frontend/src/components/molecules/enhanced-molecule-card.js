/**
 * Enhanced Molecule Card Component
 * Displays molecule with enhanced properties, scores, and real-time updates
 */
import React from 'react';
import { Card } from '../ui/card';
import { Badge } from '../ui/badge';

// Icons
const ChemicalIcon = () => (
  <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" 
    stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" 
    className="h-4 w-4">
    <path d="M9 3 5 7l4 4"/>
    <path d="m15 3 4 4-4 4"/>
    <path d="M9 19c-4.3 0-8-1.7-8-4"/>
    <path d="M15 19c4.3 0 8-1.7 8-4"/>
    <path d="M9 11h6"/>
  </svg>
);

const ScoreIcon = () => (
  <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" 
    stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" 
    className="h-4 w-4">
    <polygon points="12 2 15.09 8.26 22 9.27 17 14.14 18.18 21.02 12 17.77 5.82 21.02 7 14.14 2 9.27 8.91 8.26 12 2"/>
  </svg>
);

const PropertiesIcon = () => (
  <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" 
    stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" 
    className="h-4 w-4">
    <rect width="7" height="7" x="3" y="3" rx="1"/>
    <rect width="7" height="7" x="14" y="3" rx="1"/>
    <rect width="7" height="7" x="14" y="14" rx="1"/>
    <rect width="7" height="7" x="3" y="14" rx="1"/>
  </svg>
);

const QualityIcon = () => (
  <svg xmlns="http://www.w3.org/2000/svg" width="24" height="24" viewBox="0 0 24 24" fill="none" 
    stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round" 
    className="h-4 w-4">
    <path d="M9 12l2 2 4-4"/>
    <path d="M21 12c-1 0-3-1-3-3s2-3 3-3 3 1 3 3-2 3-3 3"/>
    <path d="M3 12c1 0 3-1 3-3s-2-3-3-3-3 1-3 3 2 3 3 3"/>
    <path d="M3 12c0 5.657 4.343 10 10 10 5.657 0 10-4.343 10-10"/>
  </svg>
);

/**
 * Get category color based on cryoprotectant category
 */
const getCategoryColor = (category) => {
  switch (category?.toLowerCase()) {
    case 'penetrating': return 'bg-blue-100 text-blue-800 border-blue-200';
    case 'non_penetrating': return 'bg-green-100 text-green-800 border-green-200';
    case 'sugar': return 'bg-purple-100 text-purple-800 border-purple-200';
    case 'polymer': return 'bg-orange-100 text-orange-800 border-orange-200';
    case 'alcohol': return 'bg-pink-100 text-pink-800 border-pink-200';
    default: return 'bg-gray-100 text-gray-800 border-gray-200';
  }
};

/**
 * Get score color based on score value
 */
const getScoreColor = (score) => {
  if (score >= 90) return 'text-green-600 font-semibold';
  if (score >= 80) return 'text-blue-600 font-semibold';
  if (score >= 70) return 'text-yellow-600 font-semibold';
  if (score >= 60) return 'text-orange-600 font-semibold';
  return 'text-red-600 font-semibold';
};

/**
 * Enhanced Molecule Card Component
 */
export function EnhancedMoleculeCard({ molecule, onClick, showDetails = false }) {
  if (!molecule) {
    return (
      <Card className="p-4 border-dashed border-gray-300">
        <div className="text-center text-gray-500">
          <ChemicalIcon className="mx-auto mb-2 opacity-50" />
          <p>Molecule data unavailable</p>
        </div>
      </Card>
    );
  }

  const handleClick = () => {
    if (onClick) {
      onClick(molecule);
    }
  };

  // Calculate completeness score
  const requiredFields = ['molecularWeight', 'logP', 'tpsa', 'exactMass', 'complexity'];
  const completedFields = requiredFields.filter(field => molecule[field] != null).length;
  const completeness = Math.round((completedFields / requiredFields.length) * 100);

  return (
    <Card 
      className={`p-4 hover:shadow-lg transition-all duration-200 cursor-pointer border ${
        molecule.cryoprotectantScore >= 80 ? 'border-blue-200 bg-blue-50/30' : 'border-gray-200'
      }`}
      onClick={handleClick}
    >
      {/* Header */}
      <div className="flex items-start justify-between mb-3">
        <div className="flex-1 min-w-0">
          <h3 className="font-semibold text-lg text-gray-900 truncate">
            {molecule.name}
          </h3>
          {molecule.formula && (
            <p className="text-sm text-gray-600 font-mono">
              {molecule.formula}
            </p>
          )}
        </div>
        
        {molecule.cryoprotectantScore && (
          <div className="flex-shrink-0 ml-2">
            <div className="flex items-center space-x-1">
              <ScoreIcon />
              <span className={`text-lg font-bold ${getScoreColor(molecule.cryoprotectantScore)}`}>
                {Math.round(molecule.cryoprotectantScore)}
              </span>
            </div>
          </div>
        )}
      </div>

      {/* Categories and Status */}
      <div className="flex flex-wrap gap-2 mb-3">
        {molecule.cryoprotectantCategory && (
          <Badge 
            variant="outline" 
            className={getCategoryColor(molecule.cryoprotectantCategory)}
          >
            {molecule.cryoprotectantCategory.replace('_', ' ')}
          </Badge>
        )}
        
        {molecule.status && (
          <Badge variant="outline">
            {molecule.status}
          </Badge>
        )}
        
        {completeness === 100 && (
          <Badge variant="outline" className="bg-green-100 text-green-800 border-green-200">
            <QualityIcon className="mr-1" />
            Complete
          </Badge>
        )}
      </div>

      {/* Key Properties */}
      <div className="grid grid-cols-2 gap-3 text-sm">
        {molecule.molecularWeight && (
          <div>
            <span className="text-gray-600">Mol. Weight:</span>
            <span className="ml-1 font-medium">
              {Math.round(molecule.molecularWeight * 100) / 100} Da
            </span>
          </div>
        )}
        
        {molecule.logP && (
          <div>
            <span className="text-gray-600">LogP:</span>
            <span className="ml-1 font-medium">
              {Math.round(molecule.logP * 100) / 100}
            </span>
          </div>
        )}
        
        {molecule.tpsa && (
          <div>
            <span className="text-gray-600">TPSA:</span>
            <span className="ml-1 font-medium">
              {Math.round(molecule.tpsa * 100) / 100} Ų
            </span>
          </div>
        )}
        
        {molecule.heavyAtomCount && (
          <div>
            <span className="text-gray-600">Heavy Atoms:</span>
            <span className="ml-1 font-medium">
              {molecule.heavyAtomCount}
            </span>
          </div>
        )}
      </div>

      {/* Enhanced Details (shown when showDetails is true) */}
      {showDetails && (
        <div className="mt-4 pt-3 border-t border-gray-200">
          <div className="grid grid-cols-1 gap-2 text-sm">
            {/* Cryoprotectant Scores */}
            {molecule.glassTempScore && (
              <div className="flex justify-between">
                <span className="text-gray-600">Glass Formation:</span>
                <span className={`font-medium ${getScoreColor(molecule.glassTempScore)}`}>
                  {Math.round(molecule.glassTempScore)}%
                </span>
              </div>
            )}
            
            {molecule.permeabilityScore && (
              <div className="flex justify-between">
                <span className="text-gray-600">Permeability:</span>
                <span className={`font-medium ${getScoreColor(molecule.permeabilityScore)}`}>
                  {Math.round(molecule.permeabilityScore)}%
                </span>
              </div>
            )}
            
            {molecule.toxicityScore && (
              <div className="flex justify-between">
                <span className="text-gray-600">Safety:</span>
                <span className={`font-medium ${getScoreColor(molecule.toxicityScore)}`}>
                  {Math.round(molecule.toxicityScore)}%
                </span>
              </div>
            )}
            
            {/* Additional Properties */}
            {molecule.complexity && (
              <div className="flex justify-between">
                <span className="text-gray-600">Complexity:</span>
                <span className="font-medium">
                  {Math.round(molecule.complexity * 100) / 100}
                </span>
              </div>
            )}
            
            {molecule.exactMass && (
              <div className="flex justify-between">
                <span className="text-gray-600">Exact Mass:</span>
                <span className="font-medium">
                  {Math.round(molecule.exactMass * 10000) / 10000} Da
                </span>
              </div>
            )}
          </div>
          
          {/* Data Quality Indicator */}
          <div className="mt-3 pt-2 border-t border-gray-100">
            <div className="flex items-center justify-between text-xs">
              <span className="text-gray-500">Data Completeness:</span>
              <div className="flex items-center space-x-1">
                <div className="w-16 h-2 bg-gray-200 rounded-full overflow-hidden">
                  <div 
                    className={`h-full transition-all duration-300 ${
                      completeness >= 80 ? 'bg-green-500' : 
                      completeness >= 60 ? 'bg-yellow-500' : 'bg-red-500'
                    }`}
                    style={{ width: `${completeness}%` }}
                  />
                </div>
                <span className="font-medium text-gray-700">
                  {completeness}%
                </span>
              </div>
            </div>
          </div>
        </div>
      )}

      {/* PubChem CID (if available) */}
      {molecule.pubchemCid && (
        <div className="mt-3 pt-2 border-t border-gray-100">
          <span className="text-xs text-gray-500">
            PubChem CID: {molecule.pubchemCid}
          </span>
        </div>
      )}
    </Card>
  );
}

/**
 * Enhanced Molecule Card Skeleton for loading states
 */
export function EnhancedMoleculeCardSkeleton() {
  return (
    <Card className="p-4 animate-pulse">
      <div className="flex items-start justify-between mb-3">
        <div className="flex-1">
          <div className="h-5 bg-gray-200 rounded w-3/4 mb-2"></div>
          <div className="h-4 bg-gray-200 rounded w-1/2"></div>
        </div>
        <div className="w-8 h-8 bg-gray-200 rounded"></div>
      </div>
      
      <div className="flex gap-2 mb-3">
        <div className="h-5 bg-gray-200 rounded w-16"></div>
        <div className="h-5 bg-gray-200 rounded w-12"></div>
      </div>
      
      <div className="grid grid-cols-2 gap-3">
        <div className="h-4 bg-gray-200 rounded"></div>
        <div className="h-4 bg-gray-200 rounded"></div>
        <div className="h-4 bg-gray-200 rounded"></div>
        <div className="h-4 bg-gray-200 rounded"></div>
      </div>
    </Card>
  );
}