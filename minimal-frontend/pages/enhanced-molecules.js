/**
 * Enhanced Molecules Page
 * Showcases the enhanced Convex integration with real-time updates
 */
import React from 'react';
import { EnhancedMoleculesList } from '../src/components/molecules/enhanced-molecules-list';
import { useSimpleQualityMetrics } from '../src/hooks/use-simple-enhanced-molecules';
import { Card } from '../src/components/ui/card';
import { Badge } from '../src/components/ui/badge';

/**
 * Quality Summary Component
 */
function QualitySummary() {
  const { data: qualityMetrics, isLoading } = useSimpleQualityMetrics();

  if (isLoading) {
    return (
      <Card className="p-4 animate-pulse">
        <div className="h-4 bg-gray-200 rounded w-1/3 mb-2"></div>
        <div className="grid grid-cols-3 gap-4">
          <div className="h-8 bg-gray-200 rounded"></div>
          <div className="h-8 bg-gray-200 rounded"></div>
          <div className="h-8 bg-gray-200 rounded"></div>
        </div>
      </Card>
    );
  }

  if (!qualityMetrics) return null;

  return (
    <Card className="p-4 bg-gradient-to-r from-blue-50 to-green-50 border-blue-200">
      <h3 className="font-semibold mb-3 text-gray-900">Database Quality Overview</h3>
      <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
        <div className="text-center">
          <div className="text-2xl font-bold text-blue-600">
            {qualityMetrics.totalMolecules}
          </div>
          <div className="text-sm text-gray-600">Total Molecules</div>
        </div>
        
        <div className="text-center">
          <div className="text-2xl font-bold text-green-600">
            {Math.round(qualityMetrics.propertiesCompletionRate)}%
          </div>
          <div className="text-sm text-gray-600">Properties Complete</div>
        </div>
        
        <div className="text-center">
          <div className="text-2xl font-bold text-purple-600">
            {Math.round(qualityMetrics.scoresCompletionRate)}%
          </div>
          <div className="text-sm text-gray-600">Scored Molecules</div>
        </div>
        
        <div className="text-center">
          <div className="text-2xl font-bold text-orange-600">
            {qualityMetrics.avgQualityScore}
          </div>
          <div className="text-sm text-gray-600">Avg Quality Score</div>
        </div>
      </div>
      
      <div className="mt-3 flex justify-center">
        {qualityMetrics.avgQualityScore >= 90 && (
          <Badge className="bg-green-100 text-green-800 border-green-200">
            🎯 Excellent Quality
          </Badge>
        )}
        {qualityMetrics.avgQualityScore >= 80 && qualityMetrics.avgQualityScore < 90 && (
          <Badge className="bg-blue-100 text-blue-800 border-blue-200">
            ✅ High Quality
          </Badge>
        )}
        {qualityMetrics.avgQualityScore < 80 && (
          <Badge className="bg-yellow-100 text-yellow-800 border-yellow-200">
            ⚡ Improving
          </Badge>
        )}
      </div>
    </Card>
  );
}

/**
 * Enhanced Molecules Page Component
 */
export default function EnhancedMoleculesPage() {
  return (
    <div className="container mx-auto py-8">
      {/* Header */}
      <div className="mb-8">
        <div className="flex items-center gap-3 mb-3">
          <h1 className="text-3xl font-bold">Enhanced Molecules</h1>
          <Badge className="bg-blue-100 text-blue-800 border-blue-200">
            🚀 Real-time Convex
          </Badge>
        </div>
        <p className="text-muted-foreground max-w-3xl">
          Experience our enhanced molecular database with real-time updates, advanced cryoprotectant scoring, 
          comprehensive molecular properties, and intelligent caching for lightning-fast performance.
        </p>
      </div>

      {/* Quality Summary */}
      <div className="mb-8">
        <QualitySummary />
      </div>

      {/* Features Overview */}
      <div className="mb-8">
        <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
          <Card className="p-4 text-center">
            <div className="text-2xl mb-2">🔄</div>
            <h4 className="font-semibold mb-1">Real-time Updates</h4>
            <p className="text-sm text-gray-600">Live data synchronization</p>
          </Card>
          
          <Card className="p-4 text-center">
            <div className="text-2xl mb-2">⚡</div>
            <h4 className="font-semibold mb-1">Smart Caching</h4>
            <p className="text-sm text-gray-600">5x faster performance</p>
          </Card>
          
          <Card className="p-4 text-center">
            <div className="text-2xl mb-2">🧪</div>
            <h4 className="font-semibold mb-1">Advanced Scoring</h4>
            <p className="text-sm text-gray-600">AI-powered cryoprotectant analysis</p>
          </Card>
          
          <Card className="p-4 text-center">
            <div className="text-2xl mb-2">📊</div>
            <h4 className="font-semibold mb-1">Rich Properties</h4>
            <p className="text-sm text-gray-600">Complete molecular descriptors</p>
          </Card>
        </div>
      </div>

      {/* Enhanced Molecules List */}
      <EnhancedMoleculesList />
    </div>
  );
}
