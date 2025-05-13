'use client'

import { Card, CardContent } from "@/components/ui/card"
import { formatNumber } from "@/lib/utils"
import { GitMerge, Beaker, LineChart } from "lucide-react"
import { Flask } from '@/components/ui/icons'

// This would typically come from an API endpoint
const mockStats = {
  molecules: {
    total: 1243,
    consolidated: 187,
    unique: 1056
  },
  mixtures: {
    total: 89,
    public: 67,
    private: 22
  },
  experiments: {
    total: 54,
    successful: 42,
    failed: 12
  },
  predictions: {
    total: 276,
    high_confidence: 168,
    medium_confidence: 87,
    low_confidence: 21
  }
}

export function DashboardStats() {
  return (
    <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-4">
      <Card>
        <CardContent className="p-6">
          <div className="flex items-center space-x-4">
            <div className="p-2 bg-blue-100 rounded-full dark:bg-blue-900">
              <Flask className="h-6 w-6 text-primary" />
            </div>
            <div>
              <p className="text-sm font-medium text-muted-foreground">Molecules</p>
              <h3 className="text-2xl font-bold">{formatNumber(mockStats.molecules.total)}</h3>
              <p className="text-xs text-muted-foreground mt-1">
                {formatNumber(mockStats.molecules.unique)} unique
              </p>
            </div>
          </div>
        </CardContent>
      </Card>
      
      <Card>
        <CardContent className="p-6">
          <div className="flex items-center space-x-4">
            <div className="p-2 bg-green-100 rounded-full dark:bg-green-900">
              <GitMerge className="h-6 w-6 text-secondary" />
            </div>
            <div>
              <p className="text-sm font-medium text-muted-foreground">Mixtures</p>
              <h3 className="text-2xl font-bold">{formatNumber(mockStats.mixtures.total)}</h3>
              <p className="text-xs text-muted-foreground mt-1">
                {formatNumber(mockStats.mixtures.public)} public
              </p>
            </div>
          </div>
        </CardContent>
      </Card>
      
      <Card>
        <CardContent className="p-6">
          <div className="flex items-center space-x-4">
            <div className="p-2 bg-yellow-100 rounded-full dark:bg-yellow-900">
              <Beaker className="h-6 w-6 text-yellow-600 dark:text-yellow-400" />
            </div>
            <div>
              <p className="text-sm font-medium text-muted-foreground">Experiments</p>
              <h3 className="text-2xl font-bold">{formatNumber(mockStats.experiments.total)}</h3>
              <p className="text-xs text-muted-foreground mt-1">
                {formatNumber(mockStats.experiments.successful)} successful
              </p>
            </div>
          </div>
        </CardContent>
      </Card>
      
      <Card>
        <CardContent className="p-6">
          <div className="flex items-center space-x-4">
            <div className="p-2 bg-purple-100 rounded-full dark:bg-purple-900">
              <LineChart className="h-6 w-6 text-purple-600 dark:text-purple-400" />
            </div>
            <div>
              <p className="text-sm font-medium text-muted-foreground">Predictions</p>
              <h3 className="text-2xl font-bold">{formatNumber(mockStats.predictions.total)}</h3>
              <p className="text-xs text-muted-foreground mt-1">
                {formatNumber(mockStats.predictions.high_confidence)} high confidence
              </p>
            </div>
          </div>
        </CardContent>
      </Card>
    </div>
  )
}