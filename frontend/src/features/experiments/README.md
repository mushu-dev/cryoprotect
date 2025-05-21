# Experimental Data Enhancement Feature

This directory contains the implementation of the Experimental Data Enhancement feature for CryoProtect.

## Feature Overview

The Experimental Data Enhancement feature provides improved visualization, analysis, and comparison capabilities for cryopreservation experiments. Key components include:

1. **Data Visualization**: Charts and graphs for visualizing experiment results
2. **Data Filtering**: Tools for filtering experiments by date, cryoprotectant, protocol, etc.
3. **Experiment Comparison**: Side-by-side comparison of multiple experiments
4. **Enhanced Detail View**: Improved experiment detail view with tabbed navigation
5. **Responsive Design**: Grid and list views that work on all device sizes

## Directory Structure

- `components/`: UI components for the feature
  - `ExperimentCard.js`: Card component for displaying experiment summaries
  - `ExperimentChart.js`: Chart component for visualizing experiment results
  - `ExperimentComparison.js`: Component for comparing multiple experiments
  - `ExperimentDetail.js`: Enhanced detail view for a single experiment
  - `ExperimentFilters.js`: Filtering interface for experiments
  - `ExperimentsList.js`: List view for experiments
  - `experiment-card.tsx`: TypeScript version of experiment card (newer implementation)
  - `experiment-results-chart.tsx`: TypeScript version of chart component
  - `experiments-list.tsx`: TypeScript version of list component

- `hooks/`: Custom React hooks for data management
  - `use-experiments.ts`: Hook for fetching and managing experiment data
  - `useExperimentData.js`: Hook for processing experiment data

- `services/`: Service layer for data fetching and processing
  - `experiment-service.ts`: Service for fetching experiment data

- `utils/`: Utility functions
  - `data-transformation.js`: Functions for transforming and analyzing experiment data

## Testing

We have comprehensive test coverage for this feature:

1. **User Flow Tests**: Basic navigation and interaction tests
   - `npm run test:user-flow`

2. **Data Integrity Tests**: Verify data consistency across pages
   - `npm run test:data-integrity`

3. **Experimental Feature Tests**: Specific tests for enhancement features
   - `npm run test:experimental-data`: Tests the UI components
   - `npm run test:data-analysis`: Tests the data processing functionality
   - `npm run test:experimental-features`: Runs both test suites above

4. **UI-based Testing**: Interactive test runner
   - `npm run test:experimental-features:ui`

5. **All Tests**: Run the complete test suite
   - `npm run test:experimental`

## Usage

The experimental data enhancement feature is accessible through the Experiments section of the application. Users can:

1. Browse experiments in grid or list view
2. Filter experiments by various criteria
3. View detailed experiment information with enhanced visualizations
4. Compare multiple experiments side-by-side
5. Export or share experiment results

## Development

When developing this feature, please ensure:

1. All components follow the established design system
2. Data fetching is done through the appropriate hooks
3. UI is responsive and works on all device sizes
4. All new code has corresponding test coverage
5. New functionality is documented in this README