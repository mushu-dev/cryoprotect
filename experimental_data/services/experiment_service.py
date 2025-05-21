#!/usr/bin/env python3
"""
Experiment Service for CryoProtect Enhanced Experimental Data System.

This module provides services for experiment management, including creating,
retrieving, updating, and analyzing experiments and their results.
"""

from typing import Dict, List, Any, Optional, Union, Tuple
import uuid
from datetime import datetime
import json

from ..models import (
    Experiment, 
    ExperimentResult, 
    Protocol, 
    TissueType,
    TimeSeries,
    TimeSeriesDataPoint,
    ValidationError,
    Uncertainty,
    Provenance
)

class ExperimentService:
    """Service for managing experiments and their results."""
    
    def __init__(self, db_adapter):
        """
        Initialize experiment service.
        
        Args:
            db_adapter: Database adapter for persistence
        """
        self.db_adapter = db_adapter
    
    async def create_experiment(self, experiment_data: Dict[str, Any]) -> Experiment:
        """
        Create a new experiment.
        
        Args:
            experiment_data: Data for the new experiment
            
        Returns:
            Created experiment
            
        Raises:
            ValidationError: If experiment data is invalid
        """
        # Create experiment model
        experiment = Experiment.from_dict(experiment_data)
        
        # Validate experiment
        experiment.validate()
        
        # Persist to database
        experiment_dict = experiment.to_dict()
        experiment_id = await self.db_adapter.create('experiments', experiment_dict)
        experiment.id = experiment_id
        
        return experiment
    
    async def get_experiment(self, experiment_id: str) -> Optional[Experiment]:
        """
        Get an experiment by ID.
        
        Args:
            experiment_id: ID of the experiment
            
        Returns:
            Experiment if found, None otherwise
        """
        experiment_data = await self.db_adapter.get('experiments', experiment_id)
        
        if not experiment_data:
            return None
        
        return Experiment.from_dict(experiment_data)
    
    async def list_experiments(
        self, 
        filters: Optional[Dict[str, Any]] = None,
        page: int = 1,
        page_size: int = 20,
        sort_by: str = 'created_at',
        sort_order: str = 'desc'
    ) -> Tuple[List[Experiment], int]:
        """
        List experiments with pagination and filtering.
        
        Args:
            filters: Filters to apply
            page: Page number (1-indexed)
            page_size: Page size
            sort_by: Field to sort by
            sort_order: Sort order ('asc' or 'desc')
            
        Returns:
            Tuple of (experiments, total_count)
        """
        # Apply pagination
        offset = (page - 1) * page_size
        limit = page_size
        
        # Fetch experiments
        experiments_data, total_count = await self.db_adapter.list(
            'experiments',
            filters=filters,
            offset=offset,
            limit=limit,
            sort_by=sort_by,
            sort_order=sort_order
        )
        
        # Convert to models
        experiments = [Experiment.from_dict(data) for data in experiments_data]
        
        return experiments, total_count
    
    async def update_experiment(self, experiment_id: str, data: Dict[str, Any]) -> Optional[Experiment]:
        """
        Update an experiment.
        
        Args:
            experiment_id: ID of the experiment to update
            data: Fields to update
            
        Returns:
            Updated experiment if found, None otherwise
            
        Raises:
            ValidationError: If updated data is invalid
        """
        # Get experiment
        experiment = await self.get_experiment(experiment_id)
        
        if not experiment:
            return None
        
        # Update fields
        experiment.update(data)
        
        # Persist to database
        experiment_dict = experiment.to_dict()
        await self.db_adapter.update('experiments', experiment_id, experiment_dict)
        
        return experiment
    
    async def delete_experiment(self, experiment_id: str) -> bool:
        """
        Delete an experiment.
        
        Args:
            experiment_id: ID of the experiment to delete
            
        Returns:
            True if deleted, False if not found
        """
        # Check if experiment exists
        experiment = await self.get_experiment(experiment_id)
        
        if not experiment:
            return False
        
        # Delete experiment
        result = await self.db_adapter.delete('experiments', experiment_id)
        
        return result
    
    async def create_experiment_result(self, result_data: Dict[str, Any]) -> ExperimentResult:
        """
        Create a new experiment result.
        
        Args:
            result_data: Data for the new result
            
        Returns:
            Created result
            
        Raises:
            ValidationError: If result data is invalid
        """
        # Create result model
        result = ExperimentResult.from_dict(result_data)
        
        # Validate result
        result.validate()
        
        # Persist to database
        result_dict = result.to_dict()
        result_id = await self.db_adapter.create('experiment_results', result_dict)
        result.id = result_id
        
        return result
    
    async def get_experiment_results(self, experiment_id: str) -> List[ExperimentResult]:
        """
        Get results for an experiment.
        
        Args:
            experiment_id: ID of the experiment
            
        Returns:
            List of experiment results
        """
        # Fetch results
        results_data = await self.db_adapter.list(
            'experiment_results',
            filters={'experiment_id': experiment_id},
            offset=0,
            limit=1000  # Arbitrary large limit
        )
        
        # Convert to models
        results = [ExperimentResult.from_dict(data) for data in results_data[0]]
        
        return results
    
    async def create_time_series(self, time_series_data: Dict[str, Any]) -> TimeSeries:
        """
        Create a new time series for an experiment.
        
        Args:
            time_series_data: Data for the new time series
            
        Returns:
            Created time series
            
        Raises:
            ValidationError: If time series data is invalid
        """
        # Create time series model
        time_series = TimeSeries.from_dict(time_series_data)
        
        # Validate time series
        time_series.validate()
        
        # Persist to database
        time_series_dict = time_series.to_dict()
        time_series_id = await self.db_adapter.create('time_series', time_series_dict)
        time_series.id = time_series_id
        
        return time_series
    
    async def add_time_series_data(
        self, 
        time_series_id: str, 
        data_points: List[Dict[str, Any]]
    ) -> List[TimeSeriesDataPoint]:
        """
        Add data points to a time series.
        
        Args:
            time_series_id: ID of the time series
            data_points: List of data points to add
            
        Returns:
            List of created data points
            
        Raises:
            ValidationError: If data points are invalid
        """
        # Create data point models
        points = []
        for point_data in data_points:
            # Ensure time_series_id is set
            point_data['time_series_id'] = time_series_id
            
            point = TimeSeriesDataPoint.from_dict(point_data)
            point.validate()
            points.append(point)
        
        # Persist to database
        created_points = []
        for point in points:
            point_dict = point.to_dict()
            point_id = await self.db_adapter.create('time_series_data', point_dict)
            point.id = point_id
            created_points.append(point)
        
        return created_points
    
    async def get_time_series_data(
        self, 
        time_series_id: str,
        start_time: Optional[datetime] = None,
        end_time: Optional[datetime] = None,
        limit: int = 1000
    ) -> List[TimeSeriesDataPoint]:
        """
        Get data points for a time series.
        
        Args:
            time_series_id: ID of the time series
            start_time: Start time for filtering
            end_time: End time for filtering
            limit: Maximum number of points to return
            
        Returns:
            List of time series data points
        """
        # Build filters
        filters = {'time_series_id': time_series_id}
        
        if start_time:
            filters['timestamp__gte'] = start_time.isoformat()
        
        if end_time:
            filters['timestamp__lte'] = end_time.isoformat()
        
        # Fetch data points
        points_data = await self.db_adapter.list(
            'time_series_data',
            filters=filters,
            offset=0,
            limit=limit,
            sort_by='timestamp',
            sort_order='asc'
        )
        
        # Convert to models
        points = [TimeSeriesDataPoint.from_dict(data) for data in points_data[0]]
        
        return points
    
    async def analyze_experiment_results(
        self, 
        experiment_id: str,
        analysis_type: str = 'basic'
    ) -> Dict[str, Any]:
        """
        Analyze results for an experiment.
        
        Args:
            experiment_id: ID of the experiment
            analysis_type: Type of analysis to perform
            
        Returns:
            Analysis results
        """
        # Get experiment results
        results = await self.get_experiment_results(experiment_id)
        
        if not results:
            return {
                'experiment_id': experiment_id,
                'analysis_type': analysis_type,
                'timestamp': datetime.now().isoformat(),
                'summary': {
                    'count': 0,
                    'message': 'No results found for this experiment'
                }
            }
        
        # Perform analysis
        if analysis_type == 'basic':
            return self._perform_basic_analysis(experiment_id, results)
        elif analysis_type == 'detailed':
            return self._perform_detailed_analysis(experiment_id, results)
        else:
            raise ValueError(f"Unsupported analysis type: {analysis_type}")
    
    def _perform_basic_analysis(
        self, 
        experiment_id: str, 
        results: List[ExperimentResult]
    ) -> Dict[str, Any]:
        """
        Perform basic analysis on experiment results.
        
        Args:
            experiment_id: ID of the experiment
            results: Experiment results
            
        Returns:
            Analysis results
        """
        # Extract numeric metrics
        viability_values = []
        recovery_values = []
        functionality_values = []
        
        for result in results:
            if result.viability_percentage is not None:
                viability_values.append(result.viability_percentage)
            
            if result.recovery_rate is not None:
                recovery_values.append(result.recovery_rate)
            
            if result.functionality_score is not None:
                functionality_values.append(result.functionality_score)
        
        # Calculate statistics
        viability_stats = self._calculate_statistics(viability_values)
        recovery_stats = self._calculate_statistics(recovery_values)
        functionality_stats = self._calculate_statistics(functionality_values)
        
        # Group by molecule/mixture
        molecule_groups = {}
        mixture_groups = {}
        
        for result in results:
            if result.molecule_id:
                if result.molecule_id not in molecule_groups:
                    molecule_groups[result.molecule_id] = []
                molecule_groups[result.molecule_id].append(result)
            
            if result.mixture_id:
                if result.mixture_id not in mixture_groups:
                    mixture_groups[result.mixture_id] = []
                mixture_groups[result.mixture_id].append(result)
        
        # Summarize groups
        molecule_summaries = {}
        for molecule_id, molecule_results in molecule_groups.items():
            molecule_summaries[molecule_id] = {
                'count': len(molecule_results),
                'viability': self._calculate_statistics([r.viability_percentage for r in molecule_results if r.viability_percentage is not None]),
                'recovery': self._calculate_statistics([r.recovery_rate for r in molecule_results if r.recovery_rate is not None]),
                'functionality': self._calculate_statistics([r.functionality_score for r in molecule_results if r.functionality_score is not None])
            }
        
        mixture_summaries = {}
        for mixture_id, mixture_results in mixture_groups.items():
            mixture_summaries[mixture_id] = {
                'count': len(mixture_results),
                'viability': self._calculate_statistics([r.viability_percentage for r in mixture_results if r.viability_percentage is not None]),
                'recovery': self._calculate_statistics([r.recovery_rate for r in mixture_results if r.recovery_rate is not None]),
                'functionality': self._calculate_statistics([r.functionality_score for r in mixture_results if r.functionality_score is not None])
            }
        
        # Compile analysis results
        analysis = {
            'experiment_id': experiment_id,
            'analysis_type': 'basic',
            'timestamp': datetime.now().isoformat(),
            'summary': {
                'count': len(results),
                'viability': viability_stats,
                'recovery': recovery_stats,
                'functionality': functionality_stats
            },
            'by_molecule': molecule_summaries,
            'by_mixture': mixture_summaries
        }
        
        return analysis
    
    def _perform_detailed_analysis(
        self, 
        experiment_id: str, 
        results: List[ExperimentResult]
    ) -> Dict[str, Any]:
        """
        Perform detailed analysis on experiment results.
        
        Args:
            experiment_id: ID of the experiment
            results: Experiment results
            
        Returns:
            Analysis results
        """
        # Start with basic analysis
        analysis = self._perform_basic_analysis(experiment_id, results)
        
        # Update analysis type
        analysis['analysis_type'] = 'detailed'
        
        # Add correlation analysis
        correlations = self._calculate_correlations(results)
        analysis['correlations'] = correlations
        
        # Add uncertainty analysis
        uncertainty_analysis = self._analyze_uncertainty(results)
        analysis['uncertainty'] = uncertainty_analysis
        
        # Add additional metrics
        analysis['metrics'] = {
            'success_rate': len([r for r in results if r.viability_percentage and r.viability_percentage > 50]) / len(results) if results else 0,
            'high_viability_count': len([r for r in results if r.viability_percentage and r.viability_percentage > 80]),
            'low_viability_count': len([r for r in results if r.viability_percentage and r.viability_percentage < 20])
        }
        
        # Add recommendations
        analysis['recommendations'] = self._generate_recommendations(results)
        
        return analysis
    
    def _calculate_statistics(self, values: List[float]) -> Dict[str, Any]:
        """
        Calculate descriptive statistics for a list of values.
        
        Args:
            values: List of numeric values
            
        Returns:
            Dictionary of statistics
        """
        if not values:
            return {
                'count': 0,
                'min': None,
                'max': None,
                'mean': None,
                'median': None,
                'std_dev': None
            }
        
        import numpy as np
        
        return {
            'count': len(values),
            'min': float(np.min(values)),
            'max': float(np.max(values)),
            'mean': float(np.mean(values)),
            'median': float(np.median(values)),
            'std_dev': float(np.std(values, ddof=1)) if len(values) > 1 else 0
        }
    
    def _calculate_correlations(self, results: List[ExperimentResult]) -> Dict[str, Any]:
        """
        Calculate correlations between metrics.
        
        Args:
            results: Experiment results
            
        Returns:
            Dictionary of correlation coefficients
        """
        import numpy as np
        from scipy.stats import pearsonr
        
        # Extract paired values for correlation
        viability_vs_recovery = []
        viability_vs_functionality = []
        recovery_vs_functionality = []
        
        for result in results:
            if result.viability_percentage is not None and result.recovery_rate is not None:
                viability_vs_recovery.append((result.viability_percentage, result.recovery_rate))
            
            if result.viability_percentage is not None and result.functionality_score is not None:
                viability_vs_functionality.append((result.viability_percentage, result.functionality_score))
            
            if result.recovery_rate is not None and result.functionality_score is not None:
                recovery_vs_functionality.append((result.recovery_rate, result.functionality_score))
        
        # Calculate correlations
        correlations = {}
        
        if len(viability_vs_recovery) > 1:
            viability = [x[0] for x in viability_vs_recovery]
            recovery = [x[1] for x in viability_vs_recovery]
            corr, p_value = pearsonr(viability, recovery)
            correlations['viability_vs_recovery'] = {
                'coefficient': float(corr),
                'p_value': float(p_value),
                'sample_size': len(viability_vs_recovery)
            }
        
        if len(viability_vs_functionality) > 1:
            viability = [x[0] for x in viability_vs_functionality]
            functionality = [x[1] for x in viability_vs_functionality]
            corr, p_value = pearsonr(viability, functionality)
            correlations['viability_vs_functionality'] = {
                'coefficient': float(corr),
                'p_value': float(p_value),
                'sample_size': len(viability_vs_functionality)
            }
        
        if len(recovery_vs_functionality) > 1:
            recovery = [x[0] for x in recovery_vs_functionality]
            functionality = [x[1] for x in recovery_vs_functionality]
            corr, p_value = pearsonr(recovery, functionality)
            correlations['recovery_vs_functionality'] = {
                'coefficient': float(corr),
                'p_value': float(p_value),
                'sample_size': len(recovery_vs_functionality)
            }
        
        return correlations
    
    def _analyze_uncertainty(self, results: List[ExperimentResult]) -> Dict[str, Any]:
        """
        Analyze uncertainty in results.
        
        Args:
            results: Experiment results
            
        Returns:
            Uncertainty analysis
        """
        # Extract uncertainty values
        viability_uncertainty = []
        recovery_uncertainty = []
        functionality_uncertainty = []
        
        for result in results:
            if 'viability_percentage' in result.uncertainty:
                u = result.uncertainty['viability_percentage']
                if isinstance(u, Uncertainty):
                    viability_uncertainty.append(u.value)
            
            if 'recovery_rate' in result.uncertainty:
                u = result.uncertainty['recovery_rate']
                if isinstance(u, Uncertainty):
                    recovery_uncertainty.append(u.value)
            
            if 'functionality_score' in result.uncertainty:
                u = result.uncertainty['functionality_score']
                if isinstance(u, Uncertainty):
                    functionality_uncertainty.append(u.value)
        
        # Calculate statistics
        viability_uncertainty_stats = self._calculate_statistics(viability_uncertainty)
        recovery_uncertainty_stats = self._calculate_statistics(recovery_uncertainty)
        functionality_uncertainty_stats = self._calculate_statistics(functionality_uncertainty)
        
        # Analyze uncertainty reporting
        uncertainty_completeness = {
            'viability': len(viability_uncertainty) / len([r for r in results if r.viability_percentage is not None]) if results else 0,
            'recovery': len(recovery_uncertainty) / len([r for r in results if r.recovery_rate is not None]) if results else 0,
            'functionality': len(functionality_uncertainty) / len([r for r in results if r.functionality_score is not None]) if results else 0
        }
        
        return {
            'viability': viability_uncertainty_stats,
            'recovery': recovery_uncertainty_stats,
            'functionality': functionality_uncertainty_stats,
            'completeness': uncertainty_completeness
        }
    
    def _generate_recommendations(self, results: List[ExperimentResult]) -> List[Dict[str, Any]]:
        """
        Generate recommendations based on analysis.
        
        Args:
            results: Experiment results
            
        Returns:
            List of recommendations
        """
        recommendations = []
        
        # Check for small sample size
        if len(results) < 3:
            recommendations.append({
                'type': 'sample_size',
                'severity': 'high',
                'message': 'Sample size is too small for reliable conclusions',
                'suggestion': 'Increase the number of experimental replicates'
            })
        
        # Check for missing uncertainty values
        viability_with_uncertainty = [r for r in results if r.viability_percentage is not None and 'viability_percentage' in r.uncertainty]
        if viability_with_uncertainty and len(viability_with_uncertainty) / len([r for r in results if r.viability_percentage is not None]) < 0.5:
            recommendations.append({
                'type': 'uncertainty',
                'severity': 'medium',
                'message': 'Many measurements lack uncertainty values',
                'suggestion': 'Report uncertainty for all measurements to improve reproducibility'
            })
        
        # Check for high variability
        viability_values = [r.viability_percentage for r in results if r.viability_percentage is not None]
        if viability_values:
            import numpy as np
            viability_std = np.std(viability_values, ddof=1)
            viability_mean = np.mean(viability_values)
            
            if viability_std / viability_mean > 0.25:  # Coefficient of variation > 25%
                recommendations.append({
                    'type': 'variability',
                    'severity': 'medium',
                    'message': 'High variability in viability measurements',
                    'suggestion': 'Review experimental conditions for sources of variability'
                })
        
        # Check for outliers
        if viability_values and len(viability_values) >= 5:
            import numpy as np
            from scipy import stats
            
            z_scores = np.abs(stats.zscore(viability_values))
            outliers = np.where(z_scores > 3)[0]
            
            if len(outliers) > 0:
                recommendations.append({
                    'type': 'outliers',
                    'severity': 'medium',
                    'message': f'Detected {len(outliers)} potential outliers in viability measurements',
                    'suggestion': 'Review outlier data points for measurement errors'
                })
        
        return recommendations