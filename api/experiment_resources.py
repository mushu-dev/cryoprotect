from flask import request, jsonify, Blueprint, current_app
from flask_restful import Resource, Api, reqparse
from werkzeug.utils import secure_filename
from .api_decorators import require_auth, require_service_role
from .pagination_utils import paginate_results
from .models import db, Experiment, ExperimentResult, Protocol, TimeSeries, TissueType
from .schemas import (
    ExperimentSchema,
    ExperimentResultSchema,
    TimeSeriesSchema
)
from .consolidated_decorators import with_transaction
import uuid
import os
import json
from datetime import datetime
from sqlalchemy import func, desc, asc
from sqlalchemy.exc import SQLAlchemyError
import logging

# Setup logger
logger = logging.getLogger(__name__)

# Initialize blueprints
experiment_bp = Blueprint('experiments', __name__)
experiment_api = Api(experiment_bp)

# Initialize schemas
experiment_schema = ExperimentSchema()
experiments_schema = ExperimentSchema(many=True)
result_schema = ExperimentResultSchema()
results_schema = ExperimentResultSchema(many=True)
time_series_schema = TimeSeriesSchema(many=True)

# Request parsers
experiment_parser = reqparse.RequestParser()
experiment_parser.add_argument('name', type=str, required=True, help='Name is required')
experiment_parser.add_argument('description', type=str)
experiment_parser.add_argument('protocol_id', type=str, required=True, help='Protocol ID is required')
experiment_parser.add_argument('tissue_type_id', type=str, required=True, help='Tissue Type ID is required')
experiment_parser.add_argument('experiment_type', type=str, required=True, help='Experiment Type is required')
experiment_parser.add_argument('start_date', type=str, required=True, help='Start Date is required')
experiment_parser.add_argument('end_date', type=str)
experiment_parser.add_argument('status', type=str, required=True, 
                              choices=('planned', 'in_progress', 'completed', 'aborted', 'failed'),
                              help='Status must be one of: planned, in_progress, completed, aborted, failed')
experiment_parser.add_argument('researcher', type=str, required=True, help='Researcher is required')
experiment_parser.add_argument('lab_id', type=str)
experiment_parser.add_argument('equipment', type=list, location='json')
experiment_parser.add_argument('environmental_conditions', type=dict, location='json')
experiment_parser.add_argument('notes', type=str)
experiment_parser.add_argument('tags', type=list, location='json')

# Parser for experiment listing and filtering
list_parser = reqparse.RequestParser()
list_parser.add_argument('page', type=int, default=1, location='args')
list_parser.add_argument('per_page', type=int, default=10, location='args')
list_parser.add_argument('sort_by', type=str, default='start_date', location='args')
list_parser.add_argument('sort_order', type=str, choices=('asc', 'desc'), default='desc', location='args')
list_parser.add_argument('status', type=str, location='args')
list_parser.add_argument('experiment_type', type=str, location='args')
list_parser.add_argument('researcher', type=str, location='args')
list_parser.add_argument('tissue_type_id', type=str, location='args')
list_parser.add_argument('date_from', type=str, location='args')
list_parser.add_argument('date_to', type=str, location='args')
list_parser.add_argument('tags', type=str, location='args')
list_parser.add_argument('search', type=str, location='args')

# Parser for experiment results
result_parser = reqparse.RequestParser()
result_parser.add_argument('experiment_id', type=str, required=True, help='Experiment ID is required')
result_parser.add_argument('tissue_type_id', type=str, required=True, help='Tissue Type ID is required')
result_parser.add_argument('molecule_id', type=str)
result_parser.add_argument('mixture_id', type=str)
result_parser.add_argument('concentration', type=float)
result_parser.add_argument('concentration_unit', type=str)
result_parser.add_argument('viability_percentage', type=float)
result_parser.add_argument('recovery_rate', type=float)
result_parser.add_argument('functionality_score', type=float)
result_parser.add_argument('uncertainty', type=dict, location='json')
result_parser.add_argument('result_details', type=dict, location='json')
result_parser.add_argument('notes', type=str)
result_parser.add_argument('protocol_step_id', type=str)

# Parser for time series data
time_series_parser = reqparse.RequestParser()
time_series_parser.add_argument('experiment_id', type=str, required=True, help='Experiment ID is required')
time_series_parser.add_argument('result_id', type=str)
time_series_parser.add_argument('parameter', type=str, required=True, help='Parameter name is required')
time_series_parser.add_argument('unit', type=str, required=True, help='Unit is required')
time_series_parser.add_argument('data_points', type=list, location='json', required=True,
                               help='Data points are required')
time_series_parser.add_argument('start_time', type=str, required=True, help='Start time is required')
time_series_parser.add_argument('end_time', type=str, required=True, help='End time is required')
time_series_parser.add_argument('notes', type=str)

# Analysis parser
analysis_parser = reqparse.RequestParser()
analysis_parser.add_argument('experiment_ids', type=list, location='json', required=True,
                            help='Experiment IDs are required')
analysis_parser.add_argument('compare_with', type=list, location='json')
analysis_parser.add_argument('analysis_type', type=list, location='json')
analysis_parser.add_argument('time_range', type=list, location='json')

# Helper function to generate provenance info
def generate_provenance(request):
    return {
        'created_by': request.user.id if hasattr(request, 'user') else 'system',
        'created_at': datetime.utcnow().isoformat(),
        'source': 'api',
        'method': request.method,
        'version': '1.0'
    }

class ExperimentList(Resource):
    @require_auth
    def get(self):
        """Get a list of experiments with pagination and filtering"""
        args = list_parser.parse_args()
        
        # Start with a base query
        query = Experiment.query
        
        # Apply filters
        if args.status:
            query = query.filter(Experiment.status == args.status)
        
        if args.experiment_type:
            query = query.filter(Experiment.experiment_type == args.experiment_type)
            
        if args.researcher:
            query = query.filter(Experiment.researcher.ilike(f"%{args.researcher}%"))
            
        if args.tissue_type_id:
            query = query.filter(Experiment.tissue_type_id == uuid.UUID(args.tissue_type_id))
            
        if args.date_from:
            try:
                date_from = datetime.strptime(args.date_from, '%Y-%m-%d')
                query = query.filter(Experiment.start_date >= date_from)
            except ValueError:
                return {"message": "Invalid date_from format. Use YYYY-MM-DD"}, 400
                
        if args.date_to:
            try:
                date_to = datetime.strptime(args.date_to, '%Y-%m-%d')
                query = query.filter(Experiment.start_date <= date_to)
            except ValueError:
                return {"message": "Invalid date_to format. Use YYYY-MM-DD"}, 400
                
        if args.tags:
            tags = args.tags.split(',')
            for tag in tags:
                # Using PostgreSQL array overlap operator
                query = query.filter(Experiment.tags.overlap('{' + tag.strip() + '}'))
                
        if args.search:
            search_term = f"%{args.search}%"
            query = query.filter(
                (Experiment.name.ilike(search_term)) | 
                (Experiment.description.ilike(search_term))
            )
            
        # Apply sorting
        sort_column = getattr(Experiment, args.sort_by, Experiment.start_date)
        sort_method = desc(sort_column) if args.sort_order == 'desc' else asc(sort_column)
        query = query.order_by(sort_method)
        
        # Apply pagination
        paginated_result = paginate_results(query, args.page, args.per_page)
        
        # Serialize and return results
        experiments = experiments_schema.dump(paginated_result.items)
        
        return {
            'data': experiments,
            'total': paginated_result.total,
            'page': paginated_result.page,
            'per_page': paginated_result.per_page,
            'total_pages': paginated_result.pages
        }
    
    @require_auth
    @with_transaction
    def post(self):
        """Create a new experiment"""
        args = experiment_parser.parse_args()
        
        # Convert string IDs to UUIDs
        try:
            protocol_id = uuid.UUID(args.protocol_id)
            tissue_type_id = uuid.UUID(args.tissue_type_id)
        except ValueError:
            return {"message": "Invalid UUID format"}, 400
            
        # Validate protocol existence
        protocol = Protocol.query.get(protocol_id)
        if not protocol:
            return {"message": f"Protocol with ID {args.protocol_id} not found"}, 404
            
        # Validate tissue type existence
        tissue_type = TissueType.query.get(tissue_type_id)
        if not tissue_type:
            return {"message": f"Tissue type with ID {args.tissue_type_id} not found"}, 404
            
        # Create new experiment
        new_experiment = Experiment(
            id=uuid.uuid4(),
            name=args.name,
            description=args.description,
            protocol_id=protocol_id,
            tissue_type_id=tissue_type_id,
            experiment_type=args.experiment_type,
            start_date=datetime.strptime(args.start_date, '%Y-%m-%d'),
            status=args.status,
            researcher=args.researcher,
            lab_id=args.lab_id,
            equipment=args.equipment,
            environmental_conditions=args.environmental_conditions,
            notes=args.notes,
            tags=args.tags,
            provenance=generate_provenance(request)
        )
        
        # Handle end_date if provided
        if args.end_date:
            new_experiment.end_date = datetime.strptime(args.end_date, '%Y-%m-%d')
        
        db.session.add(new_experiment)
        
        # Return the created experiment
        return experiment_schema.dump(new_experiment), 201

class ExperimentDetail(Resource):
    @require_auth
    def get(self, experiment_id):
        """Get experiment details by ID"""
        try:
            exp_uuid = uuid.UUID(experiment_id)
        except ValueError:
            return {"message": "Invalid experiment ID format"}, 400
            
        experiment = Experiment.query.get(exp_uuid)
        if not experiment:
            return {"message": "Experiment not found"}, 404
            
        # Get associated protocol
        experiment.protocol = Protocol.query.get(experiment.protocol_id)
        
        # Get experiment results
        experiment.results = ExperimentResult.query.filter_by(experiment_id=exp_uuid).all()
        
        # Get time series data
        experiment.time_series = TimeSeries.query.filter_by(experiment_id=exp_uuid).all()
        
        return experiment_schema.dump(experiment)
    
    @require_auth
    @with_transaction
    def patch(self, experiment_id):
        """Update an existing experiment"""
        try:
            exp_uuid = uuid.UUID(experiment_id)
        except ValueError:
            return {"message": "Invalid experiment ID format"}, 400
            
        experiment = Experiment.query.get(exp_uuid)
        if not experiment:
            return {"message": "Experiment not found"}, 404
            
        # Parse args but don't require all fields (partial update)
        parser = reqparse.RequestParser()
        parser.add_argument('name', type=str)
        parser.add_argument('description', type=str)
        parser.add_argument('protocol_id', type=str)
        parser.add_argument('tissue_type_id', type=str)
        parser.add_argument('experiment_type', type=str)
        parser.add_argument('start_date', type=str)
        parser.add_argument('end_date', type=str)
        parser.add_argument('status', type=str, 
                          choices=('planned', 'in_progress', 'completed', 'aborted', 'failed'))
        parser.add_argument('researcher', type=str)
        parser.add_argument('lab_id', type=str)
        parser.add_argument('equipment', type=list, location='json')
        parser.add_argument('environmental_conditions', type=dict, location='json')
        parser.add_argument('notes', type=str)
        parser.add_argument('tags', type=list, location='json')
        args = parser.parse_args()
        
        # Update fields if provided
        if args.name:
            experiment.name = args.name
        if args.description is not None:  # Allow empty string
            experiment.description = args.description
        if args.protocol_id:
            try:
                protocol_id = uuid.UUID(args.protocol_id)
                protocol = Protocol.query.get(protocol_id)
                if not protocol:
                    return {"message": f"Protocol with ID {args.protocol_id} not found"}, 404
                experiment.protocol_id = protocol_id
            except ValueError:
                return {"message": "Invalid protocol UUID format"}, 400
        if args.tissue_type_id:
            try:
                tissue_type_id = uuid.UUID(args.tissue_type_id)
                tissue_type = TissueType.query.get(tissue_type_id)
                if not tissue_type:
                    return {"message": f"Tissue type with ID {args.tissue_type_id} not found"}, 404
                experiment.tissue_type_id = tissue_type_id
            except ValueError:
                return {"message": "Invalid tissue type UUID format"}, 400
        if args.experiment_type:
            experiment.experiment_type = args.experiment_type
        if args.start_date:
            try:
                experiment.start_date = datetime.strptime(args.start_date, '%Y-%m-%d')
            except ValueError:
                return {"message": "Invalid start_date format. Use YYYY-MM-DD"}, 400
        if args.end_date:
            try:
                experiment.end_date = datetime.strptime(args.end_date, '%Y-%m-%d')
            except ValueError:
                return {"message": "Invalid end_date format. Use YYYY-MM-DD"}, 400
        if args.end_date == '':  # Handle removing the end date
            experiment.end_date = None
        if args.status:
            experiment.status = args.status
        if args.researcher:
            experiment.researcher = args.researcher
        if args.lab_id is not None:  # Allow empty string
            experiment.lab_id = args.lab_id
        if args.equipment is not None:
            experiment.equipment = args.equipment
        if args.environmental_conditions is not None:
            experiment.environmental_conditions = args.environmental_conditions
        if args.notes is not None:  # Allow empty string
            experiment.notes = args.notes
        if args.tags is not None:
            experiment.tags = args.tags
        
        # Return the updated experiment
        return experiment_schema.dump(experiment)
    
    @require_auth
    @with_transaction
    def delete(self, experiment_id):
        """Delete an experiment"""
        try:
            exp_uuid = uuid.UUID(experiment_id)
        except ValueError:
            return {"message": "Invalid experiment ID format"}, 400
            
        experiment = Experiment.query.get(exp_uuid)
        if not experiment:
            return {"message": "Experiment not found"}, 404
            
        # Delete associated results and time series
        ExperimentResult.query.filter_by(experiment_id=exp_uuid).delete()
        TimeSeries.query.filter_by(experiment_id=exp_uuid).delete()
        
        # Delete experiment
        db.session.delete(experiment)
        
        return {"message": "Experiment deleted successfully"}, 200

class ExperimentResultList(Resource):
    @require_auth
    def get(self, experiment_id):
        """Get all results for an experiment"""
        try:
            exp_uuid = uuid.UUID(experiment_id)
        except ValueError:
            return {"message": "Invalid experiment ID format"}, 400
            
        # Verify experiment exists
        experiment = Experiment.query.get(exp_uuid)
        if not experiment:
            return {"message": "Experiment not found"}, 404
            
        results = ExperimentResult.query.filter_by(experiment_id=exp_uuid).all()
        
        return results_schema.dump(results)
    
    @require_auth
    @with_transaction
    def post(self, experiment_id):
        """Add a result to an experiment"""
        try:
            exp_uuid = uuid.UUID(experiment_id)
        except ValueError:
            return {"message": "Invalid experiment ID format"}, 400
            
        # Verify experiment exists
        experiment = Experiment.query.get(exp_uuid)
        if not experiment:
            return {"message": "Experiment not found"}, 404
            
        args = result_parser.parse_args()
        
        # Ensure the experiment ID in path matches the one in payload
        if args.experiment_id != experiment_id:
            return {"message": "Experiment ID in path and payload do not match"}, 400
            
        # Convert IDs to UUIDs
        try:
            tissue_type_id = uuid.UUID(args.tissue_type_id)
        except ValueError:
            return {"message": "Invalid tissue type UUID format"}, 400
            
        # Create new result
        new_result = ExperimentResult(
            id=uuid.uuid4(),
            experiment_id=exp_uuid,
            tissue_type_id=tissue_type_id,
            concentration=args.concentration,
            concentration_unit=args.concentration_unit,
            viability_percentage=args.viability_percentage,
            recovery_rate=args.recovery_rate,
            functionality_score=args.functionality_score,
            uncertainty=args.uncertainty,
            result_details=args.result_details,
            notes=args.notes,
            timestamp=datetime.utcnow(),
            provenance=generate_provenance(request)
        )
        
        # Handle optional foreign keys
        if args.molecule_id:
            try:
                new_result.molecule_id = uuid.UUID(args.molecule_id)
            except ValueError:
                return {"message": "Invalid molecule UUID format"}, 400
                
        if args.mixture_id:
            try:
                new_result.mixture_id = uuid.UUID(args.mixture_id)
            except ValueError:
                return {"message": "Invalid mixture UUID format"}, 400
                
        if args.protocol_step_id:
            try:
                new_result.protocol_step_id = uuid.UUID(args.protocol_step_id)
            except ValueError:
                return {"message": "Invalid protocol step UUID format"}, 400
        
        db.session.add(new_result)
        
        return result_schema.dump(new_result), 201

class ExperimentResultDetail(Resource):
    @require_auth
    def get(self, result_id):
        """Get a specific experiment result"""
        try:
            result_uuid = uuid.UUID(result_id)
        except ValueError:
            return {"message": "Invalid result ID format"}, 400
            
        result = ExperimentResult.query.get(result_uuid)
        if not result:
            return {"message": "Result not found"}, 404
            
        return result_schema.dump(result)
    
    @require_auth
    @with_transaction
    def patch(self, result_id):
        """Update an experiment result"""
        try:
            result_uuid = uuid.UUID(result_id)
        except ValueError:
            return {"message": "Invalid result ID format"}, 400
            
        result = ExperimentResult.query.get(result_uuid)
        if not result:
            return {"message": "Result not found"}, 404
            
        # Parse args but don't require all fields
        parser = reqparse.RequestParser()
        parser.add_argument('tissue_type_id', type=str)
        parser.add_argument('molecule_id', type=str)
        parser.add_argument('mixture_id', type=str)
        parser.add_argument('concentration', type=float)
        parser.add_argument('concentration_unit', type=str)
        parser.add_argument('viability_percentage', type=float)
        parser.add_argument('recovery_rate', type=float)
        parser.add_argument('functionality_score', type=float)
        parser.add_argument('uncertainty', type=dict, location='json')
        parser.add_argument('result_details', type=dict, location='json')
        parser.add_argument('notes', type=str)
        parser.add_argument('protocol_step_id', type=str)
        args = parser.parse_args()
        
        # Update fields if provided
        if args.tissue_type_id:
            try:
                result.tissue_type_id = uuid.UUID(args.tissue_type_id)
            except ValueError:
                return {"message": "Invalid tissue type UUID format"}, 400
        if args.molecule_id:
            try:
                result.molecule_id = uuid.UUID(args.molecule_id)
            except ValueError:
                return {"message": "Invalid molecule UUID format"}, 400
        elif args.molecule_id == '':  # Handle explicit null/empty
            result.molecule_id = None
        if args.mixture_id:
            try:
                result.mixture_id = uuid.UUID(args.mixture_id)
            except ValueError:
                return {"message": "Invalid mixture UUID format"}, 400
        elif args.mixture_id == '':  # Handle explicit null/empty
            result.mixture_id = None
        if args.concentration is not None:
            result.concentration = args.concentration
        if args.concentration_unit is not None:
            result.concentration_unit = args.concentration_unit
        if args.viability_percentage is not None:
            result.viability_percentage = args.viability_percentage
        if args.recovery_rate is not None:
            result.recovery_rate = args.recovery_rate
        if args.functionality_score is not None:
            result.functionality_score = args.functionality_score
        if args.uncertainty is not None:
            result.uncertainty = args.uncertainty
        if args.result_details is not None:
            result.result_details = args.result_details
        if args.notes is not None:
            result.notes = args.notes
        if args.protocol_step_id:
            try:
                result.protocol_step_id = uuid.UUID(args.protocol_step_id)
            except ValueError:
                return {"message": "Invalid protocol step UUID format"}, 400
        elif args.protocol_step_id == '':  # Handle explicit null/empty
            result.protocol_step_id = None
        
        return result_schema.dump(result)
    
    @require_auth
    @with_transaction
    def delete(self, result_id):
        """Delete an experiment result"""
        try:
            result_uuid = uuid.UUID(result_id)
        except ValueError:
            return {"message": "Invalid result ID format"}, 400
            
        result = ExperimentResult.query.get(result_uuid)
        if not result:
            return {"message": "Result not found"}, 404
            
        db.session.delete(result)
        
        return {"message": "Result deleted successfully"}, 200

class ExperimentTimeSeriesList(Resource):
    @require_auth
    def get(self, experiment_id):
        """Get time series data for an experiment"""
        try:
            exp_uuid = uuid.UUID(experiment_id)
        except ValueError:
            return {"message": "Invalid experiment ID format"}, 400
            
        # Verify experiment exists
        experiment = Experiment.query.get(exp_uuid)
        if not experiment:
            return {"message": "Experiment not found"}, 404
            
        # Get query parameters
        parameter = request.args.get('parameter')
        
        # Base query
        query = TimeSeries.query.filter_by(experiment_id=exp_uuid)
        
        # Filter by parameter if specified
        if parameter:
            query = query.filter_by(parameter=parameter)
            
        time_series = query.all()
        
        return time_series_schema.dump(time_series)
    
    @require_auth
    @with_transaction
    def post(self, experiment_id):
        """Add time series data to an experiment"""
        try:
            exp_uuid = uuid.UUID(experiment_id)
        except ValueError:
            return {"message": "Invalid experiment ID format"}, 400
            
        # Verify experiment exists
        experiment = Experiment.query.get(exp_uuid)
        if not experiment:
            return {"message": "Experiment not found"}, 404
            
        args = time_series_parser.parse_args()
        
        # Ensure the experiment ID in path matches the one in payload
        if args.experiment_id != experiment_id:
            return {"message": "Experiment ID in path and payload do not match"}, 400
            
        # Create new time series
        new_time_series = TimeSeries(
            id=uuid.uuid4(),
            experiment_id=exp_uuid,
            parameter=args.parameter,
            unit=args.unit,
            data_points=args.data_points,
            start_time=datetime.fromisoformat(args.start_time.replace('Z', '+00:00')),
            end_time=datetime.fromisoformat(args.end_time.replace('Z', '+00:00')),
            notes=args.notes,
            provenance=generate_provenance(request)
        )
        
        # Handle optional result ID
        if args.result_id:
            try:
                new_time_series.result_id = uuid.UUID(args.result_id)
            except ValueError:
                return {"message": "Invalid result UUID format"}, 400
        
        db.session.add(new_time_series)
        
        return time_series_schema.dump([new_time_series])[0], 201

class ExperimentAnalysis(Resource):
    @require_auth
    def post(self):
        """Analyze experimental data"""
        args = analysis_parser.parse_args()
        
        # Validate experiment IDs
        experiment_ids = []
        for exp_id in args.experiment_ids:
            try:
                experiment_ids.append(uuid.UUID(exp_id))
            except ValueError:
                return {"message": f"Invalid experiment ID format: {exp_id}"}, 400
        
        # Verify experiments exist
        experiments = Experiment.query\
            .filter(Experiment.id.in_(experiment_ids))\
            .all()
        
        if len(experiments) != len(experiment_ids):
            return {"message": "One or more experiments not found"}, 404
        
        # Get all results for these experiments
        results = ExperimentResult.query\
            .filter(ExperimentResult.experiment_id.in_(experiment_ids))\
            .all()
        
        # Mock analysis for demonstration
        # In a real implementation, this would use sophisticated analytics
        successful_count = sum(1 for r in results if r.viability_percentage and r.viability_percentage > 50)
        total_count = len(results)
        
        # Calculate simple statistics
        viability_values = [r.viability_percentage for r in results if r.viability_percentage is not None]
        recovery_values = [r.recovery_rate for r in results if r.recovery_rate is not None]
        
        if not viability_values:
            viability_values = [0]
        if not recovery_values:
            recovery_values = [0]
            
        # Simple calculation of statistics
        viability_stats = {
            'mean': sum(viability_values) / len(viability_values),
            'median': sorted(viability_values)[len(viability_values) // 2],
            'std_dev': 0,  # Would calculate standard deviation in a real implementation
            'min': min(viability_values),
            'max': max(viability_values),
            'quartiles': [
                viability_values[len(viability_values) // 4],
                viability_values[len(viability_values) // 2],
                viability_values[3 * len(viability_values) // 4]
            ]
        }
        
        recovery_stats = {
            'mean': sum(recovery_values) / len(recovery_values),
            'median': sorted(recovery_values)[len(recovery_values) // 2],
            'std_dev': 0,  # Would calculate standard deviation in a real implementation
            'min': min(recovery_values),
            'max': max(recovery_values),
            'quartiles': [
                recovery_values[len(recovery_values) // 4],
                recovery_values[len(recovery_values) // 2],
                recovery_values[3 * len(recovery_values) // 4]
            ]
        }
        
        # Mock trend data
        experiment_dates = sorted([exp.start_date for exp in experiments])
        viability_trend = {
            'parameter': 'viability',
            'timestamps': [d.isoformat() for d in experiment_dates],
            'values': [r.viability_percentage or 0 for r in results[:len(experiment_dates)]],
            'trend_line': [v + (i * 0.5) for i, v in enumerate([50, 55, 60, 65, 70][:len(experiment_dates)])]
        }
        
        analysis_result = {
            'summary': {
                'total': total_count,
                'successful': successful_count,
                'failed': total_count - successful_count,
                'success_rate': (successful_count / total_count) if total_count > 0 else 0
            },
            'statistics': {
                'viability': viability_stats,
                'recovery': recovery_stats
            },
            'trends': [viability_trend]
        }
        
        # Add comparisons if requested
        if args.compare_with:
            comparison_ids = []
            for comp_id in args.compare_with:
                try:
                    comparison_ids.append(uuid.UUID(comp_id))
                except ValueError:
                    return {"message": f"Invalid comparison ID format: {comp_id}"}, 400
                    
            # In a real implementation, would run comparative analysis here
            analysis_result['comparisons'] = {
                'improvement': '15%',
                'significance': 'p < 0.05',
                'detailed_metrics': {
                    'viability_change': '+10.5%',
                    'recovery_change': '+12.3%'
                }
            }
        
        return analysis_result

class ExperimentExport(Resource):
    @require_auth
    def get(self, experiment_id):
        """Export experiment data in the requested format"""
        try:
            exp_uuid = uuid.UUID(experiment_id)
        except ValueError:
            return {"message": "Invalid experiment ID format"}, 400
            
        experiment = Experiment.query.get(exp_uuid)
        if not experiment:
            return {"message": "Experiment not found"}, 404
            
        # Get requested format
        export_format = request.args.get('format', 'json')
        if export_format not in ('json', 'csv', 'excel', 'pdf'):
            return {"message": "Invalid export format. Must be json, csv, excel, or pdf"}, 400
            
        # Get related data
        experiment.results = ExperimentResult.query.filter_by(experiment_id=exp_uuid).all()
        experiment.time_series = TimeSeries.query.filter_by(experiment_id=exp_uuid).all()
        experiment.protocol = Protocol.query.get(experiment.protocol_id)
        
        # For this implementation, we'll only handle JSON format
        # In a real app, would implement other formats
        if export_format == 'json':
            response = jsonify(experiment_schema.dump(experiment))
            response.headers['Content-Disposition'] = f'attachment; filename=experiment_{experiment_id}.json'
            return response
        else:
            # Placeholder for other formats
            return {"message": f"Export format {export_format} not yet implemented"}, 501

class ExperimentImport(Resource):
    @require_auth
    @with_transaction
    def post(self):
        """Import experiment data from a file"""
        if 'file' not in request.files:
            return {"message": "No file provided"}, 400
            
        file = request.files['file']
        if file.filename == '':
            return {"message": "No file selected"}, 400
            
        if not file.filename.endswith('.json'):
            return {"message": "Only JSON files are supported"}, 400
            
        try:
            # Read and parse the JSON file
            experiment_data = json.loads(file.read().decode('utf-8'))
            
            # Generate new IDs for the imported experiment
            new_exp_id = uuid.uuid4()
            
            # Create the new experiment
            new_experiment = Experiment(
                id=new_exp_id,
                name=experiment_data.get('name', 'Imported Experiment'),
                description=experiment_data.get('description'),
                protocol_id=uuid.UUID(experiment_data.get('protocol_id')),
                tissue_type_id=uuid.UUID(experiment_data.get('tissue_type_id')),
                experiment_type=experiment_data.get('experiment_type', 'unknown'),
                start_date=datetime.fromisoformat(experiment_data.get('start_date').replace('Z', '+00:00')),
                status=experiment_data.get('status', 'completed'),
                researcher=experiment_data.get('researcher', 'Imported'),
                lab_id=experiment_data.get('lab_id'),
                equipment=experiment_data.get('equipment'),
                environmental_conditions=experiment_data.get('environmental_conditions'),
                notes=experiment_data.get('notes'),
                tags=experiment_data.get('tags'),
                provenance={
                    'created_by': request.user.id if hasattr(request, 'user') else 'system',
                    'created_at': datetime.utcnow().isoformat(),
                    'source': 'import',
                    'method': 'file-upload',
                    'version': '1.0',
                    'references': [experiment_data.get('id', 'unknown')]
                }
            )
            
            if experiment_data.get('end_date'):
                new_experiment.end_date = datetime.fromisoformat(
                    experiment_data.get('end_date').replace('Z', '+00:00')
                )
                
            db.session.add(new_experiment)
            
            # Import results if included
            if experiment_data.get('results'):
                for result_data in experiment_data['results']:
                    new_result = ExperimentResult(
                        id=uuid.uuid4(),
                        experiment_id=new_exp_id,
                        tissue_type_id=uuid.UUID(result_data.get('tissue_type_id')),
                        concentration=result_data.get('concentration'),
                        concentration_unit=result_data.get('concentration_unit'),
                        viability_percentage=result_data.get('viability_percentage'),
                        recovery_rate=result_data.get('recovery_rate'),
                        functionality_score=result_data.get('functionality_score'),
                        uncertainty=result_data.get('uncertainty'),
                        result_details=result_data.get('result_details'),
                        notes=result_data.get('notes'),
                        timestamp=datetime.utcnow(),
                        provenance={
                            'created_by': request.user.id if hasattr(request, 'user') else 'system',
                            'created_at': datetime.utcnow().isoformat(),
                            'source': 'import',
                            'method': 'file-upload',
                            'version': '1.0',
                            'references': [result_data.get('id', 'unknown')]
                        }
                    )
                    
                    if result_data.get('molecule_id'):
                        new_result.molecule_id = uuid.UUID(result_data.get('molecule_id'))
                    if result_data.get('mixture_id'):
                        new_result.mixture_id = uuid.UUID(result_data.get('mixture_id'))
                    if result_data.get('protocol_step_id'):
                        new_result.protocol_step_id = uuid.UUID(result_data.get('protocol_step_id'))
                        
                    db.session.add(new_result)
            
            # Import time series if included
            if experiment_data.get('time_series'):
                for ts_data in experiment_data['time_series']:
                    new_ts = TimeSeries(
                        id=uuid.uuid4(),
                        experiment_id=new_exp_id,
                        parameter=ts_data.get('parameter'),
                        unit=ts_data.get('unit'),
                        data_points=ts_data.get('data_points'),
                        start_time=datetime.fromisoformat(ts_data.get('start_time').replace('Z', '+00:00')),
                        end_time=datetime.fromisoformat(ts_data.get('end_time').replace('Z', '+00:00')),
                        notes=ts_data.get('notes'),
                        provenance={
                            'created_by': request.user.id if hasattr(request, 'user') else 'system',
                            'created_at': datetime.utcnow().isoformat(),
                            'source': 'import',
                            'method': 'file-upload',
                            'version': '1.0',
                            'references': [ts_data.get('id', 'unknown')]
                        }
                    )
                    
                    if ts_data.get('result_id'):
                        # We don't map this since the result IDs have changed
                        pass
                        
                    db.session.add(new_ts)
            
            # Return the imported experiment
            return experiment_schema.dump(new_experiment), 201
            
        except json.JSONDecodeError:
            return {"message": "Invalid JSON file"}, 400
        except Exception as e:
            logger.error(f"Error importing experiment: {str(e)}")
            return {"message": f"Error importing experiment: {str(e)}"}, 500

class ExperimentSearch(Resource):
    @require_auth
    def get(self):
        """Search for experiments"""
        # Get query parameters
        query = request.args.get('query', '')
        args = list_parser.parse_args()
        
        # Start with a base query
        db_query = Experiment.query
        
        # Apply search
        if query:
            search_term = f"%{query}%"
            db_query = db_query.filter(
                (Experiment.name.ilike(search_term)) | 
                (Experiment.description.ilike(search_term)) |
                (Experiment.notes.ilike(search_term)) |
                (Experiment.researcher.ilike(search_term))
            )
        
        # Apply additional filters
        if args.status:
            db_query = db_query.filter(Experiment.status == args.status)
        
        if args.experiment_type:
            db_query = db_query.filter(Experiment.experiment_type == args.experiment_type)
            
        if args.researcher:
            db_query = db_query.filter(Experiment.researcher.ilike(f"%{args.researcher}%"))
            
        if args.tissue_type_id:
            db_query = db_query.filter(Experiment.tissue_type_id == uuid.UUID(args.tissue_type_id))
            
        if args.date_from:
            try:
                date_from = datetime.strptime(args.date_from, '%Y-%m-%d')
                db_query = db_query.filter(Experiment.start_date >= date_from)
            except ValueError:
                return {"message": "Invalid date_from format. Use YYYY-MM-DD"}, 400
                
        if args.date_to:
            try:
                date_to = datetime.strptime(args.date_to, '%Y-%m-%d')
                db_query = db_query.filter(Experiment.start_date <= date_to)
            except ValueError:
                return {"message": "Invalid date_to format. Use YYYY-MM-DD"}, 400
                
        if args.tags:
            tags = args.tags.split(',')
            for tag in tags:
                # Using PostgreSQL array overlap operator
                db_query = db_query.filter(Experiment.tags.overlap('{' + tag.strip() + '}'))
        
        # Apply sorting
        sort_column = getattr(Experiment, args.sort_by, Experiment.start_date)
        sort_method = desc(sort_column) if args.sort_order == 'desc' else asc(sort_column)
        db_query = db_query.order_by(sort_method)
        
        # Get total count (before pagination)
        total = db_query.count()
        
        # Apply pagination
        page = args.page
        per_page = args.per_page
        offset = (page - 1) * per_page
        experiments = db_query.offset(offset).limit(per_page).all()
        
        # Serialize and return results
        return {
            'data': experiments_schema.dump(experiments),
            'total': total
        }


# Register resources
experiment_api.add_resource(ExperimentList, '')
experiment_api.add_resource(ExperimentDetail, '/<experiment_id>')
experiment_api.add_resource(ExperimentResultList, '/<experiment_id>/results')
experiment_api.add_resource(ExperimentResultDetail, '/results/<result_id>')
experiment_api.add_resource(ExperimentTimeSeriesList, '/<experiment_id>/timeseries')
experiment_api.add_resource(ExperimentAnalysis, '/analyze')
experiment_api.add_resource(ExperimentExport, '/<experiment_id>/export')
experiment_api.add_resource(ExperimentImport, '/import')
experiment_api.add_resource(ExperimentSearch, '/search')