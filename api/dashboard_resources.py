"""
CryoProtect Analyzer - Dashboard Resources

This module provides API endpoints for the dashboard.
"""

from flask import request, jsonify, current_app
from api.utils import (_handle_json_serialization, handle_supabase_error,
                      handle_error, format_response, token_required, get_user_id)
from flask_restful import Resource, fields, marshal_with
from datetime import datetime, timedelta
import random
from marshmallow import ValidationError
from api.models import Molecule, Mixture, Protocol, Experiment, MolecularProperty, ProjectSchema, ProjectExperimentSchema

# Field definitions for marshal_with
dashboard_fields = {
    'keyMetrics': fields.Raw,
    'topMolecules': fields.Raw,
    'topMixtures': fields.Raw,
    'recentProtocols': fields.Raw,
    'propertyDistribution': fields.Raw,
    'recentActivity': fields.Raw,
    'comparisonData': fields.Raw
}

project_dashboard_fields = {
    'projects': fields.Raw,
    'experiments': fields.Raw,
    'recentActivity': fields.Raw,
    'frequentlyUsed': fields.Raw,
    'savedSearches': fields.Raw
}

project_list_fields = fields.List(fields.Raw)
project_detail_fields = fields.Raw
project_experiment_list_fields = fields.List(fields.Raw)
project_experiment_fields = fields.Raw

class DashboardResource(Resource):
    """Resource for dashboard data."""
    
    @token_required
    @marshal_with(dashboard_fields)
    def post(self):
        """
        Get dashboard data based on filters.
        
        Returns:
            dict: Dashboard data
        """
        try:
            # Get filter parameters from request
            data = request.get_json() or {}
            date_range = data.get('dateRange', 'all')
            category = data.get('category', 'all')
            
            # Generate dashboard data
            dashboard_data = self._generate_dashboard_data(date_range, category)
            
            return dashboard_data
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Generating dashboard data",
                log_level='error',
                return_response=True
            )
            return error_response, error_status
    
    def _generate_dashboard_data(self, date_range, category):
        """
        Generate dashboard data based on filters.

        Args:
            date_range (str): Date range filter ('all', 'today', 'week', 'month', 'year')
            category (str): Category filter ('all', 'molecules', 'mixtures', 'protocols', 'experiments')

        Returns:
            dict: Dashboard data
        """
        # Query the database for key metrics and top molecules
        try:
            key_metrics = {
                'molecules': {
                    'value': Molecule.count(),
                    'label': 'Total Molecules'
                },
                'mixtures': {
                    'value': Mixture.count(),
                    'label': 'Total Mixtures'
                },
                'protocols': {
                    'value': Protocol.count(),
                    'label': 'Total Protocols'
                },
                'experiments': {
                    'value': Experiment.count(),
                    'label': 'Total Experiments'
                }
            }
        except Exception as e:
            error_msg, _ = handle_error(
                e,
                context="Fetching key metrics",
                log_level='error',
                return_response=False
            )
            key_metrics = {
                'molecules': {'value': 0, 'label': 'Total Molecules'},
                'mixtures': {'value': 0, 'label': 'Total Mixtures'},
                'protocols': {'value': 0, 'label': 'Total Protocols'},
                'experiments': {'value': 0, 'label': 'Total Experiments'}
            }

        # Get top molecules by "Cryoprotection Score"
        try:
            # Fetch a batch of molecules (limit to 50 for performance)
            molecules = Molecule.get_all(limit=50)
            scored_molecules = []
            for mol in molecules:
                try:
                    score_prop = MolecularProperty.get_property(mol['id'], "Cryoprotection Score")
                    score = score_prop.get("numeric_value") if score_prop else None
                    if score is not None:
                        scored_molecules.append((mol, score))
                except Exception as e:
                    current_app.logger.warning(f"Error fetching score for molecule {mol.get('id')}: {str(e)}")
            # Sort by score descending and take top 5
            top_scored = sorted(scored_molecules, key=lambda x: x[1], reverse=True)[:5]
            top_molecules = []
            for mol, score in top_scored:
                # Fetch additional properties for dashboard
                try:
                    props = {}
                    # Example: fetch some properties, fallback to None if not found
                    mw_prop = MolecularProperty.get_property(mol['id'], "Molecular Weight")
                    logp_prop = MolecularProperty.get_property(mol['id'], "LogP")
                    hbd_prop = MolecularProperty.get_property(mol['id'], "Hydrogen Bond Donors")
                    hba_prop = MolecularProperty.get_property(mol['id'], "Hydrogen Bond Acceptors")
                    rb_prop = MolecularProperty.get_property(mol['id'], "Rotatable Bonds")
                    tpsa_prop = MolecularProperty.get_property(mol['id'], "TPSA")
                    # Add more as needed

                    props['Molecular Weight'] = mw_prop.get("numeric_value") if mw_prop else None
                    props['LogP'] = logp_prop.get("numeric_value") if logp_prop else None
                    props['Hydrogen Bond Donors'] = hbd_prop.get("numeric_value") if hbd_prop else None
                    props['Hydrogen Bond Acceptors'] = hba_prop.get("numeric_value") if hba_prop else None
                    props['Rotatable Bonds'] = rb_prop.get("numeric_value") if rb_prop else None
                    props['TPSA'] = tpsa_prop.get("numeric_value") if tpsa_prop else None
                    props['Total Score'] = score

                    top_molecules.append({
                        'id': mol.get('id'),
                        'name': mol.get('name'),
                        'formula': mol.get('molecular_formula'),
                        'smiles': mol.get('smiles'),
                        'cid': mol.get('cid'),
                        'properties': props
                    })
                except Exception as e:
                    current_app.logger.warning(f"Error fetching properties for molecule {mol.get('id')}: {str(e)}")
        except Exception as e:
            current_app.logger.error(f"Error fetching top molecules: {str(e)}")
            top_molecules = []

        # Top mixtures
        top_mixtures = [
            {
                'id': '1',
                'name': 'Standard Mixture',
                'components': [
                    {'name': 'Glycerol', 'concentration': 40, 'concentration_unit': '%'},
                    {'name': 'DMSO', 'concentration': 30, 'concentration_unit': '%'},
                    {'name': 'Ethylene Glycol', 'concentration': 20, 'concentration_unit': '%'},
                    {'name': 'Trehalose', 'concentration': 10, 'concentration_unit': '%'}
                ],
                'properties': {
                    'Cryoprotection Score': 85,
                    'Toxicity Score': 80,
                    'Viability Score': 75,
                    'Recovery Score': 70,
                    'Osmolarity': 1200,
                    'pH': 7.2
                }
            },
            {
                'id': '2',
                'name': 'Optimized Mixture',
                'components': [
                    {'name': 'Glycerol', 'concentration': 25, 'concentration_unit': '%'},
                    {'name': 'DMSO', 'concentration': 25, 'concentration_unit': '%'},
                    {'name': 'Propylene Glycol', 'concentration': 25, 'concentration_unit': '%'},
                    {'name': 'Trehalose', 'concentration': 25, 'concentration_unit': '%'}
                ],
                'properties': {
                    'Cryoprotection Score': 90,
                    'Toxicity Score': 85,
                    'Viability Score': 85,
                    'Recovery Score': 80,
                    'Osmolarity': 1300,
                    'pH': 7.0
                }
            },
            {
                'id': '3',
                'name': 'Experimental Mixture',
                'components': [
                    {'name': 'Glycerol', 'concentration': 50, 'concentration_unit': '%'},
                    {'name': 'DMSO', 'concentration': 10, 'concentration_unit': '%'},
                    {'name': 'Trehalose', 'concentration': 40, 'concentration_unit': '%'}
                ],
                'properties': {
                    'Cryoprotection Score': 88,
                    'Toxicity Score': 90,
                    'Viability Score': 80,
                    'Recovery Score': 75,
                    'Osmolarity': 1250,
                    'pH': 7.4
                }
            }
        ]
        
        # Recent protocols
        recent_protocols = [
            {
                'id': '1',
                'name': 'Standard Protocol',
                'steps': [
                    {'step_number': 1, 'temperature': 20, 'concentration': 10, 'duration': 10, 'description': 'Initial exposure'},
                    {'step_number': 2, 'temperature': 10, 'concentration': 20, 'duration': 15, 'description': 'Gradual cooling'},
                    {'step_number': 3, 'temperature': 0, 'concentration': 30, 'duration': 20, 'description': 'Pre-freezing'},
                    {'step_number': 4, 'temperature': -10, 'concentration': 40, 'duration': 25, 'description': 'Freezing'},
                    {'step_number': 5, 'temperature': -20, 'concentration': 40, 'duration': 30, 'description': 'Storage'}
                ],
                'properties': {
                    'Total Duration': 100,
                    'Cooling Rate': 0.4,
                    'Success Rate': 75
                }
            },
            {
                'id': '2',
                'name': 'Slow Cooling Protocol',
                'steps': [
                    {'step_number': 1, 'temperature': 20, 'concentration': 10, 'duration': 20, 'description': 'Initial exposure'},
                    {'step_number': 2, 'temperature': 15, 'concentration': 15, 'duration': 20, 'description': 'Slow cooling 1'},
                    {'step_number': 3, 'temperature': 10, 'concentration': 20, 'duration': 20, 'description': 'Slow cooling 2'},
                    {'step_number': 4, 'temperature': 5, 'concentration': 25, 'duration': 20, 'description': 'Slow cooling 3'},
                    {'step_number': 5, 'temperature': 0, 'concentration': 30, 'duration': 20, 'description': 'Pre-freezing'},
                    {'step_number': 6, 'temperature': -5, 'concentration': 35, 'duration': 20, 'description': 'Freezing 1'},
                    {'step_number': 7, 'temperature': -10, 'concentration': 40, 'duration': 20, 'description': 'Freezing 2'},
                    {'step_number': 8, 'temperature': -20, 'concentration': 40, 'duration': 30, 'description': 'Storage'}
                ],
                'properties': {
                    'Total Duration': 170,
                    'Cooling Rate': 0.2,
                    'Success Rate': 85
                }
            },
            {
                'id': '3',
                'name': 'Vitrification Protocol',
                'steps': [
                    {'step_number': 1, 'temperature': 20, 'concentration': 20, 'duration': 5, 'description': 'Initial exposure'},
                    {'step_number': 2, 'temperature': 10, 'concentration': 40, 'duration': 5, 'description': 'Equilibration'},
                    {'step_number': 3, 'temperature': 4, 'concentration': 60, 'duration': 5, 'description': 'Loading'},
                    {'step_number': 4, 'temperature': -196, 'concentration': 60, 'duration': 1, 'description': 'Plunge freezing'},
                    {'step_number': 5, 'temperature': -196, 'concentration': 60, 'duration': 60, 'description': 'Storage'}
                ],
                'properties': {
                    'Total Duration': 76,
                    'Cooling Rate': 200,
                    'Success Rate': 90
                }
            }
        ]
        
        # Property distribution
        property_distribution = {
            'Molecular Weight': {
                'values': [5, 12, 25, 18, 10, 5, 3, 2, 1, 1],
                'labels': ['0-100', '100-200', '200-300', '300-400', '400-500', '500-600', '600-700', '700-800', '800-900', '900+']
            },
            'LogP': {
                'values': [8, 15, 22, 25, 15, 8, 4, 2, 1, 0],
                'labels': ['-5 to -4', '-4 to -3', '-3 to -2', '-2 to -1', '-1 to 0', '0 to 1', '1 to 2', '2 to 3', '3 to 4', '4+']
            },
            'Hydrogen Bond Donors': {
                'values': [10, 25, 30, 20, 10, 5, 0, 0, 0, 0],
                'labels': ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9+']
            },
            'Hydrogen Bond Acceptors': {
                'values': [5, 15, 25, 30, 15, 5, 3, 1, 1, 0],
                'labels': ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9+']
            },
            'Total Score': {
                'values': [2, 5, 10, 15, 25, 20, 15, 5, 2, 1],
                'labels': ['0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '90-100']
            }
        }
        
        # Recent activity
        recent_activity = [
            {
                'title': 'New Molecule Added',
                'description': 'Glycerol phosphate (CID: 752) was added to the database',
                'time': '5 minutes ago'
            },
            {
                'title': 'Mixture Created',
                'description': 'New mixture "Optimized DMSO-Free" was created',
                'time': '1 hour ago'
            },
            {
                'title': 'Protocol Updated',
                'description': 'Protocol "Slow Cooling Protocol" was updated',
                'time': '3 hours ago'
            },
            {
                'title': 'Experiment Completed',
                'description': 'Experiment "Cell Line Viability Test" was completed',
                'time': '1 day ago'
            },
            {
                'title': 'Prediction Generated',
                'description': 'New prediction for "Experimental Mixture" was generated',
                'time': '2 days ago'
            }
        ]
        
        # Comparison data
        comparison_data = {
            'categories': ['Viability', 'Recovery', 'Functionality', 'Morphology', 'Genetic Stability', 'Long-term Storage'],
            'series': [
                {
                    'name': 'Standard Protocol',
                    'values': [75, 70, 65, 80, 90, 85]
                },
                {
                    'name': 'Slow Cooling Protocol',
                    'values': [85, 80, 75, 85, 85, 90]
                },
                {
                    'name': 'Vitrification Protocol',
                    'values': [90, 85, 80, 75, 80, 95]
                }
            ]
        }
        
        # Return dashboard data
        return {
            'keyMetrics': key_metrics,
            'topMolecules': top_molecules,
            'topMixtures': top_mixtures,
            'recentProtocols': recent_protocols,
            'propertyDistribution': property_distribution,
            'recentActivity': recent_activity,
            'comparisonData': comparison_data
        }

def register_resources(api):
    """Register dashboard resources with the API."""
    api.add_resource(DashboardResource, '/api/v1/dashboard')
    api.add_resource(ProjectDashboardResource, '/api/v1/dashboard/projects')
    api.add_resource(ProjectResource, '/api/v1/projects/<string:project_id>')
    api.add_resource(ProjectListResource, '/api/v1/projects')
    api.add_resource(ProjectExperimentResource, '/api/v1/projects/<string:project_id>/experiments')
    api.add_resource(ProjectExperimentDetailResource, '/api/v1/projects/<string:project_id>/experiments/<string:experiment_id>')


class ProjectDashboardResource(Resource):
    """Resource for project dashboard data."""
    
    @token_required
    @marshal_with(project_dashboard_fields)
    def get(self):
        """
        Get project dashboard data for the authenticated user.
        
        Returns:
            dict: Project dashboard data
        """
        try:
            user_id = get_user_id()
            if not user_id:
                error_response, error_status = handle_error(
                    "Authentication required",
                    context="Fetching project dashboard data",
                    log_level='warning',
                    return_response=True,
                    status_code=401
                )
                return error_response, error_status
            
            # Get user's projects
            projects = Project.get_user_projects(user_id)
            
            # Get recent activity across all projects
            recent_activity = []
            for project in projects[:5]:  # Limit to 5 most recent projects
                project_activity = Project.get_recent_activity(project['id'], 3)  # Get 3 most recent activities
                for activity in project_activity:
                    activity['project_name'] = project['name']
                    recent_activity.append(activity)
            
            # Sort by timestamp
            recent_activity.sort(key=lambda x: x['timestamp'], reverse=True)
            recent_activity = recent_activity[:10]  # Limit to 10 most recent activities
            
            # Get experiment statistics
            experiment_stats = {
                'total': 0,
                'by_project': [],
                'by_month': []
            }
            
            for project in projects:
                experiments = ProjectExperiment.get_project_experiments(project['id'])
                experiment_stats['total'] += len(experiments)
                experiment_stats['by_project'].append({
                    'project_name': project['name'],
                    'count': len(experiments)
                })
            
            # Get frequently used tools
            frequently_used = [
                {'name': 'Create Project', 'count': random.randint(5, 20)},
                {'name': 'Add Experiment', 'count': random.randint(10, 30)},
                {'name': 'Compare Results', 'count': random.randint(8, 25)},
                {'name': 'Export Data', 'count': random.randint(3, 15)},
                {'name': 'Share Project', 'count': random.randint(2, 10)}
            ]
            frequently_used.sort(key=lambda x: x['count'], reverse=True)
            
            return {
                'projects': {
                    'total': len(projects),
                    'recent': projects[:5],  # 5 most recent projects
                    'by_status': [
                        {'status': 'Active', 'count': random.randint(5, 15)},
                        {'status': 'Completed', 'count': random.randint(2, 10)},
                        {'status': 'On Hold', 'count': random.randint(1, 5)}
                    ]
                },
                'experiments': experiment_stats,
                'recentActivity': recent_activity,
                'frequentlyUsed': frequently_used,
                'savedSearches': [
                    {'name': 'High Viability Experiments', 'query': 'viability>80'},
                    {'name': 'Recent DMSO Mixtures', 'query': 'component:DMSO date:last30days'},
                    {'name': 'Shared Projects', 'query': 'is:shared'},
                    {'name': 'Completed Experiments', 'query': 'status:completed'}
                ]
            }
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Fetching project dashboard data",
                log_level='error',
                return_response=True
            )
            return error_response, error_status


class ProjectListResource(Resource):
    """Resource for listing and creating projects."""
    
    @token_required
    @marshal_with(project_list_fields)
    def get(self):
        """
        Get a list of projects for the authenticated user.
        
        Returns:
            list: List of projects
        """
        try:
            user_id = get_user_id()
            if not user_id:
                error_response, error_status = handle_error(
                    "Authentication required",
                    context="Fetching projects",
                    log_level='warning',
                    return_response=True,
                    status_code=401
                )
                return error_response, error_status
            
            # Get query parameters
            limit = request.args.get('limit', 100, type=int)
            offset = request.args.get('offset', 0, type=int)
            
            # Get user's projects
            projects = Project.get_user_projects(user_id, limit, offset)
            
            return projects
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Fetching projects",
                log_level='error',
                return_response=True
            )
            return error_response, error_status
    
    @token_required
    @marshal_with(project_detail_fields)
    def post(self):
        """
        Create a new project.
        
        Returns:
            dict: Created project
        """
        try:
            user_id = get_user_id()
            if not user_id:
                error_response, error_status = handle_error(
                    "Authentication required",
                    context="Creating project",
                    log_level='warning',
                    return_response=True,
                    status_code=401
                )
                return error_response, error_status
            
            # Validate request data
            schema = ProjectSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating project data",
                    log_level='warning',
                    return_response=True,
                    status_code=400
                )
                return error_response, error_status
            
            # Create project
            project_data = {
                'name': data['name'],
                'description': data.get('description', ''),
                'is_public': data.get('is_public', False),
                'created_by': user_id
            }
            
            project = Project.create(project_data)
            
            return project, 201
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Creating project",
                log_level='error',
                return_response=True
            )
            return error_response, error_status


class ProjectResource(Resource):
    """Resource for retrieving, updating, and deleting a project."""
    
    @token_required
    @marshal_with(project_detail_fields)
    def get(self, project_id):
        """
        Get a project by ID.
        
        Args:
            project_id: ID of the project
            
        Returns:
            dict: Project data
        """
        try:
            user_id = get_user_id()
            if not user_id:
                error_response, error_status = handle_error(
                    "Authentication required",
                    context="Fetching project",
                    log_level='warning',
                    return_response=True,
                    status_code=401
                )
                return error_response, error_status
            
            # Get project with experiment count
            project = Project.get_with_experiment_count(project_id)
            
            if not project:
                error_response, error_status = handle_error(
                    f"Project with ID {project_id} not found",
                    context="Fetching project",
                    log_level='warning',
                    return_response=True,
                    status_code=404
                )
                return error_response, error_status
            
            # Check if user has access to the project
            if project['created_by'] != user_id and not project['is_public']:
                error_response, error_status = handle_error(
                    "You do not have access to this project",
                    context="Fetching project",
                    log_level='warning',
                    return_response=True,
                    status_code=403
                )
                return error_response, error_status
            
            return project
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Fetching project",
                log_level='error',
                return_response=True
            )
            return error_response, error_status
    
    @token_required
    @marshal_with(project_detail_fields)
    def put(self, project_id):
        """
        Update a project.
        
        Args:
            project_id: ID of the project
            
        Returns:
            dict: Updated project
        """
        try:
            user_id = get_user_id()
            if not user_id:
                error_response, error_status = handle_error(
                    "Authentication required",
                    context="Updating project",
                    log_level='warning',
                    return_response=True,
                    status_code=401
                )
                return error_response, error_status
            
            # Check if project exists and user has access
            project = Project.get(project_id)
            
            if not project:
                error_response, error_status = handle_error(
                    f"Project with ID {project_id} not found",
                    context="Updating project",
                    log_level='warning',
                    return_response=True,
                    status_code=404
                )
                return error_response, error_status
            
            if project['created_by'] != user_id:
                error_response, error_status = handle_error(
                    "You do not have permission to update this project",
                    context="Updating project",
                    log_level='warning',
                    return_response=True,
                    status_code=403
                )
                return error_response, error_status
            
            # Validate request data
            schema = ProjectSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating project data",
                    log_level='warning',
                    return_response=True,
                    status_code=400
                )
                return error_response, error_status
            
            # Update project
            project_data = {
                'name': data['name'],
                'description': data.get('description', project['description']),
                'is_public': data.get('is_public', project['is_public'])
            }
            
            updated_project = Project.update(project_id, project_data)
            
            return updated_project
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Updating project",
                log_level='error',
                return_response=True
            )
            return error_response, error_status
    
    @token_required
    @marshal_with(fields.Raw)
    def delete(self, project_id):
        """
        Delete a project.
        
        Args:
            project_id: ID of the project
            
        Returns:
            dict: Success message
        """
        try:
            user_id = get_user_id()
            if not user_id:
                error_response, error_status = handle_error(
                    "Authentication required",
                    context="Deleting project",
                    log_level='warning',
                    return_response=True,
                    status_code=401
                )
                return error_response, error_status
            
            # Check if project exists and user has access
            project = Project.get(project_id)
            
            if not project:
                error_response, error_status = handle_error(
                    f"Project with ID {project_id} not found",
                    context="Deleting project",
                    log_level='warning',
                    return_response=True,
                    status_code=404
                )
                return error_response, error_status
            
            if project['created_by'] != user_id:
                error_response, error_status = handle_error(
                    "You do not have permission to delete this project",
                    context="Deleting project",
                    log_level='warning',
                    return_response=True,
                    status_code=403
                )
                return error_response, error_status
            
            # Delete project
            Project.delete(project_id)
            
            return {'message': 'Project deleted successfully'}
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Deleting project",
                log_level='error',
                return_response=True
            )
            return error_response, error_status


class ProjectExperimentResource(Resource):
    """Resource for managing experiments in a project."""
    
    @token_required
    @marshal_with(project_experiment_list_fields)
    def get(self, project_id):
        """
        Get all experiments for a project.
        
        Args:
            project_id: ID of the project
            
        Returns:
            list: List of experiments
        """
        try:
            user_id = get_user_id()
            if not user_id:
                error_response, error_status = handle_error(
                    "Authentication required",
                    context="Fetching project experiments",
                    log_level='warning',
                    return_response=True,
                    status_code=401
                )
                return error_response, error_status
            
            # Check if project exists and user has access
            project = Project.get(project_id)
            
            if not project:
                error_response, error_status = handle_error(
                    f"Project with ID {project_id} not found",
                    context="Fetching project experiments",
                    log_level='warning',
                    return_response=True,
                    status_code=404
                )
                return error_response, error_status
            
            if project['created_by'] != user_id and not project['is_public']:
                error_response, error_status = handle_error(
                    "You do not have access to this project",
                    context="Fetching project experiments",
                    log_level='warning',
                    return_response=True,
                    status_code=403
                )
                return error_response, error_status
            
            # Get experiments for the project
            experiments = ProjectExperiment.get_project_experiments(project_id)
            
            return experiments
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Fetching project experiments",
                log_level='error',
                return_response=True
            )
            return error_response, error_status
    
    @token_required
    @marshal_with(project_experiment_fields)
    def post(self, project_id):
        """
        Add an experiment to a project.
        
        Args:
            project_id: ID of the project
            
        Returns:
            dict: Added experiment
        """
        try:
            user_id = get_user_id()
            if not user_id:
                error_response, error_status = handle_error(
                    "Authentication required",
                    context="Adding experiment to project",
                    log_level='warning',
                    return_response=True,
                    status_code=401
                )
                return error_response, error_status
            
            # Check if project exists and user has access
            project = Project.get(project_id)
            
            if not project:
                error_response, error_status = handle_error(
                    f"Project with ID {project_id} not found",
                    context="Adding experiment to project",
                    log_level='warning',
                    return_response=True,
                    status_code=404
                )
                return error_response, error_status
            
            if project['created_by'] != user_id:
                error_response, error_status = handle_error(
                    "You do not have permission to add experiments to this project",
                    context="Adding experiment to project",
                    log_level='warning',
                    return_response=True,
                    status_code=403
                )
                return error_response, error_status
            
            # Validate request data
            schema = ProjectExperimentSchema()
            try:
                data = schema.load(request.json)
            except ValidationError as err:
                error_response, error_status = handle_error(
                    err,
                    context="Validating experiment data",
                    log_level='warning',
                    return_response=True,
                    status_code=400
                )
                return error_response, error_status
            
            # Check if experiment exists
            experiment = Experiment.get(data['experiment_id'])
            if not experiment:
                error_response, error_status = handle_error(
                    f"Experiment with ID {data['experiment_id']} not found",
                    context="Adding experiment to project",
                    log_level='warning',
                    return_response=True,
                    status_code=404
                )
                return error_response, error_status
            
            # Add experiment to project
            project_experiment = ProjectExperiment.add_experiment_to_project(
                project_id,
                data['experiment_id'],
                data.get('notes')
            )
            
            return project_experiment, 201
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Adding experiment to project",
                log_level='error',
                return_response=True
            )
            return error_response, error_status


class ProjectExperimentDetailResource(Resource):
    """Resource for managing a specific experiment in a project."""
    
    @token_required
    @marshal_with(fields.Raw)
    def delete(self, project_id, experiment_id):
        """
        Remove an experiment from a project.
        
        Args:
            project_id: ID of the project
            experiment_id: ID of the experiment
            
        Returns:
            dict: Success message
        """
        try:
            user_id = get_user_id()
            if not user_id:
                error_response, error_status = handle_error(
                    "Authentication required",
                    context="Removing experiment from project",
                    log_level='warning',
                    return_response=True,
                    status_code=401
                )
                return error_response, error_status
            
            # Check if project exists and user has access
            project = Project.get(project_id)
            
            if not project:
                error_response, error_status = handle_error(
                    f"Project with ID {project_id} not found",
                    context="Removing experiment from project",
                    log_level='warning',
                    return_response=True,
                    status_code=404
                )
                return error_response, error_status
            
            if project['created_by'] != user_id:
                error_response, error_status = handle_error(
                    "You do not have permission to remove experiments from this project",
                    context="Removing experiment from project",
                    log_level='warning',
                    return_response=True,
                    status_code=403
                )
                return error_response, error_status
            
            # Remove experiment from project
            ProjectExperiment.remove_experiment_from_project(project_id, experiment_id)
            
            return {'message': 'Experiment removed from project successfully'}
        except Exception as e:
            error_response, error_status = handle_error(
                e,
                context="Removing experiment from project",
                log_level='error',
                return_response=True
            )
            return error_response, error_status