#!/usr/bin/env python3
"""
Enhance OpenAPI YAML file for the CryoProtect API.

This script enhances an existing OpenAPI YAML file with additional endpoints
and schemas based on the API resources. It doesn't require the full application
to run, making it more reliable in different environments.

Usage:
    python enhance_openapi.py [input_path] [output_path]

Arguments:
    input_path: Path to the input YAML file (default: docs/api/openapi.yaml)
    output_path: Path to save the enhanced YAML file (default: same as input_path)
"""

import os
import sys
import logging
import yaml
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def load_yaml(file_path):
    """Load YAML file."""
    try:
        with open(file_path, 'r') as f:
            return yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Error loading YAML file: {str(e)}")
        return None

def save_yaml(data, file_path):
    """Save YAML file."""
    try:
        with open(file_path, 'w') as f:
            yaml.dump(data, f, sort_keys=False, default_flow_style=False)
        return True
    except Exception as e:
        logger.error(f"Error saving YAML file: {str(e)}")
        return False

def enhance_openapi_spec(spec):
    """Enhance OpenAPI specification with additional endpoints and schemas."""
    # Ensure all required components are present
    if 'components' not in spec:
        spec['components'] = {}
    if 'schemas' not in spec['components']:
        spec['components']['schemas'] = {}
    if 'paths' not in spec:
        spec['paths'] = {}
    
    # Add or update mixture analysis endpoints
    add_mixture_analysis_endpoints(spec)
    
    # Add or update predictive models endpoints
    add_predictive_models_endpoints(spec)
    
    # Add or update scoring endpoints
    add_scoring_endpoints(spec)
    
    # Add or update schemas
    add_schemas(spec)
    
    return spec

def add_mixture_analysis_endpoints(spec):
    """Add mixture analysis endpoints to the specification."""
    paths = spec['paths']
    
    # Mixture Properties endpoint
    if '/api/v1/mixture-analysis/properties/{mixture_id}' not in paths:
        paths['/api/v1/mixture-analysis/properties/{mixture_id}'] = {
            'get': {
                'summary': 'Get predicted properties for a mixture',
                'description': 'Calculates various physicochemical properties of a mixture based on its components and their concentrations.',
                'operationId': 'getMixtureProperties',
                'tags': ['Mixture Analysis'],
                'parameters': [
                    {
                        'name': 'mixture_id',
                        'in': 'path',
                        'required': True,
                        'schema': {'type': 'string'},
                        'description': 'ID of the mixture to analyze'
                    }
                ],
                'security': [{'bearerAuth': []}],
                'responses': {
                    '200': {
                        'description': 'Predicted mixture properties',
                        'content': {
                            'application/json': {
                                'schema': {
                                    'type': 'object',
                                    'properties': {
                                        'mixture_id': {'type': 'string'},
                                        'mixture_name': {'type': 'string'},
                                        'properties': {'type': 'object', 'additionalProperties': True},
                                        'raw_properties': {'type': 'object', 'additionalProperties': True}
                                    },
                                    'required': ['mixture_id', 'properties']
                                }
                            }
                        }
                    },
                    '400': {'$ref': '#/components/responses/BadRequest'},
                    '401': {'$ref': '#/components/responses/Unauthorized'},
                    '404': {'$ref': '#/components/responses/NotFound'},
                    '500': {'$ref': '#/components/responses/ServerError'}
                }
            }
        }
    
    # Mixture Compatibility endpoint
    if '/api/v1/mixture-analysis/compatibility/{mixture_id}' not in paths:
        paths['/api/v1/mixture-analysis/compatibility/{mixture_id}'] = {
            'get': {
                'summary': 'Analyze compatibility between mixture components',
                'description': 'Evaluates chemical and physical compatibility between all components in the mixture.',
                'operationId': 'getMixtureCompatibility',
                'tags': ['Mixture Analysis'],
                'parameters': [
                    {
                        'name': 'mixture_id',
                        'in': 'path',
                        'required': True,
                        'schema': {'type': 'string'},
                        'description': 'ID of the mixture to analyze'
                    }
                ],
                'security': [{'bearerAuth': []}],
                'responses': {
                    '200': {
                        'description': 'Compatibility analysis results',
                        'content': {
                            'application/json': {
                                'schema': {
                                    'type': 'object',
                                    'properties': {
                                        'mixture_id': {'type': 'string'},
                                        'mixture_name': {'type': 'string'},
                                        'compatibility': {'type': 'object', 'additionalProperties': True}
                                    },
                                    'required': ['mixture_id', 'compatibility']
                                }
                            }
                        }
                    },
                    '400': {'$ref': '#/components/responses/BadRequest'},
                    '401': {'$ref': '#/components/responses/Unauthorized'},
                    '404': {'$ref': '#/components/responses/NotFound'},
                    '500': {'$ref': '#/components/responses/ServerError'}
                }
            }
        }
    
    # Mixture Synergy endpoint
    if '/api/v1/mixture-analysis/synergy/{mixture_id}' not in paths:
        paths['/api/v1/mixture-analysis/synergy/{mixture_id}'] = {
            'get': {
                'summary': 'Analyze synergistic or antagonistic effects in a mixture',
                'description': 'Evaluates how components interact to produce effects that are greater than or less than the sum of their individual effects.',
                'operationId': 'getMixtureSynergy',
                'tags': ['Mixture Analysis'],
                'parameters': [
                    {
                        'name': 'mixture_id',
                        'in': 'path',
                        'required': True,
                        'schema': {'type': 'string'},
                        'description': 'ID of the mixture to analyze'
                    }
                ],
                'security': [{'bearerAuth': []}],
                'responses': {
                    '200': {
                        'description': 'Synergy analysis results',
                        'content': {
                            'application/json': {
                                'schema': {
                                    'type': 'object',
                                    'properties': {
                                        'mixture_id': {'type': 'string'},
                                        'mixture_name': {'type': 'string'},
                                        'synergy': {'type': 'object', 'additionalProperties': True}
                                    },
                                    'required': ['mixture_id', 'synergy']
                                }
                            }
                        }
                    },
                    '400': {'$ref': '#/components/responses/BadRequest'},
                    '401': {'$ref': '#/components/responses/Unauthorized'},
                    '404': {'$ref': '#/components/responses/NotFound'},
                    '500': {'$ref': '#/components/responses/ServerError'}
                }
            }
        }
    
    # Mixture Optimization endpoint
    if '/api/v1/mixture-analysis/optimize/{mixture_id}' not in paths:
        paths['/api/v1/mixture-analysis/optimize/{mixture_id}'] = {
            'post': {
                'summary': 'Optimize the composition of a mixture',
                'description': 'Adjusts component concentrations to achieve a target property value or maximize overall effectiveness.',
                'operationId': 'optimizeMixture',
                'tags': ['Mixture Analysis'],
                'parameters': [
                    {
                        'name': 'mixture_id',
                        'in': 'path',
                        'required': True,
                        'schema': {'type': 'string'},
                        'description': 'ID of the mixture to optimize'
                    }
                ],
                'requestBody': {
                    'required': True,
                    'content': {
                        'application/json': {
                            'schema': {
                                'type': 'object',
                                'properties': {
                                    'target_property': {
                                        'type': 'string',
                                        'description': 'Property to optimize for',
                                        'default': 'Cryoprotection Score'
                                    },
                                    'target_value': {
                                        'type': 'number',
                                        'description': 'Target value for the property'
                                    },
                                    'constraints': {
                                        'type': 'object',
                                        'description': 'Constraints for optimization',
                                        'additionalProperties': True
                                    }
                                }
                            }
                        }
                    }
                },
                'security': [{'bearerAuth': []}],
                'responses': {
                    '200': {
                        'description': 'Optimized mixture composition',
                        'content': {
                            'application/json': {
                                'schema': {
                                    'type': 'object',
                                    'properties': {
                                        'mixture_id': {'type': 'string'},
                                        'mixture_name': {'type': 'string'},
                                        'original_composition': {'type': 'object', 'additionalProperties': True},
                                        'optimized_composition': {'type': 'object', 'additionalProperties': True},
                                        'target_property': {'type': 'string'},
                                        'target_value': {'type': 'number'},
                                        'achieved_value': {'type': 'number'},
                                        'improvement': {'type': 'number'}
                                    },
                                    'required': ['mixture_id', 'optimized_composition']
                                }
                            }
                        }
                    },
                    '400': {'$ref': '#/components/responses/BadRequest'},
                    '401': {'$ref': '#/components/responses/Unauthorized'},
                    '404': {'$ref': '#/components/responses/NotFound'},
                    '500': {'$ref': '#/components/responses/ServerError'}
                }
            }
        }

def add_predictive_models_endpoints(spec):
    """Add predictive models endpoints to the specification."""
    paths = spec['paths']
    
    # List available models endpoint
    if '/api/v1/predictive-models' not in paths:
        paths['/api/v1/predictive-models'] = {
            'get': {
                'summary': 'Get a list of available predictive models',
                'description': 'Retrieves information about all available predictive models.',
                'operationId': 'getAvailableModels',
                'tags': ['Predictive Models'],
                'responses': {
                    '200': {
                        'description': 'List of available models',
                        'content': {
                            'application/json': {
                                'schema': {
                                    'type': 'array',
                                    'items': {
                                        'type': 'object',
                                        'properties': {
                                            'property_name': {'type': 'string'},
                                            'algorithm': {'type': 'string'},
                                            'algorithm_name': {'type': 'string'},
                                            'hyperparameters': {'type': 'object', 'additionalProperties': True},
                                            'trained_date': {'type': 'string', 'format': 'date-time'},
                                            'metrics': {'type': 'object', 'additionalProperties': True},
                                            'feature_importance': {'type': 'object', 'additionalProperties': True}
                                        }
                                    }
                                }
                            }
                        }
                    },
                    '500': {'$ref': '#/components/responses/ServerError'}
                }
            },
            'post': {
                'summary': 'Train a new predictive model',
                'description': 'Trains a machine learning model to predict the specified property.',
                'operationId': 'trainModel',
                'tags': ['Predictive Models'],
                'requestBody': {
                    'required': True,
                    'content': {
                        'application/json': {
                            'schema': {
                                'type': 'object',
                                'properties': {
                                    'property_name': {
                                        'type': 'string',
                                        'description': 'Property to predict'
                                    },
                                    'algorithm': {
                                        'type': 'string',
                                        'description': 'Machine learning algorithm to use',
                                        'default': 'random_forest',
                                        'enum': ['random_forest', 'gradient_boosting', 'svm', 'neural_network']
                                    },
                                    'hyperparameters': {
                                        'type': 'object',
                                        'description': 'Algorithm-specific hyperparameters',
                                        'additionalProperties': True
                                    }
                                },
                                'required': ['property_name']
                            }
                        }
                    }
                },
                'security': [{'bearerAuth': []}],
                'responses': {
                    '201': {
                        'description': 'Trained model information',
                        'content': {
                            'application/json': {
                                'schema': {
                                    'type': 'object',
                                    'properties': {
                                        'property_name': {'type': 'string'},
                                        'algorithm': {'type': 'string'},
                                        'algorithm_name': {'type': 'string'},
                                        'hyperparameters': {'type': 'object', 'additionalProperties': True},
                                        'trained_date': {'type': 'string', 'format': 'date-time'},
                                        'metrics': {'type': 'object', 'additionalProperties': True},
                                        'feature_importance': {'type': 'object', 'additionalProperties': True}
                                    }
                                }
                            }
                        }
                    },
                    '400': {'$ref': '#/components/responses/BadRequest'},
                    '401': {'$ref': '#/components/responses/Unauthorized'},
                    '500': {'$ref': '#/components/responses/ServerError'}
                }
            }
        }
    
    # Mixture prediction endpoint
    if '/api/v1/predictive-models/predict/{mixture_id}' not in paths:
        paths['/api/v1/predictive-models/predict/{mixture_id}'] = {
            'get': {
                'summary': 'Get a prediction for a mixture',
                'description': 'Predicts the cryoprotection effectiveness of a mixture using the specified machine learning algorithm.',
                'operationId': 'predictMixture',
                'tags': ['Predictive Models'],
                'parameters': [
                    {
                        'name': 'mixture_id',
                        'in': 'path',
                        'required': True,
                        'schema': {'type': 'string'},
                        'description': 'ID of the mixture to predict for'
                    },
                    {
                        'name': 'algorithm',
                        'in': 'query',
                        'schema': {
                            'type': 'string',
                            'enum': ['random_forest', 'gradient_boosting', 'svm', 'neural_network'],
                            'default': 'random_forest'
                        },
                        'description': 'Algorithm to use for prediction'
                    }
                ],
                'responses': {
                    '200': {
                        'description': 'Prediction results',
                        'content': {
                            'application/json': {
                                'schema': {
                                    'type': 'object',
                                    'properties': {
                                        'property_name': {'type': 'string'},
                                        'prediction': {'type': 'number'},
                                        'confidence': {'type': 'number'},
                                        'confidence_interval': {
                                            'type': 'object',
                                            'properties': {
                                                'lower': {'type': 'number'},
                                                'upper': {'type': 'number'}
                                            }
                                        },
                                        'algorithm': {'type': 'string'},
                                        'algorithm_name': {'type': 'string'},
                                        'feature_importance': {'type': 'object', 'additionalProperties': True}
                                    }
                                }
                            }
                        }
                    },
                    '400': {'$ref': '#/components/responses/BadRequest'},
                    '500': {'$ref': '#/components/responses/ServerError'}
                }
            }
        }

def add_scoring_endpoints(spec):
    """Add scoring endpoints to the specification."""
    paths = spec['paths']
    
    # Molecule score endpoint
    if '/api/v1/scoring/molecule' not in paths:
        paths['/api/v1/scoring/molecule'] = {
            'post': {
                'summary': 'Calculate cryoprotection score for a molecule',
                'description': 'Evaluates a molecule\'s potential effectiveness as a cryoprotectant based on its physicochemical properties and structural features.',
                'operationId': 'scoreMolecule',
                'tags': ['Scoring'],
                'requestBody': {
                    'required': True,
                    'content': {
                        'application/json': {
                            'schema': {
                                'type': 'object',
                                'properties': {
                                    'molecule_data': {
                                        'type': 'string',
                                        'description': 'Molecule data in the specified format'
                                    },
                                    'input_format': {
                                        'type': 'string',
                                        'description': 'Format of the input data',
                                        'default': 'smiles',
                                        'enum': ['smiles', 'inchi', 'mol']
                                    },
                                    'store_result': {
                                        'type': 'boolean',
                                        'description': 'Whether to store the result in the database',
                                        'default': False
                                    }
                                },
                                'required': ['molecule_data']
                            }
                        }
                    }
                },
                'responses': {
                    '200': {
                        'description': 'Cryoprotection score results',
                        'content': {
                            'application/json': {
                                'schema': {
                                    'type': 'object',
                                    'properties': {
                                        'overall_score': {'type': 'integer'},
                                        'component_scores': {
                                            'type': 'object',
                                            'properties': {
                                                'hydrogen_bonding': {'type': 'integer'},
                                                'logp': {'type': 'integer'},
                                                'molecular_size': {'type': 'integer'},
                                                'tpsa': {'type': 'integer'},
                                                'functional_groups': {'type': 'integer'},
                                                'permeability': {'type': 'integer'}
                                            }
                                        },
                                        'properties': {'type': 'object', 'additionalProperties': True}
                                    },
                                    'required': ['overall_score', 'component_scores']
                                }
                            }
                        }
                    },
                    '400': {'$ref': '#/components/responses/BadRequest'},
                    '500': {'$ref': '#/components/responses/ServerError'}
                }
            }
        }
    
    # Mixture score endpoint
    if '/api/v1/scoring/mixture/{mixture_id}' not in paths:
        paths['/api/v1/scoring/mixture/{mixture_id}'] = {
            'post': {
                'summary': 'Calculate cryoprotection score for a mixture',
                'description': 'Evaluates a mixture\'s potential effectiveness as a cryoprotectant based on its components and their properties.',
                'operationId': 'scoreMixture',
                'tags': ['Scoring'],
                'parameters': [
                    {
                        'name': 'mixture_id',
                        'in': 'path',
                        'required': True,
                        'schema': {'type': 'string'},
                        'description': 'ID of the mixture to score'
                    }
                ],
                'requestBody': {
                    'content': {
                        'application/json': {
                            'schema': {
                                'type': 'object',
                                'properties': {
                                    'store_result': {
                                        'type': 'boolean',
                                        'description': 'Whether to store the result in the database',
                                        'default': True
                                    }
                                }
                            }
                        }
                    }
                },
                'security': [{'bearerAuth': []}],
                'responses': {
                    '200': {
                        'description': 'Mixture cryoprotection score results',
                        'content': {
                            'application/json': {
                                'schema': {
                                    'type': 'object',
                                    'properties': {
                                        'mixture_id': {'type': 'string'},
                                        'name': {'type': 'string'},
                                        'overall_score': {'type': 'integer'},
                                        'component_scores': {
                                            'type': 'array',
                                            'items': {
                                                'type': 'object',
                                                'properties': {
                                                    'molecule_id': {'type': 'string'},
                                                    'name': {'type': 'string'},
                                                    'concentration': {'type': 'number'},
                                                    'concentration_unit': {'type': 'string'},
                                                    'score': {'type': 'integer'}
                                                }
                                            }
                                        }
                                    },
                                    'required': ['mixture_id', 'overall_score', 'component_scores']
                                }
                            }
                        }
                    },
                    '400': {'$ref': '#/components/responses/BadRequest'},
                    '401': {'$ref': '#/components/responses/Unauthorized'},
                    '404': {'$ref': '#/components/responses/NotFound'},
                    '500': {'$ref': '#/components/responses/ServerError'}
                }
            }
        }

def add_schemas(spec):
    """Add or update schemas in the specification."""
    schemas = spec['components']['schemas']
    
    # Add common response schemas if not present
    if 'components' not in spec:
        spec['components'] = {}
    
    if 'responses' not in spec['components']:
        spec['components']['responses'] = {
            'BadRequest': {
                'description': 'Invalid request',
                'content': {
                    'application/json': {
                        'schema': {'$ref': '#/components/schemas/Error'}
                    }
                }
            },
            'Unauthorized': {
                'description': 'Authentication required',
                'content': {
                    'application/json': {
                        'schema': {'$ref': '#/components/schemas/Error'}
                    }
                }
            },
            'Forbidden': {
                'description': 'Permission denied',
                'content': {
                    'application/json': {
                        'schema': {'$ref': '#/components/schemas/Error'}
                    }
                }
            },
            'NotFound': {
                'description': 'Resource not found',
                'content': {
                    'application/json': {
                        'schema': {'$ref': '#/components/schemas/Error'}
                    }
                }
            },
            'ServerError': {
                'description': 'Server error',
                'content': {
                    'application/json': {
                        'schema': {'$ref': '#/components/schemas/Error'}
                    }
                }
            }
        }

def main():
    """Enhance the OpenAPI YAML file."""
    # Get the input and output paths from command line arguments or use defaults
    input_path = sys.argv[1] if len(sys.argv) > 1 else 'docs/api/openapi.yaml'
    output_path = sys.argv[2] if len(sys.argv) > 2 else input_path
    
    logger.info(f"Loading OpenAPI YAML file from {input_path}...")
    spec = load_yaml(input_path)
    
    if not spec:
        logger.error("Failed to load OpenAPI YAML file")
        return 1
    
    logger.info("Enhancing OpenAPI specification...")
    enhanced_spec = enhance_openapi_spec(spec)
    
    logger.info(f"Saving enhanced OpenAPI YAML file to {output_path}...")
    success = save_yaml(enhanced_spec, output_path)
    
    if success:
        logger.info(f"OpenAPI YAML file enhanced successfully at {output_path}")
        
        # Get the file size
        file_size = Path(output_path).stat().st_size
        logger.info(f"File size: {file_size / 1024:.2f} KB")
        
        return 0
    else:
        logger.error("Failed to save enhanced OpenAPI YAML file")
        return 1

if __name__ == '__main__':
    sys.exit(main())