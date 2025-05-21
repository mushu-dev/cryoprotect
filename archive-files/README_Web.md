# CryoProtect Analyzer Web Interface

This document provides information about the CryoProtect Analyzer Web Interface, which provides a user-friendly way to interact with the CryoProtect Analyzer API.

## Overview

The CryoProtect Analyzer Web Interface is a responsive web application built with modern web technologies that allows users to:

- View molecules and their properties
- View mixtures and their components
- Create new mixtures
- Add predictions for mixtures
- Record experimental results
- Compare predictions with experimental results

The web interface integrates with the Flask API and uses Supabase for authentication and data storage.

## Technologies Used

- **Frontend**:
  - HTML5, CSS3, JavaScript (ES6+)
  - Bootstrap 5 for responsive design
  - Chart.js for data visualization
  - Supabase JavaScript client for authentication

- **Backend**:
  - Flask API (see README_API.md for details)
  - Supabase for authentication and data storage

## Directory Structure

```
CryoProtect Analyzer/
├── static/
│   ├── css/
│   │   └── styles.css         # Custom CSS styles
│   └── js/
│       ├── app.js             # Main JavaScript file
│       ├── api.js             # API interaction module
│       ├── auth.js            # Authentication module
│       └── charts.js          # Data visualization module
├── templates/
│   ├── index.html             # Main page template
│   ├── login.html             # Login page template
│   ├── molecules.html         # Molecules page template
│   ├── mixtures.html          # Mixtures page template
│   ├── predictions.html       # Predictions page template
│   ├── experiments.html       # Experiments page template
│   └── comparisons.html       # Comparisons page template
└── app.py                     # Flask application (includes both API and web interface)
```

## Features

### Authentication

The web interface uses Supabase Authentication to secure access to certain features. Users need to log in to:

- Import molecules from PubChem
- Create new mixtures
- Add predictions
- Record experiments
- Compare results

### Molecules

- View a list of all molecules with their basic properties
- View detailed information about a specific molecule
- Import molecules from PubChem by CID (Compound ID)
- Visualize molecule properties with charts

### Mixtures

- View a list of all mixtures
- View detailed information about a specific mixture, including its components
- Create new mixtures by combining multiple molecules
- Visualize mixture composition with pie charts

### Predictions

- View predictions for a specific mixture
- Add new predictions for various properties
- Specify confidence levels and calculation methods

### Experiments

- View experimental results for a specific mixture
- Record new experimental results
- Document experimental conditions and dates

### Comparisons

- Compare predictions with experimental results for a specific property
- Visualize the comparison with charts
- Calculate difference and percent error
- Color-coded results based on accuracy

## JavaScript Modules

### app.js

The main JavaScript file that initializes the application and handles page-specific functionality.

Key functions:
- `initializeAuth()`: Initializes authentication
- `initializeCurrentPage()`: Initializes page-specific functionality
- `showToast()`: Displays toast notifications
- Page-specific initialization functions

### api.js

Handles communication with the Flask API.

Key functions:
- `getMolecules()`: Get all molecules
- `getMolecule(moleculeId)`: Get a specific molecule
- `importMolecule(cid)`: Import a molecule from PubChem
- `getMixtures()`: Get all mixtures
- `createMixture(mixtureData)`: Create a new mixture
- `getPredictions(mixtureId)`: Get predictions for a mixture
- `createPrediction(mixtureId, predictionData)`: Add a prediction
- `getExperiments(mixtureId)`: Get experiments for a mixture
- `createExperiment(mixtureId, experimentData)`: Record an experiment
- `getComparison(mixtureId, propertyName)`: Compare prediction with experiment

### auth.js

Handles user authentication with Supabase.

Key functions:
- `signIn(email, password)`: Sign in with email and password
- `signOut()`: Sign out the current user
- `checkSession()`: Check if the user is signed in
- `getToken()`: Get the current authentication token

### charts.js

Handles data visualization using Chart.js.

Key functions:
- `createMoleculePropertiesChart()`: Create a chart for molecule properties
- `createPropertyDistributionChart()`: Create a chart for property distribution
- `createMixtureCompositionChart()`: Create a chart for mixture composition
- `createComparisonChart()`: Create a chart comparing prediction with experiment
- `createMoleculeRadarChart()`: Create a radar chart for molecule properties

## Integration with Flask API

The web interface integrates with the Flask API to fetch and manipulate data. The API endpoints are defined in the `api/__init__.py` file and implemented in the `api/resources.py` file.

See README_API.md for more details about the API.

## Authentication Flow

1. User enters email and password on the login page
2. The auth.js module sends the credentials to Supabase
3. If authentication is successful, Supabase returns a JWT token
4. The token is stored in localStorage
5. The token is included in the Authorization header for API requests
6. The API validates the token and grants access to protected endpoints

## Responsive Design

The web interface is built with Bootstrap 5 and is fully responsive, providing an optimal viewing experience across a wide range of devices:

- Desktop computers
- Laptops
- Tablets
- Mobile phones

## Data Visualization

The web interface uses Chart.js to visualize data in various formats:

- Bar charts for molecule properties
- Pie charts for mixture composition
- Radar charts for molecule property profiles
- Comparison charts for predictions vs. experiments

## Error Handling

The web interface includes comprehensive error handling:

- API request errors are caught and displayed to the user
- Authentication errors are handled gracefully
- Form validation prevents invalid data submission
- Toast notifications provide feedback to the user

## Browser Compatibility

The web interface is compatible with modern browsers:

- Google Chrome (latest)
- Mozilla Firefox (latest)
- Microsoft Edge (latest)
- Safari (latest)

## Future Enhancements

Potential future enhancements for the web interface:

1. User registration and profile management
2. Advanced search and filtering options
3. Batch operations for molecules and mixtures
4. Export data to CSV or Excel
5. Interactive molecular structure visualization
6. Real-time collaboration features
7. Dark mode theme option
8. Mobile app version using Progressive Web App (PWA) technology

## Development

### Prerequisites

- Node.js and npm (for frontend development)
- Python 3.6+ (for Flask API)
- Supabase project with the CryoProtect schema applied

### Setup

1. Clone the repository
2. Install the required Python packages:
   ```
   pip install -r requirements.txt
   ```
3. Create a `.env` file with your Supabase credentials
4. Run the Flask application:
   ```
   python app.py
   ```
5. Access the web interface at http://localhost:5000

### Adding New Pages

1. Create a new HTML template in the `templates` directory
2. Add the necessary JavaScript code to `app.js`
3. Update the navigation menu in all templates
4. Add any required API endpoints to the Flask API

## License

This project is licensed under the MIT License - see the LICENSE file for details.