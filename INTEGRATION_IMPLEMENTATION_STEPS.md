# CryoProtect Integration Implementation Steps

Based on validation testing, here are the specific steps needed to complete the backend-frontend integration:

## 1. Implement API Endpoints on Heroku

The main API on Heroku needs the following endpoints:

### Create a simple `app.py` for testing:

```python
from flask import Flask, jsonify
from flask_cors import CORS

app = Flask(__name__)

# Configure CORS to allow requests from all origins
CORS(app, resources={
    r"/*": {
        "origins": "*",
        "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"],
        "allow_headers": ["Content-Type", "Authorization", "X-API-Key"]
    }
})

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        "status": "ok",
        "timestamp": "now",
        "service": "cryoprotect-api"
    })

@app.route('/test-cors', methods=['GET', 'OPTIONS'])
def test_cors():
    """Endpoint for testing CORS configuration"""
    return jsonify({
        "success": True,
        "message": "CORS is configured properly",
        "service": "cryoprotect-api"
    })

@app.route('/api/v1/health/connectivity', methods=['GET'])
def connectivity():
    """API connectivity test endpoint"""
    return jsonify({
        "status": "ok",
        "connectivity": "established",
        "service": "cryoprotect-api"
    })

if __name__ == '__main__':
    app.run(debug=True)
```

### Create a `requirements.txt` file:

```
flask==2.0.1
flask-cors==3.0.10
gunicorn==20.1.0
```

### Create a `Procfile`:

```
web: gunicorn app:app
```

### Deploy to Heroku:

```bash
# Create a new directory for deployment
mkdir -p heroku-deploy
cd heroku-deploy

# Create app.py, requirements.txt, and Procfile

# Initialize git repository
git init
git add .
git commit -m "Add basic API with CORS support"

# Set up Heroku remote
heroku git:remote -a cryoprotect

# Deploy to Heroku
git push heroku main:main -f

# Set CORS environment variables
heroku config:set FRONTEND_URL="https://cryoprotect.netlify.app" -a cryoprotect
heroku config:set ALLOWED_ORIGINS="https://cryoprotect.netlify.app,https://rdkit.cryoprotect.app,https://dynamic-mink-63.convex.cloud" -a cryoprotect

# Verify deployment
curl https://cryoprotect.herokuapp.com/health
```

## 2. Deploy Netlify Changes

Deploy the updated netlify.toml to enable proper redirects:

```bash
cd frontend
netlify deploy --prod
```

After deployment, verify the redirects are working:

```bash
curl https://cryoprotect.netlify.app/api/health
```

## 3. Set Up RDKit Service (when fly.io CLI is available)

Follow these steps to set up the RDKit service:

```bash
# Install fly.io CLI
curl -L https://fly.io/install.sh | sh

# Add fly.io to your PATH (may vary based on installation)
export PATH=$HOME/.fly/bin:$PATH

# Deploy RDKit service
./deploy-rdkit-service.sh
```

## 4. Run Full Validation

After implementing all components, run the full validation test:

```bash
./test-connection.sh
```

## Troubleshooting

If specific tests are still failing:

1. **API connectivity issues**:
   - Check Heroku logs: `heroku logs --tail -a cryoprotect`
   - Verify the app is running: `heroku ps -a cryoprotect`

2. **Netlify redirect issues**:
   - Check netlify.toml formatting
   - Verify netlify.toml was deployed correctly
   - Try clearing Netlify cache

3. **CORS issues**:
   - Verify CORS headers using browser developer tools
   - Check for proper CORS configuration in Flask app
   - Test with curl including Origin header

Refer to the full `INTEGRATION_TROUBLESHOOTING_GUIDE.md` for more detailed troubleshooting steps.