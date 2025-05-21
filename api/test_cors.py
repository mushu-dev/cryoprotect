"""
Simple Flask app to test CORS configuration
"""

from flask import Flask, jsonify
from cors_config import configure_cors

app = Flask(__name__)
configure_cors(app)

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        "status": "ok",
        "timestamp": "now",
        "database": "connected"
    })

@app.route('/test-cors', methods=['GET', 'OPTIONS'])
def test_cors():
    """Endpoint for testing CORS configuration"""
    return jsonify({
        "success": True,
        "message": "CORS is configured properly",
        "service": "main-api"
    })

@app.route('/api/v1/health/connectivity', methods=['GET'])
def connectivity():
    """API connectivity test endpoint"""
    return jsonify({
        "status": "ok",
        "connectivity": "established",
        "service": "main-api"
    })

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8080)