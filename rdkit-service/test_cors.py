"""
Simple Flask app to test CORS configuration for RDKit service
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
        "rdkit": "available"
    })

@app.route('/test-cors', methods=['GET', 'OPTIONS'])
def cors_test():
    """Endpoint for testing CORS configuration"""
    return jsonify({
        "success": True,
        "message": "CORS is configured properly",
        "service": "rdkit-service"
    })

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=8080)