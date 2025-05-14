#!/usr/bin/env python3
"""
Simplified Flask application for CryoProtect to test Heroku deployment.
"""

import os
from flask import Flask, jsonify

# Create the app
app = Flask(__name__)

# Add a simple health check endpoint
@app.route('/')
def index():
    return jsonify({
        'status': 'ok', 
        'message': 'CryoProtect API is running', 
        'version': '1.0.0'
    })

@app.route('/health')
def health_check():
    return jsonify({
        'status': 'ok',
        'version': '1.0.0',
        'message': 'Health check passed'
    })

if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port)