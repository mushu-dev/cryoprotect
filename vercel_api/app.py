from http.server import BaseHTTPRequestHandler
import sys
import os

# Add the parent directory to sys.path to import app.py
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from app import app
except ImportError as e:
    print(f"Error importing app: {e}")

def handler(request):
    # Get the Flask app response
    from io import BytesIO
    import urllib.parse
    
    # Extract path and method
    path = request.get('path', '/')
    method = request.get('method', 'GET')
    
    # Query string parameters
    query = request.get('query', {})
    query_string = urllib.parse.urlencode(query)
    
    # Prepare the WSGI environment
    env = {
        "wsgi.input": BytesIO(),
        "wsgi.errors": sys.stderr,
        "wsgi.version": (1, 0),
        "wsgi.multithread": False,
        "wsgi.multiprocess": False,
        "wsgi.run_once": False,
        "wsgi.url_scheme": "https",
        "REQUEST_METHOD": method,
        "PATH_INFO": path,
        "QUERY_STRING": query_string,
        "SERVER_PROTOCOL": "HTTP/1.1",
        "SERVER_NAME": "vercel-serverless",
        "SERVER_PORT": "443",
        "REMOTE_ADDR": request.get('headers', {}).get('x-forwarded-for', '127.0.0.1'),
    }
    
    # Add HTTP headers to the environment
    for name, value in request.get('headers', {}).items():
        name = name.upper().replace("-", "_")
        if name not in ("CONTENT_TYPE", "CONTENT_LENGTH"):
            name = "HTTP_" + name
        env[name] = value
    
    # Response data
    response_body = []
    status = "200 OK"
    headers = []
    
    def start_response(status_line, response_headers):
        nonlocal status, headers
        status = status_line
        headers = response_headers
    
    # Call the Flask app
    try:
        result = app(env, start_response)
        for data in result:
            if data:
                response_body.append(data)
    except Exception as e:
        print(f"Error calling app: {e}")
        status = "500 Internal Server Error"
        response_body = [str(e).encode('utf-8')]
    
    # Parse the status code
    status_code = int(status.split(' ')[0])
    
    # Convert headers to a dictionary
    headers_dict = {}
    for key, value in headers:
        headers_dict[key] = value
    
    # Return the response in Vercel serverless format
    return {
        'statusCode': status_code,
        'headers': headers_dict,
        'body': b''.join(response_body).decode('utf-8', errors='replace')
    }