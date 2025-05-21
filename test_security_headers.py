import requests

def test_security_headers(base_url="http://127.0.0.1:5000/health"):
    required_headers = [
        "Content-Security-Policy",
        "Strict-Transport-Security",
        "X-Content-Type-Options",
        "X-Frame-Options",
        "X-XSS-Protection",
        "Referrer-Policy",
        "Feature-Policy"
    ]
    response = requests.get(base_url)
    print(f"Status code: {response.status_code}")
    missing = []
    for header in required_headers:
        if header in response.headers:
            print(f"{header}: {response.headers[header]}")
        else:
            print(f"Missing header: {header}")
            missing.append(header)
    if missing:
        print("Test failed: Missing headers:", missing)
    else:
        print("All security headers are present. Test passed.")

if __name__ == "__main__":
    test_security_headers()