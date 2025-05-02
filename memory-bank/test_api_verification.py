import json
import requests

API_SPEC_PATH = "memory-bank/api_endpoints.json"
RESULTS_PATH = "memory-bank/verification_results.json"

AUTH_EMAIL = "eluecheelip@gmail.com"
AUTH_PASSWORD = "LDHt$rkaM&Gmf3X@LQ37"

def get_api_spec():
    with open(API_SPEC_PATH, "r", encoding="utf-8") as f:
        return json.load(f)

def get_auth_token(base_url):
    login_url = f"{base_url}/auth/login"
    data = {"email": AUTH_EMAIL, "password": AUTH_PASSWORD}
    try:
        resp = requests.post(login_url, json=data)
        resp.raise_for_status()
        token = resp.json().get("token")
        return token, resp.status_code, resp.json()
    except Exception as e:
        return None, getattr(e.response, "status_code", None), str(e)

def fill_path(path, example_id="550e8400-e29b-41d4-a716-446655440000"):
    # Replace {id} with example_id
    return path.replace("{id}", example_id)

def test_endpoint(base_url, endpoint, method, auth_token=None, query_params=None, body=None, expected_status=200, headers=None):
    url = base_url + endpoint
    if query_params:
        url += "?" + "&".join(f"{k}={v}" for k, v in query_params.items())
    req_headers = headers.copy() if headers else {}
    if auth_token:
        req_headers["Authorization"] = f"Bearer {auth_token}"
    try:
        if method == "GET":
            resp = requests.get(url, headers=req_headers)
        elif method == "POST":
            resp = requests.post(url, headers=req_headers, json=body)
        else:
            return {"status": "ERROR", "error": f"Unsupported method {method}"}
        try:
            resp_json = resp.json()
        except Exception:
            resp_json = resp.text
        return {
            "status_code": resp.status_code,
            "expected_status": expected_status,
            "response": resp_json,
            "success": resp.status_code == expected_status
        }
    except Exception as e:
        return {"status": "ERROR", "error": str(e)}

def main():
    spec = get_api_spec()
    base_url = spec.get("base_url", "http://localhost:5000")
    endpoints = spec["endpoints"]
    results = {}
    # Get auth token for authenticated endpoints
    token, auth_status, auth_resp = get_auth_token(base_url)
    results["auth"] = {
        "status_code": auth_status,
        "response": auth_resp,
        "success": bool(token)
    }
    for name, ep in endpoints.items():
        path = ep["path"]
        method = ep["method"]
        requires_auth = ep.get("authentication_required", False)
        # Fill in path parameters if needed
        if "{id}" in path:
            path = fill_path(path)
        # Prepare query parameters
        query_params = {}
        if "query_parameters" in ep:
            for k, v in ep["query_parameters"].items():
                query_params[k] = v.get("default", "")
        # Prepare request body for POST
        body = None
        if method == "POST" and "request_body" in ep:
            body = {}
            for k, v in ep["request_body"].items():
                if v.get("required", False):
                    # Use example or dummy value
                    if k == "cid":
                        body[k] = 702
                    elif k == "smiles":
                        body[k] = "CCO"
                    elif k == "name":
                        body[k] = "Ethanol"
                    else:
                        body[k] = "test"
                else:
                    # Optional fields: skip or use example
                    pass
        # Set expected status
        expected_status = 200
        # Run test
        result = test_endpoint(
            base_url,
            path,
            method,
            auth_token=token if requires_auth else None,
            query_params=query_params,
            body=body,
            expected_status=expected_status
        )
        results[name] = result
    # Write results
    with open(RESULTS_PATH, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

if __name__ == "__main__":
    main()