# Supabase API Integration Test Report

**Date:** 2025-04-18 20:52:39

## Summary

- **Total Tests:** 19
- **Passed:** 1
- **Failed:** 0
- **Errors:** 18
- **Success Rate:** 5.26%

## Errors

### 1. test_admin_user_can_add_prediction_to_any_mixture (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_admin_user_can_add_prediction_to_any_mixture)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 585, in test_admin_user_can_add_prediction_to_any_mixture
    self.mock_auth_user(mock_client, self.admin_user_id, self.admin_auth_id, is_admin=True)
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 172, in mock_auth_user
    mock_client.auth.get_user.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

### 2. test_admin_user_can_update_any_mixture (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_admin_user_can_update_any_mixture)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 539, in test_admin_user_can_update_any_mixture
    self.mock_auth_user(mock_client, self.admin_user_id, self.admin_auth_id, is_admin=True)
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 172, in mock_auth_user
    mock_client.auth.get_user.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

### 3. test_authenticated_user_add_experiment (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_authenticated_user_add_experiment)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 408, in test_authenticated_user_add_experiment
    self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 172, in mock_auth_user
    mock_client.auth.get_user.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

### 4. test_authenticated_user_add_prediction (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_authenticated_user_add_prediction)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 364, in test_authenticated_user_add_prediction
    self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 172, in mock_auth_user
    mock_client.auth.get_user.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

### 5. test_authenticated_user_cannot_add_experiment_to_others_mixture (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_authenticated_user_cannot_add_experiment_to_others_mixture)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 511, in test_authenticated_user_cannot_add_experiment_to_others_mixture
    self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 172, in mock_auth_user
    mock_client.auth.get_user.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

### 6. test_authenticated_user_cannot_add_prediction_to_others_mixture (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_authenticated_user_cannot_add_prediction_to_others_mixture)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 485, in test_authenticated_user_cannot_add_prediction_to_others_mixture
    self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 172, in mock_auth_user
    mock_client.auth.get_user.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

### 7. test_authenticated_user_cannot_update_others_mixture (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_authenticated_user_cannot_update_others_mixture)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 454, in test_authenticated_user_cannot_update_others_mixture
    self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 172, in mock_auth_user
    mock_client.auth.get_user.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

### 8. test_authenticated_user_create_mixture (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_authenticated_user_create_mixture)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 263, in test_authenticated_user_create_mixture
    self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 172, in mock_auth_user
    mock_client.auth.get_user.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

### 9. test_authenticated_user_update_own_mixture (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_authenticated_user_update_own_mixture)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 318, in test_authenticated_user_update_own_mixture
    self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 172, in mock_auth_user
    mock_client.auth.get_user.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

### 10. test_not_found_error (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_not_found_error)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 734, in test_not_found_error
    mock_client.table('molecule_with_properties').select().eq().execute.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TypeError: MockQueryBuilder.eq() missing 2 required positional arguments: 'column' and 'value'

```

### 11. test_search_molecules_by_name (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_search_molecules_by_name)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 714, in test_search_molecules_by_name
    mock_client.table('molecules').select().ilike().limit().execute.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TypeError: MockQueryBuilder.ilike() missing 2 required positional arguments: 'column' and 'value'

```

### 12. test_server_error (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_server_error)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 772, in test_server_error
    mock_client.table('molecule_with_properties').select().eq().execute.side_effect = Exception("Database error")
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TypeError: MockQueryBuilder.eq() missing 2 required positional arguments: 'column' and 'value'

```

### 13. test_service_role_can_access_any_data (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_service_role_can_access_any_data)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1387, in patched
    with self.decoration_helper(patched,
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\contextlib.py", line 137, in __enter__
    return next(self.gen)
           ^^^^^^^^^^^^^^
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1369, in decoration_helper
    arg = exit_stack.enter_context(patching)
          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\contextlib.py", line 526, in enter_context
    result = _enter(cm)
             ^^^^^^^^^^
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1458, in __enter__
    original, local = self.get_original()
                      ^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1431, in get_original
    raise AttributeError(
AttributeError: <module 'api.utils' from 'C:\\Users\\1edwa\\Documents\\CryoProtect v2\\api\\utils.py'> does not have the attribute 'USER_ID'

```

### 14. test_unauthenticated_cannot_create_mixture (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_unauthenticated_cannot_create_mixture)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 239, in test_unauthenticated_cannot_create_mixture
    response = self.client.post(
               ^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\werkzeug\test.py", line 1167, in post
    return self.open(*args, **kw)
           ^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask\testing.py", line 235, in open
    response = super().open(
               ^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\werkzeug\test.py", line 1116, in open
    response_parts = self.run_wsgi_app(request.environ, buffered=buffered)
                     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\werkzeug\test.py", line 988, in run_wsgi_app
    rv = run_wsgi_app(self.application, environ, buffered=buffered)
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\werkzeug\test.py", line 1264, in run_wsgi_app
    app_rv = app(environ, start_response)
             ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask\app.py", line 1488, in __call__
    return self.wsgi_app(environ, start_response)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask\app.py", line 1466, in wsgi_app
    response = self.handle_exception(e)
               ^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask_restful\__init__.py", line 298, in error_router
    return original_handler(e)
           ^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask_cors\extension.py", line 176, in wrapped_function
    return cors_after_request(app.make_response(f(*args, **kwargs)))
                                                ^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask_restful\__init__.py", line 295, in error_router
    return self.handle_error(e)
           ^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask_restful\__init__.py", line 310, in handle_error
    _handle_flask_propagate_exceptions_config(current_app, e)
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask\app.py", line 1463, in wsgi_app
    response = self.full_dispatch_request()
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask\app.py", line 872, in full_dispatch_request
    rv = self.handle_user_exception(e)
         ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask_restful\__init__.py", line 298, in error_router
    return original_handler(e)
           ^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask_cors\extension.py", line 176, in wrapped_function
    return cors_after_request(app.make_response(f(*args, **kwargs)))
                                                ^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask_restful\__init__.py", line 295, in error_router
    return self.handle_error(e)
           ^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask_restful\__init__.py", line 310, in handle_error
    _handle_flask_propagate_exceptions_config(current_app, e)
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask\app.py", line 870, in full_dispatch_request
    rv = self.dispatch_request()
         ^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask\app.py", line 855, in dispatch_request
    return self.ensure_sync(self.view_functions[rule.endpoint])(**view_args)  # type: ignore[no-any-return]
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask_restful\__init__.py", line 493, in wrapper
    return self.make_response(data, code, headers=headers)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask_restful\__init__.py", line 522, in make_response
    resp = self.representations[mediatype](data, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\.venv\Lib\site-packages\flask_restful\representations\json.py", line 21, in output_json
    dumped = dumps(data, **settings) + "\n"
             ^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\json\__init__.py", line 238, in dumps
    **kw).encode(obj)
          ^^^^^^^^^^^
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\json\encoder.py", line 202, in encode
    chunks = list(chunks)
             ^^^^^^^^^^^^
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\json\encoder.py", line 439, in _iterencode
    o = _default(o)
        ^^^^^^^^^^^
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\json\encoder.py", line 180, in default
    raise TypeError(f'Object of type {o.__class__.__name__} '
TypeError: Object of type Response is not JSON serializable

```

### 15. test_unauthenticated_get_mixtures (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_unauthenticated_get_mixtures)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 222, in test_unauthenticated_get_mixtures
    mock_client.table('mixture_with_components').select('*').execute.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

### 16. test_unauthenticated_get_molecule_detail (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_unauthenticated_get_molecule_detail)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 206, in test_unauthenticated_get_molecule_detail
    mock_client.table('molecule_with_properties').select('*').eq().execute.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
TypeError: MockQueryBuilder.eq() missing 2 required positional arguments: 'column' and 'value'

```

### 17. test_unauthenticated_get_molecules (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_unauthenticated_get_molecules)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 188, in test_unauthenticated_get_molecules
    mock_client.table('molecules').select('*').execute.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

### 18. test_validation_error (tests.test_supabase_api_integration.TestSupabaseAPIIntegration.test_validation_error)

```
Traceback (most recent call last):
  File "C:\Users\1edwa\AppData\Local\Programs\Python\Python312\Lib\unittest\mock.py", line 1390, in patched
    return func(*newargs, **newkeywargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\mock_supabase\helpers.py", line 35, in wrapper
    return test_method(self, mock_client, *args, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 750, in test_validation_error
    self.mock_auth_user(mock_client, self.regular_user_id, self.regular_auth_id)
  File "C:\Users\1edwa\Documents\CryoProtect v2\tests\test_supabase_api_integration.py", line 172, in mock_auth_user
    mock_client.auth.get_user.return_value = MagicMock(
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'method' object has no attribute 'return_value'

```

## Authentication Test Results

