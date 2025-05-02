# CryoProtect Analyzer - Authentication System

This document provides an overview of the authentication system implemented in the CryoProtect Analyzer application.

## Overview

The CryoProtect Analyzer uses Supabase for authentication, providing a complete authentication flow that allows users to:

1. Register new accounts
2. Log in with existing accounts
3. Log out
4. Reset passwords
5. View and edit their profile information

The authentication system is integrated with the Flask backend and Supabase database, with a frontend built using Bootstrap 5.

## Architecture

### Backend Components

- **Flask Routes**: Web routes for authentication pages (`/login`, `/register`, `/profile`, `/reset-password`)
- **API Endpoints**: REST endpoints for authentication operations (`/auth/login`, `/auth/register`, etc.)
- **Authentication Middleware**: Protects routes that require authentication
- **Supabase Integration**: Uses the Supabase Python client for authentication operations

### Frontend Components

- **HTML Templates**: Login, registration, profile, and password reset pages
- **Authentication Module (auth.js)**: Core authentication functionality
- **UI Module (auth-ui.js)**: Handles UI updates based on authentication state
- **Supabase JS Client**: Communicates with Supabase for authentication operations

## Authentication Flow

### Registration

1. User navigates to `/register`
2. User fills out the registration form with email and password
3. Frontend sends registration request to `/auth/register`
4. Backend creates a new user in Supabase
5. User is redirected to the login page or automatically logged in

### Login

1. User navigates to `/login`
2. User enters email and password
3. Frontend authenticates with Supabase
4. On successful authentication:
   - Authentication token is stored in localStorage
   - User data is stored in localStorage
   - User is redirected to the dashboard

### Logout

1. User clicks the logout button
2. Frontend sends logout request to Supabase
3. Authentication token and user data are cleared from localStorage
4. User is redirected to the login page

### Password Reset

1. User navigates to `/reset-password`
2. User enters their email address
3. Backend sends a password reset link via Supabase
4. User receives an email with a reset link
5. User clicks the link and is directed to the password reset page
6. User enters a new password
7. Password is updated in Supabase

### Profile Management

1. User navigates to `/profile`
2. User can view and edit their profile information
3. User can change their password
4. Changes are sent to the backend and updated in Supabase

## Route Protection

Routes that require authentication are protected by middleware that checks if the user is authenticated. If not, the user is redirected to the login page.

Protected routes include:
- `/profile`
- `/molecules`
- `/mixtures`
- `/predictions`
- `/experiments`
- `/comparisons`

## Session Management

- Authentication tokens are stored in localStorage
- User data is stored in localStorage
- The application checks for valid sessions on page load
- Sessions are verified with Supabase to ensure they are still valid

## Security Considerations

- Passwords are never stored in plaintext
- Authentication tokens are stored securely in localStorage
- Password reset links are time-limited
- Supabase handles password hashing and token generation
- Row Level Security (RLS) policies in Supabase control data access

## Implementation Details

### Backend (Flask)

The authentication system is implemented in the following files:
- `app.py`: Contains routes and API endpoints
- `api/utils.py`: Contains authentication utilities

### Frontend

The authentication system is implemented in the following files:
- `static/js/auth.js`: Core authentication functionality
- `static/js/auth-ui.js`: UI updates based on authentication state
- `templates/login.html`: Login page
- `templates/register.html`: Registration page
- `templates/profile.html`: Profile management page
- `templates/reset_password.html`: Password reset page

## Configuration

Authentication configuration is stored in the `.env` file:

```
SUPABASE_URL=your-supabase-url
SUPABASE_KEY=your-supabase-key
```

These values are loaded into the Flask application via `config.py` and passed to the frontend templates.

## Troubleshooting

### Common Issues

1. **Authentication Failed**: Check that the Supabase URL and key are correct in the `.env` file.
2. **Password Reset Email Not Received**: Check spam folder or verify the email address is correct.
3. **Session Expired**: The user may need to log in again if their session has expired.

### Debugging

- Check browser console for JavaScript errors
- Check Flask logs for backend errors
- Verify Supabase configuration in the `.env` file