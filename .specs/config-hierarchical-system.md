# CryoProtect v2: Hierarchical Environment Configuration System

## Overview

This specification defines a robust, hierarchical configuration system for CryoProtect v2, supporting:
- A `BaseConfig` with environment-specific subclasses (`DevelopmentConfig`, `TestingConfig`, `StagingConfig`, `ProductionConfig`)
- Type validation, required/optional enforcement, and runtime error handling
- Fallback and override precedence (env vars, .env, Docker/CI/CD secrets, hardcoded defaults)
- Per-environment overrides and extensibility
- Automated and human-readable documentation requirements
- Compatibility with Docker, CI/CD, and security best practices

---

## 1. Class Hierarchy & Environment Selection

```mermaid
classDiagram
    class BaseConfig {
        +validate()
        +from_env()
        +as_dict()
        # _get_env_var(name, default, required, type)
    }
    class DevelopmentConfig
    class TestingConfig
    class StagingConfig
    class ProductionConfig

    BaseConfig <|-- DevelopmentConfig
    BaseConfig <|-- TestingConfig
    BaseConfig <|-- StagingConfig
    BaseConfig <|-- ProductionConfig
```

- **BaseConfig**: Defines all common configuration variables, their types, required/optional status, and default values.
- **Environment subclasses**: Override or extend variables as needed for each environment.
- **Selection**: The active config class is selected at runtime based on `FLASK_ENV` (or a similar variable), with a default fallback.

---

## 2. Type Validation & Error Handling

- **Type enforcement**: Each config variable must declare its type (e.g., str, int, bool, list, dict).
- **Required/optional**: Each variable is marked as required or optional. If required and missing, raise a clear error at startup.
- **Format validation**: For variables with specific formats (e.g., URLs, JSON, email), validate format and raise errors if invalid.
- **Error reporting**: On misconfiguration, raise a `ConfigurationError` with a helpful message (including variable name, expected type, and source).
- **Validation timing**: Validation occurs at application startup, before the app is fully initialized.

---

## 3. Fallback & Override Precedence

- **Precedence order** (highest to lowest):
    1. **Environment variables** (including those injected by Docker, CI/CD, or Kubernetes secrets)
    2. **.env file** (loaded via python-dotenv or similar)
    3. **Hardcoded defaults** (in BaseConfig or subclasses)
- **Per-environment overrides**: Environment-specific variables (e.g., `PRODUCTION_SUPABASE_URL`) override base variables if present.
- **Secret management**: For sensitive variables (e.g., `SECRET_KEY`, database credentials), support loading from Docker secrets or external secret managers (see secret management spec).

---

## 4. Per-Environment Overrides

- **Subclassing**: Each environment subclass can override any variable from BaseConfig.
- **Variable naming**: Use clear, environment-prefixed variable names for overrides (e.g., `PRODUCTION_SUPABASE_URL`).
- **Fallback**: If an environment-specific variable is missing, fall back to the base variable or raise an error if required.

---

## 5. Validation & Documentation Requirements

- **Validation**:
    - All required variables must be present and valid at startup.
    - Types and formats must be checked.
    - Provide a `validate()` method on each config class.
    - On error, print a summary of missing/invalid variables and exit.
- **Documentation**:
    - Each config variable must have a docstring or comment describing its purpose, type, required/optional status, and example value.
    - The `.env.template` must be kept in sync with config classes, grouping variables by function and marking required/optional.
    - Generate or maintain a human-readable configuration guide (see `docs/configuration_guide.md`).

---

## 6. Docker, CI/CD, and Security Compatibility

- **Docker**:
    - Support loading secrets from Docker secrets (e.g., `/run/secrets/SECRET_KEY`).
    - Do not hardcode secrets in images or code.
    - Validate configuration at container startup; fail fast on misconfiguration.
- **CI/CD**:
    - Allow environment variables to be injected at build/deploy time.
    - Support .env overrides for local development/testing.
    - Ensure validation works in CI pipelines (fail on missing/invalid config).
- **Security**:
    - Never commit real secrets to version control.
    - Mark all sensitive variables in documentation.
    - Provide warnings for insecure defaults (e.g., `SECRET_KEY`).

---

## 7. Example: Variable Definition Table

| Variable Name           | Type   | Required | Default | Description                                 | Example Value                  | Environment Override         |
|------------------------ |--------|----------|---------|---------------------------------------------|-------------------------------|-----------------------------|
| SUPABASE_URL            | str    | Yes      | None    | Supabase project URL                        | https://xyz.supabase.co        | STAGING_SUPABASE_URL, etc.  |
| SUPABASE_KEY            | str    | Yes      | None    | Supabase service role key                   | (long JWT string)              | PRODUCTION_SUPABASE_KEY     |
| SECRET_KEY              | str    | Yes      | dev-key | Flask secret key (must be strong in prod)   | (hex string)                   |                             |
| REDIS_URL               | str    | No       | None    | Redis connection URL                        | redis://user:pass@host:6379/0  |                             |
| LOG_LEVEL               | str    | No       | INFO    | Logging level                               | INFO, WARNING, ERROR           |                             |

---

## 8. Acceptance Criteria

- All configuration variables are validated for presence, type, and format at startup.
- Misconfiguration results in clear, actionable error messages and prevents app startup.
- Per-environment overrides are supported and documented.
- The system is compatible with Docker, CI/CD, and secret management best practices.
- The `.env.template` and documentation are kept in sync with config classes.
- Security best practices are enforced (no secrets in code, warnings for insecure defaults).

---

## 9. Implementation Notes

- Consider using `pydantic.BaseSettings`, `dataclasses`, or a custom metaclass for type enforcement and validation.
- Provide a utility to print all config variables and their sources (env, .env, default).
- Ensure backward compatibility with existing config usage patterns.
- Document all changes and update the configuration guide.

---