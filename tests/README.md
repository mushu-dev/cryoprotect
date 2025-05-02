# Tests Directory

This README provides an overview of the test structure, types of tests, instructions for running tests, and guidelines for contributing new tests to the project.

---

## Table of Contents

- [Introduction](#introduction)
- [Test Directory Structure](#test-directory-structure)
- [Test Types](#test-types)
- [How to Run Tests](#how-to-run-tests)
- [Test Utilities](#test-utilities)
- [Adding New Tests](#adding-new-tests)
- [Troubleshooting](#troubleshooting)
- [References](#references)

---

## Introduction

The `tests/` directory contains all automated tests for the CryoProtect project. These tests ensure the correctness, reliability, and maintainability of the codebase. All contributors are encouraged to run the tests regularly and add new tests when submitting code changes.

---

## Test Directory Structure

```
tests/
├── unit/                # Unit tests for individual modules/functions
├── integration/         # Integration tests for combined components
├── e2e/                 # End-to-end tests simulating real user workflows
├── utils/               # Shared test utilities and fixtures
└── README.md            # This documentation file
```

---

## Test Types

- **Unit Tests**: Test individual functions or classes in isolation.
- **Integration Tests**: Test interactions between multiple modules or components.
- **End-to-End (E2E) Tests**: Simulate real-world scenarios and user workflows.
- **Utility Modules**: Provide reusable helpers, fixtures, and mock data for tests.

---

## How to Run Tests

### Prerequisites

- Ensure all dependencies are installed:
  ```sh
  pip install -r requirements.txt
  ```

### Running All Tests

From the project root directory, run:
```sh
python -m unittest discover tests
```
or, if using `pytest`:
```sh
pytest tests
```

### Running Specific Test Modules

```sh
python -m unittest tests.unit.test_example
```

---

## Test Utilities

Common utilities and fixtures are located in `tests/utils/`. These modules help reduce code duplication and simplify test setup.

---

## Adding New Tests

1. Place new unit tests in `tests/unit/`, integration tests in `tests/integration/`, and E2E tests in `tests/e2e/`.
2. Name test files with the `test_*.py` pattern.
3. Use descriptive test function names and docstrings.
4. Reuse utilities from `tests/utils/` where possible.
5. Run all tests before submitting a pull request.

---

## Troubleshooting

- Ensure your environment matches the project requirements.
- Check for missing dependencies or incorrect paths.
- Review error messages for hints on failing tests.
- For environment-specific issues, see the main project README or contact the maintainers.

---

## References

- [Python unittest documentation](https://docs.python.org/3/library/unittest.html)
- [pytest documentation](https://docs.pytest.org/en/stable/)
- [Project Contribution Guide](../README_Developer_Guide.md)