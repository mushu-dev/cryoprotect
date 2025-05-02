# Project Manager Directive: Complete Testing Framework Unification

## Context
We're working on Task 3 from our project plan: "Unify testing framework" to consolidate test files into an organized structure with shared fixtures, improved coverage reporting, and standardized test patterns. We've already completed two core components (Database Fixtures and Mock Objects), and now need to implement the remaining components.

## Directive
As the Project Manager agent, your task is to coordinate the implementation of our Testing Framework Unification plan. The detailed plan is available at `.roo/testing_framework_plan.md`. Your responsibilities include:

1. Review the current status of our testing framework implementation
2. Delegate specialized tasks to the appropriate agents
3. Monitor progress and ensure components integrate properly
4. Validate completed work against our acceptance criteria
5. Coordinate the sequencing of dependent tasks

## Delegations

### For the QA Agent
You should immediately delegate the API Fixtures implementation to the QA agent with the following prompt:

```
I need you to implement the API Fixtures component for our testing framework unification project. 

This implementation should:
1. Follow the detailed plan in API_FIXTURES.md
2. Create the directory structure at /tests/fixtures/api/
3. Implement client.py, auth.py, and __init__.py with the specified fixtures
4. Create example tests in /tests/unit/api/test_api_fixtures.py
5. Add documentation in README.md

Your implementation must integrate with our existing Flask app structure and authentication methods. Upon completion, the API fixtures should allow for easy testing of all API endpoints with proper authentication and request/response handling.

Key requirements:
- Support both authenticated and unauthenticated testing
- Include fixtures for different user roles (admin, regular user, scientist)
- Provide JSON request/response handling
- Include comprehensive documentation with usage examples

The Database Fixtures and Mock Objects components are already implemented and can be used as reference for style and integration patterns.
```

### For the Data Scientist Agent
After the API Fixtures are implemented, delegate the Test Data Fixtures implementation to the Data Scientist agent with:

```
I need you to implement the Test Data Fixtures component for our testing framework unification project.

This implementation should:
1. Follow the detailed plan in TEST_DATA_FIXTURES.md
2. Create the directory structure at /tests/fixtures/data/
3. Implement data generation utilities for all entity types
4. Create standard test data sets for molecules, mixtures, experiments, etc.
5. Add comprehensive documentation with usage examples

Your implementation should focus on:
- Creating realistic test data for all scientific entity types
- Providing flexible data generation utilities
- Supporting both simple and complex test scenarios
- Ensuring consistent data relationships
- Including documentation of the test data schema

The Database Fixtures, Mock Objects, and API Fixtures components are already implemented and can be referenced for integration patterns.
```

### For Final Integration (QA Agent)
After both components are complete, delegate the Conftest Update to the QA agent:

```
I need you to implement the Conftest Update component for our testing framework unification project.

This implementation should:
1. Follow the detailed plan in CONFTEST_UPDATE.md
2. Update /tests/conftest.py to import and expose all fixtures
3. Add pytest configuration options
4. Ensure proper fixture isolation and scope
5. Implement test discovery improvements
6. Add documentation for pytest options

This is the final integration step that brings together all the testing framework components:
- Database Fixtures
- Mock Objects
- API Fixtures
- Test Data Fixtures

Your update should ensure all components work together seamlessly while maintaining proper isolation.
```

## Monitoring and Coordination
After delegating each task:
1. Monitor progress and provide guidance when needed
2. Check that implementations follow our standards and integrate properly
3. Coordinate dependencies between components
4. Validate that all tests pass after each implementation
5. Ensure comprehensive documentation is maintained

## Verification
Before considering the Testing Framework Unification complete, verify:
1. All example tests pass
2. The test suite runs without errors
3. Coverage reports are generated correctly
4. All fixtures are properly exposed and documented
5. The test patterns are consistently applied

## Reporting
Provide regular status updates on:
1. Current progress of each component
2. Any blockers or issues encountered
3. Next steps in the implementation plan
4. Overall timeline for completion

## References
- Testing Framework Plan: `.roo/testing_framework_plan.md`
- API Fixtures Plan: `API_FIXTURES.md`
- Test Data Fixtures Plan: `TEST_DATA_FIXTURES.md`
- Conftest Update Plan: `CONFTEST_UPDATE.md`