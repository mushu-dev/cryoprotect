# ROO TASK 1.3: Database Health Check Utilities Implementation Prompt

## CONTEXT
Tasks 1.1 (Foreign Key Relationship Fixes) and 1.2 (RLS Implementation Tools) have been successfully completed. We now need to implement Task 1.3: Database Health Check Utilities as part of the Essential Maintenance Utilities component of Phase 2. These utilities will ensure database integrity, performance, and schema validity, providing a foundation for monitoring database health in production.

## OBJECTIVE
Create a comprehensive set of database health check utilities that can:
1. Validate database schema against expected structure
2. Verify data integrity and relationships
3. Identify performance issues and bottlenecks
4. Generate detailed health reports
5. Provide automated notifications for critical issues

## FILES TO CREATE/MODIFY

### Core Implementation
1. `/database/utils/health_check.py` - Main health check implementation with core functions
2. `/database/utils/schema_validator.py` - Schema validation utilities
3. `/database/utils/integrity_checker.py` - Data integrity verification
4. `/database/utils/performance_analyzer.py` - Performance testing and analysis

### API Integration
5. `/api/system_resources.py` - Update to add health check API endpoints

### Testing
6. `/tests/test_database_health_check.py` - Tests for health check utilities

### Documentation
7. `/docs/database_health_check.md` - Comprehensive documentation

## IMPLEMENTATION REQUIREMENTS

### Database Schema Validation
- Verify all required tables exist
- Check column definitions and types
- Validate constraints (primary keys, foreign keys, unique)
- Ensure views are properly defined
- Verify RLS policies are in place

### Data Integrity Verification
- Check for orphaned records violating foreign key constraints
- Verify unique constraints are respected
- Check for data consistency across related tables
- Validate required fields are properly populated
- Identify potential data corruption

### Performance Analysis
- Measure query execution times for common operations
- Check index usage and effectiveness
- Identify slow queries and bottlenecks
- Monitor connection pool utilization
- Analyze table statistics and growth patterns

### Health Reporting
- Generate comprehensive health reports
- Support different output formats (JSON, Markdown, HTML)
- Include severity levels for issues
- Provide remediation recommendations
- Track historical health metrics

### API Integration
- Create `/health` endpoint for basic uptime monitoring
- Add `/health/database` endpoint for database checks
- Implement `/health/performance` for performance metrics
- Include authentication for detailed health information
- Support different verbosity levels

## CODE PATTERNS AND EXAMPLES

### Health Check Module Structure
```python
class DatabaseHealthCheck:
    """Database health check implementation."""
    
    def __init__(self, connection_pool=None, config=None):
        """Initialize health check with connection pool and configuration."""
        self.connection_pool = connection_pool
        self.config = config or {}
        self.validators = []
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "status": "unknown",
            "checks": {}
        }
        
        # Register validators
        self._register_validators()
    
    def _register_validators(self):
        """Register all available validators."""
        self.validators = [
            SchemaValidator(self.config.get("schema", {})),
            IntegrityChecker(self.config.get("integrity", {})),
            PerformanceAnalyzer(self.config.get("performance", {}))
        ]
    
    def run_health_check(self, categories=None):
        """Run health check for specified categories or all."""
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "status": "pending",
            "checks": {}
        }
        
        try:
            for validator in self.validators:
                if categories and validator.category not in categories:
                    continue
                
                results = validator.run_checks(self.connection_pool)
                self.results["checks"][validator.category] = results
            
            # Determine overall status
            self._calculate_status()
            
            return self.results
        except Exception as e:
            self.results["status"] = "error"
            self.results["error"] = str(e)
            return self.results
    
    def _calculate_status(self):
        """Calculate overall health status based on check results."""
        has_critical = False
        has_warning = False
        
        for category, results in self.results["checks"].items():
            if results.get("status") == "critical":
                has_critical = True
            elif results.get("status") == "warning":
                has_warning = True
        
        if has_critical:
            self.results["status"] = "critical"
        elif has_warning:
            self.results["status"] = "warning"
        else:
            self.results["status"] = "healthy"
    
    def generate_report(self, format="json", output_file=None):
        """Generate health check report in the specified format."""
        if format == "json":
            report = json.dumps(self.results, indent=2)
        elif format == "markdown":
            report = self._generate_markdown_report()
        elif format == "html":
            report = self._generate_html_report()
        else:
            raise ValueError(f"Unsupported format: {format}")
        
        if output_file:
            with open(output_file, "w") as f:
                f.write(report)
        
        return report
```

### API Integration Example
```python
@app.route('/health', methods=['GET'])
def health_check():
    """Basic health check endpoint."""
    return jsonify({
        "status": "ok",
        "timestamp": datetime.now().isoformat(),
        "version": app.config.get('VERSION', 'unknown')
    })

@app.route('/health/database', methods=['GET'])
@jwt_required
def database_health_check():
    """Detailed database health check endpoint (authenticated)."""
    # Check if user has admin role
    if not current_user.has_role('admin'):
        return jsonify({
            "status": "ok",
            "timestamp": datetime.now().isoformat(),
            "message": "Authenticated, but detailed health check requires admin role"
        }), 403
    
    # Run health check
    health_checker = DatabaseHealthCheck(app.db_pool, app.config.get('HEALTH_CHECK', {}))
    results = health_checker.run_health_check()
    
    # Return results with appropriate status code
    status_code = 200
    if results["status"] == "critical":
        status_code = 500
    elif results["status"] == "warning":
        status_code = 299  # Custom code for warning
    
    return jsonify(results), status_code
```

## IMPLEMENTATION APPROACH
1. First develop the core health check utilities
2. Create comprehensive tests for each component
3. Integrate with API endpoints
4. Generate documentation and example reports
5. Validate against database inconsistencies

## SUCCESS CRITERIA
- All health check utilities are implemented and working correctly
- Tests demonstrate proper identification of database issues
- API endpoints return correct data and status codes
- Documentation is clear and includes examples
- The solution can detect common database problems

## DELIVERABLES
1. Complete implementation of health check modules
2. API endpoints for health monitoring
3. Comprehensive test suite
4. Documentation and sample reports
5. Integration with existing monitoring system

This implementation will provide a solid foundation for maintaining database health in the CryoProtect application, ensuring data integrity and optimal performance.