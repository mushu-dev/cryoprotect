# CryoProtect v2 Updated Master Plan

## Current Status Overview
The CryoProtect v2 project has made significant progress (approximately 75% complete) across several critical areas:

1. ✅ **Database Module Consolidation**: Complete with centralized operations, unified connection management, and migration scripts.

2. ✅ **Testing Framework**: Near completion with:
   - Database fixtures fully implemented
   - Mock objects implementation complete
   - API fixtures testing finalized (excluding health endpoint)
   - Test Data Fixtures implementation in progress

3. 🔄 **API Architecture Standardization**: In progress with:
   - `api/dashboard_resources.py` standardized with consistent error handling, authentication, and response formatting
   - Documentation started in `README_API_FIXES.md`
   - Full audit of remaining API resources completed

4. ⚠️ **Missing Critical Functionality**: Lab verification workflow for validating predictions against actual lab results.

## Revised Priority Tasks for Project Completion

### 1. Implement Lab Verification Workflow (HIGHEST PRIORITY)
**Progress**: 0% complete (not started)

**Next Steps**:
- Design and implement the `LabVerification` class in `api/models.py`
- Add verification status tracking (pending, verified, rejected)
- Create API endpoints for recording and retrieving verification data
- Enhance comparison functionality to include verification information
- Implement basic reporting for verification analysis

**Standards to Apply**:
- Consistent error handling (`handle_error`)
- Authentication enforcement (`@token_required`)
- Response formatting (`@marshal_with` with field schemas)
- Request validation
- Comprehensive docstrings

**Delegation**: Backend Agent with review from Data Scientist Agent

### 2. Complete API Architecture Standardization (HIGH PRIORITY)
**Progress**: 10% complete (1 of ~10 resource files standardized)

**Next Steps**:
- Standardize resource files in priority order: 
  1. `api/resources.py` (core functionality)
  2. `api/mixture_analysis_resources.py` (high usage)
  3. `api/rdkit_resources.py` (critical functionality)
  4. Remaining resource files in usage frequency order

**Delegation**: Backend Agent with PR reviews from QA Agent

### 3. Finish Testing Framework (HIGH PRIORITY)
**Progress**: 75% complete

**Next Steps**:
- Complete Test Data Fixtures implementation (in progress)
- Implement Conftest Update to expose all fixtures
- Run full test suite and verify coverage
- Address any remaining test failures or warnings

**Delegation**: QA Agent with review from Data Scientist Agent for data fixtures

### 4. Implement Essential Maintenance Utilities (MEDIUM PRIORITY)
**Progress**: 30% complete

**Next Steps**:
- Complete foreign key relationship fixes
- Finish RLS implementation tools
- Add database health check utilities
- Implement basic backup/restore tools

**Delegation**: DBA Agent

### 5. Complete Predictive Models (MEDIUM PRIORITY)
**Progress**: 25% complete

**Next Steps**:
- Implement core prediction algorithms
- Create validation framework
- Add visualization components
- Implement model training/evaluation

**Delegation**: Data Scientist Agent

### 6. Minimal Documentation (MEDIUM PRIORITY)
**Progress**: 40% complete

**Next Steps**:
- Standardize README formats
- Generate API documentation
- Create essential user documentation
- Develop basic deployment guides

**Delegation**: PM Agent with input from all other agents

### 7. Minimal Viable Deployment (LOW PRIORITY)
**Progress**: 20% complete

**Next Steps**:
- Implement basic CI/CD pipeline
- Configure essential environment variables
- Create deployment verification
- Set up minimal monitoring

**Delegation**: DevOps Agent

### 8. Deferred to Post-Launch
- Enhanced RDKit Integration optimizations
- Advanced monitoring features
- Blue/green deployment
- Advanced security features
- Knowledge transfer materials
- Video tutorials

## Revised Implementation Timeline

### Phase 1: Essential Core Functionality (4 Weeks)
- Implement Lab Verification Workflow (1.5 weeks)
- Complete API Architecture Standardization (1 week)
- Finish Testing Framework (1 week)
- Fix critical bugs (0.5 weeks)

### Phase 2: Minimal Viable Deployment (3 Weeks)
- Implement Essential Maintenance Utilities (1 week)
- Complete Predictive Models (1 week)
- Set up basic CI/CD pipeline (0.5 weeks)
- Configure minimal monitoring (0.5 weeks)

### Phase 3: Essential Documentation (2 Weeks)
- Standardize README formats (0.5 weeks)
- Generate API documentation (0.5 weeks)
- Create essential user documentation (0.5 weeks)
- Create basic deployment guides (0.5 weeks)

### Phase 4: Post-Launch Enhancements (Deferred)
- Enhanced RDKit Integration
- Advanced monitoring
- Blue/green deployment
- Advanced security features
- Knowledge transfer

## Success Metrics

1. **Core Functionality**:
   - Lab verification workflow fully implemented and tested
   - All APIs standardized with consistent patterns
   - 80%+ test coverage for critical components
   - Predictive models functional and validated

2. **Deployment Readiness**:
   - Basic CI/CD pipeline operational
   - Environment configuration documented
   - Essential monitoring in place
   - Backup/restore procedures tested

3. **Documentation**:
   - API endpoints fully documented
   - User documentation covers essential workflows
   - Deployment process clearly documented
   - Database schema documented

## Detailed Lab Verification Implementation Plan

### Data Model (File: api/models.py)
```python
class LabVerification(BaseModel):
    """Model for lab verification data."""
    
    table_name = 'lab_verifications'
    
    # Verification states
    PENDING = 'pending'
    VERIFIED = 'verified'
    REJECTED = 'rejected'
    
    @classmethod
    def record_verification(cls, experiment_id: str, 
                          verification_status: str,
                          verifier: str, 
                          equipment_used: str,
                          comments: str = None) -> Dict[str, Any]:
        """Record verification for an experiment."""
        # Implementation details here
        
    @classmethod
    def get_verification(cls, experiment_id: str) -> Dict[str, Any]:
        """Get verification for an experiment."""
        # Implementation details here
        
    @classmethod
    def update_verification_status(cls, verification_id: str, 
                                  new_status: str,
                                  comments: str = None) -> Dict[str, Any]:
        """Update verification status."""
        # Implementation details here
```

### Database Migration (File: migrations/013_lab_verification_schema.sql)
```sql
-- Create lab_verifications table
CREATE TABLE IF NOT EXISTS public.lab_verifications (
    id UUID PRIMARY KEY DEFAULT gen_random_uuid(),
    experiment_id UUID NOT NULL REFERENCES public.experiments(id) ON DELETE CASCADE,
    verification_status VARCHAR(20) NOT NULL CHECK (verification_status IN ('pending', 'verified', 'rejected')),
    verifier VARCHAR(255) NOT NULL,
    equipment_used VARCHAR(255) NOT NULL,
    comments TEXT,
    created_at TIMESTAMP WITH TIME ZONE DEFAULT NOW(),
    updated_at TIMESTAMP WITH TIME ZONE DEFAULT NOW()
);

-- Add RLS policies
ALTER TABLE public.lab_verifications ENABLE ROW LEVEL SECURITY;

CREATE POLICY lab_verifications_select ON public.lab_verifications
    FOR SELECT USING (auth.uid() = verifier OR auth.uid() IN (
        SELECT user_id FROM user_team_membership utm
        JOIN team_roles tr ON utm.role_id = tr.id
        WHERE tr.can_view_verifications = true
    ));

CREATE POLICY lab_verifications_insert ON public.lab_verifications
    FOR INSERT WITH CHECK (auth.uid() = verifier);

CREATE POLICY lab_verifications_update ON public.lab_verifications
    FOR UPDATE USING (auth.uid() = verifier OR auth.uid() IN (
        SELECT user_id FROM user_team_membership utm
        JOIN team_roles tr ON utm.role_id = tr.id
        WHERE tr.can_update_verifications = true
    ));
```

### API Resource (File: api/lab_verification_resources.py)
```python
from flask import request
from flask_restful import Resource, reqparse
from api.utils import handle_error, token_required, marshal_with
from api.models import LabVerification, Experiment
from api.schemas import lab_verification_fields, verification_stats_fields

class LabVerificationResource(Resource):
    """Resource for lab verification operations."""
    
    @token_required
    @marshal_with(lab_verification_fields)
    def get(self, experiment_id):
        """Get verification for an experiment."""
        try:
            return LabVerification.get_verification(experiment_id)
        except Exception as e:
            return handle_error(e)
    
    @token_required
    @marshal_with(lab_verification_fields)
    def post(self, experiment_id):
        """Record verification for an experiment."""
        try:
            parser = reqparse.RequestParser()
            parser.add_argument('verification_status', type=str, required=True,
                              help='Verification status is required')
            parser.add_argument('verifier', type=str, required=True,
                              help='Verifier is required')
            parser.add_argument('equipment_used', type=str, required=True,
                              help='Equipment used is required')
            parser.add_argument('comments', type=str)
            args = parser.parse_args()
            
            return LabVerification.record_verification(
                experiment_id=experiment_id,
                verification_status=args['verification_status'],
                verifier=args['verifier'],
                equipment_used=args['equipment_used'],
                comments=args['comments']
            )
        except Exception as e:
            return handle_error(e)
    
    @token_required
    @marshal_with(lab_verification_fields)
    def put(self, verification_id):
        """Update verification status."""
        try:
            parser = reqparse.RequestParser()
            parser.add_argument('verification_status', type=str, required=True,
                              help='Verification status is required')
            parser.add_argument('comments', type=str)
            args = parser.parse_args()
            
            return LabVerification.update_verification_status(
                verification_id=verification_id,
                new_status=args['verification_status'],
                comments=args['comments']
            )
        except Exception as e:
            return handle_error(e)

class VerificationStatsResource(Resource):
    """Resource for verification statistics."""
    
    @token_required
    @marshal_with(verification_stats_fields)
    def get(self):
        """Get verification statistics."""
        try:
            # Implementation details for getting statistics
            pass
        except Exception as e:
            return handle_error(e)
```

## Risk Management

### Technical Risks
- **Database schema changes**: Ensure proper migration testing
- **Integration with existing models**: Verify experiment model links correctly
- **API consistency**: Ensure new endpoints follow standardization

### Timeline Risks
- **Scope creep**: Focus strictly on essential functionality
- **Parallel work coordination**: Clear task boundaries and dependencies
- **Resource constraints**: Prioritize lab verification above other enhancements

## Next Immediate Steps

1. Implement the `LabVerification` class in `api/models.py`
2. Create the database migration script for lab verification schema
3. Implement the API resources for lab verification
4. Add lab verification components to the frontend
5. Create comprehensive tests for the verification workflow