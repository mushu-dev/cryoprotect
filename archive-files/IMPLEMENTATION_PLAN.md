# CryoProtect v2 Phase 1: Lab Verification Workflow Implementation Plan

## 1. Overview

This plan details the step-by-step implementation of the essential lab verification workflow in CryoProtect v2. This feature is the highest priority for the project as it enables validation of prediction results against actual lab measurements.

## 2. Current State Summary

- **Experiments Module:** 
  - `api/models.py` contains the `Experiment` class for recording experimental results
  - Experiments are linked to mixtures but lack verification capability
  - No way to track lab verification status or record validation

- **Comparison Module:**
  - `api/comparisons.py` provides comparison between different mixtures
  - No capability to compare predicted vs actual lab-verified values
  - Missing statistical analysis for verification reporting

- **Frontend:**
  - Basic experiment viewing is available
  - No UI components for verification workflow
  - Missing visualization for comparing predicted vs actual results

- **Testing:**
  - Tests exist for experiments and comparisons
  - No tests for verification process

## 3. Implementation Tasks

### 3.1 Data Model Implementation

- [ ] **Create LabVerification model class in `api/models.py`**
  - Implement verification states (pending, verified, rejected)
  - Add methods for recording, retrieving, and updating verifications
  - Create relationships with Experiment model
  - Add authorization checks and error handling

- [ ] **Enhance Experiment model in `api/models.py`**
  - Add verification-related fields and methods
  - Update retrieval methods to include verification status
  - Maintain backward compatibility

- [ ] **Enhance Comparison model in `api/comparisons.py`**
  - Add methods for comparing predicted vs actual values
  - Implement statistical analysis for verification reporting
  - Add visualization data preparation

### 3.2 Database Schema Implementation

- [ ] **Create migration script `migrations/013_lab_verification_schema.sql`**
  - Define lab_verifications table with necessary fields
  - Set up foreign key relationships to experiments
  - Implement RLS policies for proper authorization
  - Add indexes for performance

- [ ] **Update existing tables if needed**
  - Add any needed fields to experiments table
  - Ensure backward compatibility

### 3.3 API Implementation

- [ ] **Create `api/lab_verification_resources.py`**
  - Implement RESTful endpoints for verification operations
  - Add proper authentication and authorization
  - Implement consistent error handling
  - Follow API standardization patterns

- [ ] **Update `app.py` to register new endpoints**
  - Add routes for verification resources
  - Set up appropriate URL patterns

- [ ] **Update API schemas in `api/schemas.py`**
  - Define field schemas for verification resources
  - Set up proper marshalling

### 3.4 Frontend Implementation

- [ ] **Create `static/js/verification.js`**
  - Implement verification UI components
  - Add forms for recording verification
  - Create visualization for comparing predictions vs actuals
  - Add status tracking and filtering

- [ ] **Update templates to include verification UI**
  - Add verification section to experiment details
  - Create verification dashboard view
  - Implement notification for pending verifications

### 3.5 Testing Implementation

- [ ] **Create `tests/test_lab_verification.py`**
  - Implement unit tests for LabVerification model
  - Add API endpoint tests for verification resources
  - Create integration tests for the complete workflow
  - Test authorization and validation logic

- [ ] **Update `tests/test_models.py` and `tests/test_comparisons.py`**
  - Add tests for enhanced experiment functionality
  - Test comparison methods with verification data

## 4. Task Dependencies & Timeline

| Task | Dependency | Timeline |
|------|------------|----------|
| LabVerification model | None | Days 1-2 |
| Database migration | LabVerification model | Day 3 |
| API resources | LabVerification model, Database migration | Days 4-5 |
| Frontend components | API resources | Days 6-8 |
| Testing | All of the above | Days 9-10 |

## 5. Implementation Details

### 5.1 LabVerification Model

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
        supabase = get_supabase_client()
        
        experiment = Experiment.get_by_id(experiment_id)
        if not experiment:
            raise ValueError(f"Experiment with ID {experiment_id} not found")
            
        if verification_status not in [cls.PENDING, cls.VERIFIED, cls.REJECTED]:
            raise ValueError(f"Invalid verification status: {verification_status}")
        
        verification_data = {
            "experiment_id": experiment_id,
            "verification_status": verification_status,
            "verifier": verifier,
            "equipment_used": equipment_used,
            "comments": comments,
        }
        
        try:
            result = supabase.table(cls.table_name).insert(verification_data).execute()
            return result.data[0] if result.data else None
        except Exception as e:
            logger.error(f"Error recording verification: {str(e)}")
            raise
    
    @classmethod
    def get_verification(cls, experiment_id: str) -> Dict[str, Any]:
        """Get verification for an experiment."""
        supabase = get_supabase_client()
        
        try:
            result = supabase.table(cls.table_name) \
                .select("*") \
                .eq("experiment_id", experiment_id) \
                .execute()
            return result.data[0] if result.data else None
        except Exception as e:
            logger.error(f"Error retrieving verification: {str(e)}")
            raise
    
    @classmethod
    def update_verification_status(cls, verification_id: str, 
                                  new_status: str,
                                  comments: str = None) -> Dict[str, Any]:
        """Update verification status."""
        supabase = get_supabase_client()
        
        if new_status not in [cls.PENDING, cls.VERIFIED, cls.REJECTED]:
            raise ValueError(f"Invalid verification status: {new_status}")
        
        update_data = {
            "verification_status": new_status,
            "updated_at": datetime.now().isoformat(),
        }
        
        if comments:
            update_data["comments"] = comments
        
        try:
            result = supabase.table(cls.table_name) \
                .update(update_data) \
                .eq("id", verification_id) \
                .execute()
            return result.data[0] if result.data else None
        except Exception as e:
            logger.error(f"Error updating verification status: {str(e)}")
            raise
```

### 5.2 Database Migration

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

-- Add index for performance
CREATE INDEX idx_lab_verifications_experiment_id ON public.lab_verifications(experiment_id);
CREATE INDEX idx_lab_verifications_status ON public.lab_verifications(verification_status);
```

### 5.3 API Resources

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
            # Implementation for retrieving verification statistics
            # - Count by status
            # - Verification accuracy metrics
            # - Prediction vs actual comparison stats
            pass
        except Exception as e:
            return handle_error(e)
```

## 6. Success Criteria

- Lab verification workflow fully implemented and tested
- Users can record verification status for experiments
- System tracks verification status (pending, verified, rejected)
- Users can compare predicted vs actual values
- Statistical analysis provides insight into prediction accuracy
- All functionality properly secured with appropriate authorization

## 7. Test Plan

- **Unit Tests:** Test each model method and API endpoint individually
- **Integration Tests:** Test the complete verification workflow
- **Authorization Tests:** Verify proper access controls
- **Edge Cases:** Test validation, error handling, and boundary conditions
- **Performance:** Ensure efficiency with large datasets

## 8. Risk Mitigation

- **Database Migration:** Backup before applying, test thoroughly
- **API Consistency:** Follow established patterns exactly
- **Authorization:** Double-check RLS policies and access controls
- **Integration:** Ensure compatibility with existing components

---

*This plan will be updated as implementation progresses.*