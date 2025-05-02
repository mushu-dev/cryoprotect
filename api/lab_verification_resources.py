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
            records = LabVerification.get_all()
            total_count = len(records)
            verified_count = sum(1 for r in records if r.get("verification_status") == LabVerification.VERIFIED)
            pending_count = sum(1 for r in records if r.get("verification_status") == LabVerification.PENDING)
            rejected_count = sum(1 for r in records if r.get("verification_status") == LabVerification.REJECTED)
            verification_rate = (verified_count / total_count) if total_count > 0 else 0.0

            # Breakdown by equipment
            by_equipment = {}
            for r in records:
                equipment = r.get("equipment_used", "Unknown")
                status = r.get("verification_status")
                if equipment not in by_equipment:
                    by_equipment[equipment] = {"total": 0, "verified": 0, "pending": 0, "rejected": 0}
                by_equipment[equipment]["total"] += 1
                if status == LabVerification.VERIFIED:
                    by_equipment[equipment]["verified"] += 1
                elif status == LabVerification.PENDING:
                    by_equipment[equipment]["pending"] += 1
                elif status == LabVerification.REJECTED:
                    by_equipment[equipment]["rejected"] += 1

            # Breakdown by verifier
            by_verifier = {}
            for r in records:
                verifier = r.get("verifier", "Unknown")
                status = r.get("verification_status")
                if verifier not in by_verifier:
                    by_verifier[verifier] = {"total": 0, "verified": 0, "pending": 0, "rejected": 0}
                by_verifier[verifier]["total"] += 1
                if status == LabVerification.VERIFIED:
                    by_verifier[verifier]["verified"] += 1
                elif status == LabVerification.PENDING:
                    by_verifier[verifier]["pending"] += 1
                elif status == LabVerification.REJECTED:
                    by_verifier[verifier]["rejected"] += 1

            return {
                "total_count": total_count,
                "verified_count": verified_count,
                "pending_count": pending_count,
                "rejected_count": rejected_count,
                "verification_rate": verification_rate,
                "by_equipment": by_equipment,
                "by_verifier": by_verifier
            }
        except Exception as e:
            return handle_error(e)