/**
 * Lab Verification System Component
 * 
 * This component provides a user interface for requesting, reviewing, and managing
 * lab verifications for experiments. It includes forms for submitting verification requests,
 * displays for showing verification status, and interfaces for verifiers to review and
 * rate experiments.
 */

import React, { useState, useEffect } from 'react';
import { useQuery, useMutation } from 'convex/react';
import { api } from '../../../convex/_generated/api';
import { Id } from '../../../convex/_generated/dataModel';
import { useToast } from '../../ui/use-toast';
import { useErrorHandler } from '../../errors/use-error-handler';

// Types
export type VerificationStatus = 'pending' | 'verified' | 'rejected' | 'needs_revision';

interface LabVerificationProps {
  experimentId: Id<"enhancedExperiments">;
  isVerifier?: boolean;
  userId?: Id<"users">;
}

/**
 * Main Lab Verification component
 */
export function LabVerificationSystem({
  experimentId,
  isVerifier = false,
  userId
}: LabVerificationProps) {
  // Get the verification for this experiment
  const verification = useQuery(
    api.labVerifications.getVerificationByExperiment,
    { experimentId, includeExperiment: true }
  );
  
  // Get experiment details
  const experiment = useQuery(
    api.experiments.enhanced_experiments.getEnhancedExperiment,
    { experimentId }
  );
  
  // States
  const [showRequestForm, setShowRequestForm] = useState(false);
  const [showVerifyForm, setShowVerifyForm] = useState(false);
  
  // Toast notifications
  const { toast } = useToast();
  const { handleError } = useErrorHandler();
  
  // If data is loading
  if (verification === undefined || experiment === undefined) {
    return <div className="p-4 bg-gray-100 rounded-md">Loading verification data...</div>;
  }
  
  // Determine if the current user can request verification
  const canRequestVerification = 
    !verification && 
    experiment && 
    userId && 
    (experiment.conductedBy?.equals(userId) || isVerifier);
  
  // Determine if the current user can verify this experiment
  const canVerify = 
    verification && 
    verification.verificationStatus === 'pending' && 
    isVerifier && 
    userId;
  
  // Determine if the current user can revise this experiment's verification
  const canRevise = 
    verification && 
    verification.verificationStatus === 'needs_revision' && 
    experiment && 
    experiment.conductedBy?.equals(userId);
  
  return (
    <div className="bg-white rounded-lg shadow-md p-6 mb-8">
      <h2 className="text-2xl font-bold mb-4">Lab Verification</h2>
      
      {/* Display verification status */}
      {verification ? (
        <div className="mb-6">
          <VerificationStatusBadge status={verification.verificationStatus} />
          
          <div className="mt-4 grid grid-cols-1 md:grid-cols-2 gap-4">
            <div>
              <h3 className="text-lg font-semibold mb-2">Verification Details</h3>
              <div className="text-sm space-y-1">
                <p><span className="font-medium">Requested:</span> {formatDate(verification.requestDate)}</p>
                {verification.verificationDate && (
                  <p><span className="font-medium">Verified:</span> {formatDate(verification.verificationDate)}</p>
                )}
                {verification.equipmentUsed && (
                  <p><span className="font-medium">Equipment Used:</span> {Array.isArray(verification.equipmentUsed) ? verification.equipmentUsed.join(', ') : verification.equipmentUsed}</p>
                )}
              </div>
              
              {verification.requestNotes && (
                <div className="mt-3">
                  <h4 className="font-medium">Request Notes</h4>
                  <p className="text-sm mt-1 bg-gray-50 p-2 rounded">{verification.requestNotes}</p>
                </div>
              )}
              
              {verification.verifierNotes && (
                <div className="mt-3">
                  <h4 className="font-medium">Verifier Notes</h4>
                  <p className="text-sm mt-1 bg-gray-50 p-2 rounded">{verification.verifierNotes}</p>
                </div>
              )}
            </div>
            
            {/* Display ratings if verified */}
            {verification.verificationStatus === 'verified' && (
              <div>
                <h3 className="text-lg font-semibold mb-2">Quality Ratings</h3>
                <div className="grid grid-cols-1 gap-2">
                  {verification.reproducibilityRating && (
                    <RatingDisplay 
                      label="Reproducibility" 
                      value={verification.reproducibilityRating} 
                    />
                  )}
                  {verification.qualityRating && (
                    <RatingDisplay 
                      label="Quality" 
                      value={verification.qualityRating} 
                    />
                  )}
                  {verification.documentationRating && (
                    <RatingDisplay 
                      label="Documentation" 
                      value={verification.documentationRating} 
                    />
                  )}
                  {verification.overallRating && (
                    <RatingDisplay 
                      label="Overall" 
                      value={verification.overallRating} 
                      isOverall
                    />
                  )}
                </div>
              </div>
            )}
          </div>
          
          {/* Verification evidence */}
          {verification.evidenceUrls && verification.evidenceUrls.length > 0 && (
            <div className="mt-4">
              <h3 className="text-lg font-semibold mb-2">Evidence</h3>
              <div className="flex flex-wrap gap-2">
                {verification.evidenceUrls.map((url, index) => (
                  <a
                    key={index}
                    href={url}
                    target="_blank"
                    rel="noopener noreferrer"
                    className="px-3 py-1 bg-blue-100 text-blue-700 rounded-full text-sm hover:bg-blue-200 transition-colors"
                  >
                    Evidence #{index + 1}
                  </a>
                ))}
              </div>
            </div>
          )}
          
          {/* Action buttons */}
          <div className="mt-6 flex gap-4">
            {canVerify && (
              <button
                onClick={() => setShowVerifyForm(true)}
                className="px-4 py-2 bg-green-600 text-white rounded-md hover:bg-green-700 transition-colors"
              >
                Review Verification
              </button>
            )}
            
            {canRevise && (
              <button
                onClick={() => setShowRequestForm(true)}
                className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
              >
                Submit Revision
              </button>
            )}
          </div>
          
          {/* Verification dialog */}
          {showVerifyForm && verification && (
            <VerificationReviewForm
              verificationId={verification._id}
              onClose={() => setShowVerifyForm(false)}
              onSuccess={() => {
                setShowVerifyForm(false);
                toast({ 
                  title: "Verification Updated", 
                  description: "The verification has been successfully updated." 
                });
              }}
              onError={(error) => {
                handleError(error, "Failed to update verification");
                toast({
                  title: "Update Failed",
                  description: "There was an error updating the verification. Please try again.",
                  variant: "destructive"
                });
              }}
            />
          )}
        </div>
      ) : (
        <div className="mb-6">
          <div className="p-4 bg-gray-100 rounded-md">
            <p className="text-gray-700">
              This experiment has not been verified yet. Verification ensures the reproducibility and quality of experimental procedures and results.
            </p>
            
            {canRequestVerification && (
              <button
                onClick={() => setShowRequestForm(true)}
                className="mt-3 px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
              >
                Request Verification
              </button>
            )}
          </div>
        </div>
      )}
      
      {/* Request form dialog */}
      {showRequestForm && (
        <VerificationRequestForm
          experimentId={experimentId}
          isRevision={!!verification && verification.verificationStatus === 'needs_revision'}
          onClose={() => setShowRequestForm(false)}
          onSuccess={() => {
            setShowRequestForm(false);
            toast({
              title: "Verification Requested", 
              description: "Your verification request has been submitted successfully."
            });
          }}
          onError={(error) => {
            handleError(error, "Failed to request verification");
            toast({
              title: "Request Failed",
              description: "There was an error submitting your verification request. Please try again.",
              variant: "destructive"
            });
          }}
        />
      )}
    </div>
  );
}

/**
 * Component to display verification status
 */
function VerificationStatusBadge({ status }: { status: VerificationStatus }) {
  let bgColor;
  let textColor;
  let label;
  
  switch (status) {
    case 'pending':
      bgColor = 'bg-yellow-100';
      textColor = 'text-yellow-800';
      label = 'Pending Verification';
      break;
    case 'verified':
      bgColor = 'bg-green-100';
      textColor = 'text-green-800';
      label = 'Verified';
      break;
    case 'rejected':
      bgColor = 'bg-red-100';
      textColor = 'text-red-800';
      label = 'Verification Rejected';
      break;
    case 'needs_revision':
      bgColor = 'bg-orange-100';
      textColor = 'text-orange-800';
      label = 'Needs Revision';
      break;
    default:
      bgColor = 'bg-gray-100';
      textColor = 'text-gray-800';
      label = 'Unknown Status';
  }
  
  return (
    <div className={`inline-flex items-center px-3 py-1 rounded-full text-sm font-medium ${bgColor} ${textColor}`}>
      {status === 'pending' && <span className="w-2 h-2 bg-yellow-500 rounded-full mr-1.5 animate-pulse"></span>}
      {status === 'verified' && <span className="w-2 h-2 bg-green-500 rounded-full mr-1.5"></span>}
      {status === 'rejected' && <span className="w-2 h-2 bg-red-500 rounded-full mr-1.5"></span>}
      {status === 'needs_revision' && <span className="w-2 h-2 bg-orange-500 rounded-full mr-1.5"></span>}
      {label}
    </div>
  );
}

/**
 * Component for requesting verification
 */
function VerificationRequestForm({
  experimentId,
  isRevision = false,
  onClose,
  onSuccess,
  onError
}: {
  experimentId: Id<"enhancedExperiments">;
  isRevision?: boolean;
  onClose: () => void;
  onSuccess: () => void;
  onError: (error: Error) => void;
}) {
  // State for form fields
  const [equipmentUsed, setEquipmentUsed] = useState<string[]>([]);
  const [equipmentInput, setEquipmentInput] = useState('');
  const [requestNotes, setRequestNotes] = useState('');
  const [methodologyDescription, setMethodologyDescription] = useState('');
  const [controlProcedures, setControlProcedures] = useState('');
  const [evidenceUrls, setEvidenceUrls] = useState<string[]>([]);
  const [evidenceInput, setEvidenceInput] = useState('');
  
  // Mutation for requesting verification
  const requestVerification = useMutation(
    api.labVerifications.requestVerification
  );
  
  // Mutation for updating verification
  const updateVerification = useMutation(
    api.labVerifications.updateVerification
  );
  
  // Get the existing verification if this is a revision
  const verification = useQuery(
    api.labVerifications.getVerificationByExperiment,
    isRevision ? { experimentId } : "skip"
  );
  
  // Prefill form data if this is a revision
  useEffect(() => {
    if (isRevision && verification) {
      if (verification.equipmentUsed) {
        setEquipmentUsed(Array.isArray(verification.equipmentUsed) 
          ? verification.equipmentUsed 
          : [verification.equipmentUsed]);
      }
      
      if (verification.requestNotes) {
        setRequestNotes(verification.requestNotes);
      }
      
      if (verification.methodologyDescription) {
        setMethodologyDescription(verification.methodologyDescription);
      }
      
      if (verification.controlProcedures) {
        setControlProcedures(verification.controlProcedures);
      }
      
      if (verification.evidenceUrls) {
        setEvidenceUrls(verification.evidenceUrls);
      }
    }
  }, [isRevision, verification]);
  
  // Add equipment
  const addEquipment = () => {
    if (equipmentInput.trim()) {
      setEquipmentUsed([...equipmentUsed, equipmentInput.trim()]);
      setEquipmentInput('');
    }
  };
  
  // Remove equipment
  const removeEquipment = (index: number) => {
    setEquipmentUsed(equipmentUsed.filter((_, i) => i !== index));
  };
  
  // Add evidence URL
  const addEvidenceUrl = () => {
    if (evidenceInput.trim()) {
      try {
        // Basic URL validation
        new URL(evidenceInput.trim());
        setEvidenceUrls([...evidenceUrls, evidenceInput.trim()]);
        setEvidenceInput('');
      } catch (error) {
        alert('Please enter a valid URL');
      }
    }
  };
  
  // Remove evidence URL
  const removeEvidenceUrl = (index: number) => {
    setEvidenceUrls(evidenceUrls.filter((_, i) => i !== index));
  };
  
  // Submit the verification request
  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    
    if (equipmentUsed.length === 0) {
      alert('Please add at least one piece of equipment');
      return;
    }
    
    try {
      if (isRevision && verification) {
        // Update existing verification
        await updateVerification({
          verificationId: verification._id,
          update: {
            equipmentUsed,
            requestNotes,
            methodologyDescription,
            controlProcedures,
            evidenceUrls: evidenceUrls.length > 0 ? evidenceUrls : undefined
          }
        });
      } else {
        // Create new verification request
        await requestVerification({
          verification: {
            experimentId,
            equipmentUsed,
            requestNotes: requestNotes || undefined,
            methodologyDescription: methodologyDescription || undefined,
            controlProcedures: controlProcedures || undefined,
            evidenceUrls: evidenceUrls.length > 0 ? evidenceUrls : undefined
          }
        });
      }
      
      onSuccess();
    } catch (error) {
      onError(error instanceof Error ? error : new Error('Failed to request verification'));
    }
  };
  
  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center p-4 z-50">
      <div className="bg-white rounded-lg shadow-xl p-6 max-w-2xl w-full max-h-[90vh] overflow-y-auto">
        <h2 className="text-xl font-bold mb-4">
          {isRevision ? 'Submit Verification Revision' : 'Request Verification'}
        </h2>
        
        <form onSubmit={handleSubmit}>
          {/* Equipment Used */}
          <div className="mb-4">
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Equipment Used <span className="text-red-500">*</span>
            </label>
            <div className="flex mb-2">
              <input
                type="text"
                value={equipmentInput}
                onChange={(e) => setEquipmentInput(e.target.value)}
                placeholder="Enter equipment name"
                className="flex-1 px-3 py-2 border border-gray-300 rounded-l-md focus:outline-none focus:ring-1 focus:ring-blue-500"
              />
              <button
                type="button"
                onClick={addEquipment}
                className="px-4 py-2 bg-blue-600 text-white rounded-r-md hover:bg-blue-700 transition-colors"
              >
                Add
              </button>
            </div>
            {equipmentUsed.length > 0 && (
              <div className="flex flex-wrap gap-2 mt-2">
                {equipmentUsed.map((equipment, index) => (
                  <div 
                    key={index}
                    className="px-3 py-1 bg-blue-100 text-blue-800 rounded-full text-sm flex items-center"
                  >
                    {equipment}
                    <button
                      type="button"
                      onClick={() => removeEquipment(index)}
                      className="ml-1.5 text-blue-500 hover:text-blue-700"
                    >
                      &times;
                    </button>
                  </div>
                ))}
              </div>
            )}
            {equipmentUsed.length === 0 && (
              <p className="text-sm text-red-500">At least one piece of equipment is required</p>
            )}
          </div>
          
          {/* Methodology Description */}
          <div className="mb-4">
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Methodology Description
            </label>
            <textarea
              value={methodologyDescription}
              onChange={(e) => setMethodologyDescription(e.target.value)}
              rows={3}
              placeholder="Describe the methodology used in this experiment"
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-1 focus:ring-blue-500"
            />
          </div>
          
          {/* Control Procedures */}
          <div className="mb-4">
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Control Procedures
            </label>
            <textarea
              value={controlProcedures}
              onChange={(e) => setControlProcedures(e.target.value)}
              rows={3}
              placeholder="Describe any control procedures used in this experiment"
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-1 focus:ring-blue-500"
            />
          </div>
          
          {/* Request Notes */}
          <div className="mb-4">
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Request Notes
            </label>
            <textarea
              value={requestNotes}
              onChange={(e) => setRequestNotes(e.target.value)}
              rows={3}
              placeholder="Additional notes about this verification request"
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-1 focus:ring-blue-500"
            />
          </div>
          
          {/* Evidence URLs */}
          <div className="mb-4">
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Evidence URLs
            </label>
            <div className="flex mb-2">
              <input
                type="text"
                value={evidenceInput}
                onChange={(e) => setEvidenceInput(e.target.value)}
                placeholder="https://example.com/evidence"
                className="flex-1 px-3 py-2 border border-gray-300 rounded-l-md focus:outline-none focus:ring-1 focus:ring-blue-500"
              />
              <button
                type="button"
                onClick={addEvidenceUrl}
                className="px-4 py-2 bg-blue-600 text-white rounded-r-md hover:bg-blue-700 transition-colors"
              >
                Add
              </button>
            </div>
            {evidenceUrls.length > 0 && (
              <div className="flex flex-wrap gap-2 mt-2">
                {evidenceUrls.map((url, index) => (
                  <div 
                    key={index}
                    className="px-3 py-1 bg-blue-100 text-blue-800 rounded-full text-sm flex items-center"
                  >
                    <a 
                      href={url} 
                      target="_blank"
                      rel="noopener noreferrer"
                      className="hover:underline"
                    >
                      Evidence #{index + 1}
                    </a>
                    <button
                      type="button"
                      onClick={() => removeEvidenceUrl(index)}
                      className="ml-1.5 text-blue-500 hover:text-blue-700"
                    >
                      &times;
                    </button>
                  </div>
                ))}
              </div>
            )}
          </div>
          
          {/* Form Actions */}
          <div className="flex justify-end gap-2 mt-6">
            <button
              type="button"
              onClick={onClose}
              className="px-4 py-2 border border-gray-300 text-gray-700 rounded-md hover:bg-gray-100 transition-colors"
            >
              Cancel
            </button>
            <button
              type="submit"
              className="px-4 py-2 bg-blue-600 text-white rounded-md hover:bg-blue-700 transition-colors"
              disabled={equipmentUsed.length === 0}
            >
              {isRevision ? 'Submit Revision' : 'Request Verification'}
            </button>
          </div>
        </form>
      </div>
    </div>
  );
}

/**
 * Component for verifying an experiment
 */
function VerificationReviewForm({
  verificationId,
  onClose,
  onSuccess,
  onError
}: {
  verificationId: Id<"labVerifications">;
  onClose: () => void;
  onSuccess: () => void;
  onError: (error: Error) => void;
}) {
  // State for form fields
  const [status, setStatus] = useState<VerificationStatus>('verified');
  const [verifierNotes, setVerifierNotes] = useState('');
  const [reproducibilityRating, setReproducibilityRating] = useState(5);
  const [qualityRating, setQualityRating] = useState(5);
  const [documentationRating, setDocumentationRating] = useState(5);
  const [evidenceUrls, setEvidenceUrls] = useState<string[]>([]);
  const [evidenceInput, setEvidenceInput] = useState('');
  
  // Get the verification
  const verification = useQuery(
    api.labVerifications.getVerification,
    { verificationId }
  );
  
  // Mutation for verifying an experiment
  const verifyExperiment = useMutation(
    api.labVerifications.verifyExperiment
  );
  
  // Prefill form data
  useEffect(() => {
    if (verification) {
      if (verification.verifierNotes) {
        setVerifierNotes(verification.verifierNotes);
      }
      
      if (verification.reproducibilityRating) {
        setReproducibilityRating(verification.reproducibilityRating);
      }
      
      if (verification.qualityRating) {
        setQualityRating(verification.qualityRating);
      }
      
      if (verification.documentationRating) {
        setDocumentationRating(verification.documentationRating);
      }
      
      if (verification.evidenceUrls) {
        setEvidenceUrls(verification.evidenceUrls);
      }
    }
  }, [verification]);
  
  // Add evidence URL
  const addEvidenceUrl = () => {
    if (evidenceInput.trim()) {
      try {
        // Basic URL validation
        new URL(evidenceInput.trim());
        setEvidenceUrls([...evidenceUrls, evidenceInput.trim()]);
        setEvidenceInput('');
      } catch (error) {
        alert('Please enter a valid URL');
      }
    }
  };
  
  // Remove evidence URL
  const removeEvidenceUrl = (index: number) => {
    setEvidenceUrls(evidenceUrls.filter((_, i) => i !== index));
  };
  
  // Calculate overall rating
  const overallRating = Math.round(
    (reproducibilityRating * 0.4) +
    (qualityRating * 0.4) +
    (documentationRating * 0.2)
  );
  
  // Submit the verification
  const handleSubmit = async (e: React.FormEvent) => {
    e.preventDefault();
    
    // Check if ratings are required for 'verified' status
    if (status === 'verified' && (
      reproducibilityRating < 1 || reproducibilityRating > 10 ||
      qualityRating < 1 || qualityRating > 10 ||
      documentationRating < 1 || documentationRating > 10
    )) {
      alert('All ratings must be between 1 and 10 for verified status');
      return;
    }
    
    try {
      await verifyExperiment({
        verificationId,
        status,
        verifierNotes: verifierNotes || undefined,
        reproducibilityRating: status === 'verified' ? reproducibilityRating : undefined,
        qualityRating: status === 'verified' ? qualityRating : undefined,
        documentationRating: status === 'verified' ? documentationRating : undefined,
        evidenceUrls: evidenceUrls.length > 0 ? evidenceUrls : undefined
      });
      
      onSuccess();
    } catch (error) {
      onError(error instanceof Error ? error : new Error('Failed to verify experiment'));
    }
  };
  
  // If verification is loading
  if (verification === undefined) {
    return (
      <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center p-4 z-50">
        <div className="bg-white rounded-lg shadow-xl p-6">
          <p>Loading verification data...</p>
        </div>
      </div>
    );
  }
  
  return (
    <div className="fixed inset-0 bg-black bg-opacity-50 flex items-center justify-center p-4 z-50">
      <div className="bg-white rounded-lg shadow-xl p-6 max-w-2xl w-full max-h-[90vh] overflow-y-auto">
        <h2 className="text-xl font-bold mb-4">Review Verification</h2>
        
        <form onSubmit={handleSubmit}>
          {/* Verification Status */}
          <div className="mb-4">
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Verification Status <span className="text-red-500">*</span>
            </label>
            <select
              value={status}
              onChange={(e) => setStatus(e.target.value as VerificationStatus)}
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-1 focus:ring-blue-500"
            >
              <option value="verified">Verified ✓</option>
              <option value="rejected">Rejected ✗</option>
              <option value="needs_revision">Needs Revision ⚠</option>
            </select>
          </div>
          
          {/* Ratings (only shown for 'verified' status) */}
          {status === 'verified' && (
            <div className="mb-4 space-y-3">
              <h3 className="font-medium text-gray-700">Experiment Ratings</h3>
              
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-1">
                  Reproducibility (40%) <span className="text-red-500">*</span>
                </label>
                <div className="flex items-center space-x-2">
                  <input
                    type="range"
                    min="1"
                    max="10"
                    value={reproducibilityRating}
                    onChange={(e) => setReproducibilityRating(parseInt(e.target.value))}
                    className="flex-1"
                  />
                  <span className="w-8 text-center font-medium">{reproducibilityRating}</span>
                </div>
              </div>
              
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-1">
                  Quality (40%) <span className="text-red-500">*</span>
                </label>
                <div className="flex items-center space-x-2">
                  <input
                    type="range"
                    min="1"
                    max="10"
                    value={qualityRating}
                    onChange={(e) => setQualityRating(parseInt(e.target.value))}
                    className="flex-1"
                  />
                  <span className="w-8 text-center font-medium">{qualityRating}</span>
                </div>
              </div>
              
              <div>
                <label className="block text-sm font-medium text-gray-700 mb-1">
                  Documentation (20%) <span className="text-red-500">*</span>
                </label>
                <div className="flex items-center space-x-2">
                  <input
                    type="range"
                    min="1"
                    max="10"
                    value={documentationRating}
                    onChange={(e) => setDocumentationRating(parseInt(e.target.value))}
                    className="flex-1"
                  />
                  <span className="w-8 text-center font-medium">{documentationRating}</span>
                </div>
              </div>
              
              <div className="mt-3 pt-3 border-t border-gray-200">
                <label className="block text-sm font-medium text-gray-700 mb-1">
                  Overall Rating
                </label>
                <div className="flex items-center space-x-2">
                  <div className="flex-1 h-4 bg-gray-200 rounded-full overflow-hidden">
                    <div 
                      className="h-full bg-blue-600" 
                      style={{ width: `${(overallRating / 10) * 100}%` }}
                    ></div>
                  </div>
                  <span className="w-8 text-center font-medium">{overallRating}</span>
                </div>
              </div>
            </div>
          )}
          
          {/* Verifier Notes */}
          <div className="mb-4">
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Verifier Notes {status !== 'verified' && <span className="text-red-500">*</span>}
            </label>
            <textarea
              value={verifierNotes}
              onChange={(e) => setVerifierNotes(e.target.value)}
              rows={4}
              placeholder={status === 'needs_revision' 
                ? "Explain what needs to be revised" 
                : (status === 'rejected' 
                  ? "Explain why the verification was rejected" 
                  : "Additional notes about the verification")}
              className="w-full px-3 py-2 border border-gray-300 rounded-md focus:outline-none focus:ring-1 focus:ring-blue-500"
              required={status !== 'verified'}
            />
          </div>
          
          {/* Evidence URLs */}
          <div className="mb-4">
            <label className="block text-sm font-medium text-gray-700 mb-1">
              Evidence URLs
            </label>
            <div className="flex mb-2">
              <input
                type="text"
                value={evidenceInput}
                onChange={(e) => setEvidenceInput(e.target.value)}
                placeholder="https://example.com/evidence"
                className="flex-1 px-3 py-2 border border-gray-300 rounded-l-md focus:outline-none focus:ring-1 focus:ring-blue-500"
              />
              <button
                type="button"
                onClick={addEvidenceUrl}
                className="px-4 py-2 bg-blue-600 text-white rounded-r-md hover:bg-blue-700 transition-colors"
              >
                Add
              </button>
            </div>
            {evidenceUrls.length > 0 && (
              <div className="flex flex-wrap gap-2 mt-2">
                {evidenceUrls.map((url, index) => (
                  <div 
                    key={index}
                    className="px-3 py-1 bg-blue-100 text-blue-800 rounded-full text-sm flex items-center"
                  >
                    <a 
                      href={url} 
                      target="_blank"
                      rel="noopener noreferrer"
                      className="hover:underline"
                    >
                      Evidence #{index + 1}
                    </a>
                    <button
                      type="button"
                      onClick={() => removeEvidenceUrl(index)}
                      className="ml-1.5 text-blue-500 hover:text-blue-700"
                    >
                      &times;
                    </button>
                  </div>
                ))}
              </div>
            )}
          </div>
          
          {/* Form Actions */}
          <div className="flex justify-end gap-2 mt-6">
            <button
              type="button"
              onClick={onClose}
              className="px-4 py-2 border border-gray-300 text-gray-700 rounded-md hover:bg-gray-100 transition-colors"
            >
              Cancel
            </button>
            <button
              type="submit"
              className={`px-4 py-2 text-white rounded-md transition-colors ${
                status === 'verified' 
                  ? 'bg-green-600 hover:bg-green-700' 
                  : (status === 'rejected' 
                    ? 'bg-red-600 hover:bg-red-700' 
                    : 'bg-orange-600 hover:bg-orange-700')
              }`}
            >
              {status === 'verified' 
                ? 'Verify Experiment' 
                : (status === 'rejected' 
                  ? 'Reject Verification' 
                  : 'Request Revisions')}
            </button>
          </div>
        </form>
      </div>
    </div>
  );
}

/**
 * Component to display a rating
 */
function RatingDisplay({
  label,
  value,
  isOverall = false
}: {
  label: string;
  value: number;
  isOverall?: boolean;
}) {
  // Get color based on rating
  const getColor = (rating: number) => {
    if (rating >= 8) return 'bg-green-500';
    if (rating >= 6) return 'bg-blue-500';
    if (rating >= 4) return 'bg-yellow-500';
    return 'bg-red-500';
  };
  
  return (
    <div>
      <div className="flex justify-between mb-1">
        <span className="text-sm font-medium">{label}</span>
        <span className={`text-sm font-medium ${isOverall ? 'font-bold' : ''}`}>{value}/10</span>
      </div>
      <div className="w-full bg-gray-200 rounded-full h-2.5">
        <div
          className={`h-2.5 rounded-full ${getColor(value)} ${isOverall ? 'h-3' : ''}`}
          style={{ width: `${(value / 10) * 100}%` }}
        ></div>
      </div>
    </div>
  );
}

/**
 * Helper to format dates
 */
function formatDate(timestamp: number): string {
  return new Date(timestamp).toLocaleString();
}