/**
 * Lab Verification Module
 * 
 * This module handles the user interface for lab verification workflow,
 * including recording verification results, updating verification status,
 * and visualizing verification statistics.
 */

const VerificationModule = (function() {
    // Cache DOM elements
    let $verificationForm;
    let $statusSelect;
    let $equipmentInput;
    let $commentsInput;
    let $verificationList;
    let $statsContainer;
    
    // API endpoints
    const API_ENDPOINTS = {
        getVerification: (experimentId) => `/api/experiments/${experimentId}/verification`,
        updateVerification: (verificationId) => `/api/verifications/${verificationId}`,
        getStats: () => '/api/verification/stats'
    };
    
    // Verification status constants
    const STATUS = {
        PENDING: 'pending',
        VERIFIED: 'verified',
        REJECTED: 'rejected'
    };
    
    // Initialize the module
    function init() {
        // Initialize DOM elements
        cacheDomElements();
        bindEvents();
        
        // Load verification stats on initial load
        loadVerificationStats();
    }
    
    // Cache DOM elements for later use
    function cacheDomElements() {
        $verificationForm = $('#verification-form');
        $statusSelect = $('#verification-status');
        $equipmentInput = $('#verification-equipment');
        $commentsInput = $('#verification-comments');
        $verificationList = $('#verification-list');
        $statsContainer = $('#verification-stats');
    }
    
    // Bind event handlers
    function bindEvents() {
        // Form submission for recording verification
        $verificationForm.on('submit', handleVerificationSubmit);
        
        // Status update in verification list
        $verificationList.on('change', '.status-select', handleStatusUpdate);
        
        // Refresh stats button
        $('#refresh-stats').on('click', loadVerificationStats);
    }
    
    // Handle verification form submission
    function handleVerificationSubmit(e) {
        e.preventDefault();
        
        const experimentId = $verificationForm.data('experiment-id');
        const data = {
            verification_status: $statusSelect.val(),
            verifier: getCurrentUserId(),
            equipment_used: $equipmentInput.val(),
            comments: $commentsInput.val()
        };
        
        // Submit verification via API
        $.ajax({
            url: API_ENDPOINTS.getVerification(experimentId),
            method: 'POST',
            data: JSON.stringify(data),
            contentType: 'application/json',
            success: function(response) {
                showNotification('Verification recorded successfully', 'success');
                $verificationForm[0].reset();
                
                // Refresh verification list and stats
                loadExperimentVerification(experimentId);
                loadVerificationStats();
            },
            error: function(err) {
                showNotification('Error recording verification: ' + (err.responseJSON?.message || 'Unknown error'), 'error');
            }
        });
    }
    
    // Handle status update from verification list
    function handleStatusUpdate(e) {
        const $select = $(e.target);
        const verificationId = $select.data('verification-id');
        const newStatus = $select.val();
        const comments = prompt('Add comments for this status update (optional):');
        
        // Update verification status via API
        $.ajax({
            url: API_ENDPOINTS.updateVerification(verificationId),
            method: 'PUT',
            data: JSON.stringify({
                verification_status: newStatus,
                comments: comments
            }),
            contentType: 'application/json',
            success: function(response) {
                showNotification('Verification status updated successfully', 'success');
                
                // Refresh stats
                loadVerificationStats();
            },
            error: function(err) {
                showNotification('Error updating verification status: ' + (err.responseJSON?.message || 'Unknown error'), 'error');
                // Reset select to previous value
                $select.val($select.data('original-value'));
            }
        });
    }
    
    // Load verification for a specific experiment
    function loadExperimentVerification(experimentId) {
        $.ajax({
            url: API_ENDPOINTS.getVerification(experimentId),
            method: 'GET',
            success: function(verification) {
                if (verification) {
                    renderVerification(verification, experimentId);
                } else {
                    $verificationList.html('<p>No verification recorded for this experiment.</p>');
                }
            },
            error: function(err) {
                showNotification('Error loading verification: ' + (err.responseJSON?.message || 'Unknown error'), 'error');
            }
        });
    }
    
    // Load verification statistics
    function loadVerificationStats() {
        $.ajax({
            url: API_ENDPOINTS.getStats(),
            method: 'GET',
            success: function(stats) {
                renderStats(stats);
            },
            error: function(err) {
                showNotification('Error loading verification statistics: ' + (err.responseJSON?.message || 'Unknown error'), 'error');
            }
        });
    }
    
    // Render verification details
    function renderVerification(verification, experimentId) {
        const statusClass = {
            [STATUS.PENDING]: 'bg-warning',
            [STATUS.VERIFIED]: 'bg-success',
            [STATUS.REJECTED]: 'bg-danger'
        };
        
        const html = `
            <div class="card mt-3">
                <div class="card-header ${statusClass[verification.verification_status] || ''} text-white">
                    Verification Status: ${verification.verification_status.toUpperCase()}
                </div>
                <div class="card-body">
                    <h5 class="card-title">Verified by: ${verification.verifier}</h5>
                    <p><strong>Equipment Used:</strong> ${verification.equipment_used}</p>
                    <p><strong>Date:</strong> ${new Date(verification.created_at).toLocaleString()}</p>
                    ${verification.comments ? `<p><strong>Comments:</strong> ${verification.comments}</p>` : ''}
                    
                    <div class="form-group mt-3">
                        <label for="status-update-${verification.id}">Update Status:</label>
                        <select id="status-update-${verification.id}" 
                                class="form-control status-select" 
                                data-verification-id="${verification.id}"
                                data-original-value="${verification.verification_status}">
                            <option value="${STATUS.PENDING}" ${verification.verification_status === STATUS.PENDING ? 'selected' : ''}>Pending</option>
                            <option value="${STATUS.VERIFIED}" ${verification.verification_status === STATUS.VERIFIED ? 'selected' : ''}>Verified</option>
                            <option value="${STATUS.REJECTED}" ${verification.verification_status === STATUS.REJECTED ? 'selected' : ''}>Rejected</option>
                        </select>
                    </div>
                </div>
            </div>
        `;
        
        $verificationList.html(html);
    }
    
    // Render verification statistics
    function renderStats(stats) {
        if (!stats) {
            $statsContainer.html('<p>No verification statistics available.</p>');
            return;
        }
        
        const html = `
            <div class="card mt-3">
                <div class="card-header bg-primary text-white">
                    Verification Statistics
                </div>
                <div class="card-body">
                    <div class="row">
                        <div class="col-md-3">
                            <div class="stats-box">
                                <h5>Total Verifications</h5>
                                <div class="stats-value">${stats.total_count}</div>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stats-box bg-success text-white">
                                <h5>Verified</h5>
                                <div class="stats-value">${stats.verified_count}</div>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stats-box bg-warning text-white">
                                <h5>Pending</h5>
                                <div class="stats-value">${stats.pending_count}</div>
                            </div>
                        </div>
                        <div class="col-md-3">
                            <div class="stats-box bg-danger text-white">
                                <h5>Rejected</h5>
                                <div class="stats-value">${stats.rejected_count}</div>
                            </div>
                        </div>
                    </div>
                    
                    <div class="mt-4">
                        <h5>Verification Rate</h5>
                        <div class="progress">
                            <div class="progress-bar bg-success" role="progressbar" 
                                style="width: ${stats.verification_rate * 100}%" 
                                aria-valuenow="${stats.verification_rate * 100}" 
                                aria-valuemin="0" 
                                aria-valuemax="100">
                                ${Math.round(stats.verification_rate * 100)}%
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        `;
        
        $statsContainer.html(html);
    }
    
    // Helper function to get current user ID
    function getCurrentUserId() {
        return window.currentUser?.id || 'unknown_user';
    }
    
    // Helper function to show notifications
    function showNotification(message, type) {
        // Use application's notification system
        if (window.app && window.app.showNotification) {
            window.app.showNotification(message, type);
        } else {
            alert(message);
        }
    }
    
    // Public API
    return {
        init: init,
        loadExperimentVerification: loadExperimentVerification,
        loadVerificationStats: loadVerificationStats,
        // Export constants
        STATUS: STATUS
    };
})();

// Initialize on document ready
$(document).ready(function() {
    VerificationModule.init();
});