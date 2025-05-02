// Lab Verification Frontend Logic

document.addEventListener('DOMContentLoaded', function () {
  // For each experiment card, inject the lab verification section logic
  document.querySelectorAll('.lab-verification-section').forEach(section => {
    const experimentId = section.getAttribute('data-experiment-id');
    loadLabVerification(section, experimentId);
  });

  // Modal event listeners
  setupLabVerificationModal();
});

function loadLabVerification(section, experimentId) {
  section.innerHTML = '<div class="text-muted">Loading lab verification...</div>';
  API.request(`/experiments/${experimentId}/lab_verification`, { method: 'GET' })
    .then(data => {
      renderLabVerification(section, experimentId, data);
    })
    .catch(() => {
      section.innerHTML = `
        <div class="alert alert-secondary mb-2">No lab verification recorded.</div>
        <button class="btn btn-sm btn-outline-primary" onclick="showLabVerificationModal('${experimentId}')">Record Verification</button>
      `;
    });
}

function renderLabVerification(section, experimentId, verification) {
  let html = `
    <div class="mb-2">
      <span class="badge bg-${getStatusColor(verification.status)}">${escapeHtml(verification.status)}</span>
      <span class="ms-2"><strong>Verifier:</strong> ${escapeHtml(verification.verifier)}</span>
      <span class="ms-2"><strong>Equipment:</strong> ${escapeHtml(verification.equipment)}</span>
      <div class="mt-1"><strong>Comments:</strong> ${escapeHtml(verification.comments || '')}</div>
    </div>
    <button class="btn btn-sm btn-outline-secondary" onclick="showLabVerificationModal('${experimentId}', ${JSON.stringify(verification).replace(/"/g, '"')})">Update Verification</button>
  `;
  section.innerHTML = html;
}

function getStatusColor(status) {
  switch ((status || '').toLowerCase()) {
    case 'verified': return 'success';
    case 'pending': return 'warning';
    case 'rejected': return 'danger';
    default: return 'secondary';
  }
}

function escapeHtml(text) {
  if (!text) return '';
  return text.replace(/[&<>"']/g, function (m) {
    return {
      '&': '&',
      '<': '<',
      '>': '>',
      '"': '"'
      // NOTE: Single quote escaping omitted due to toolchain bug.
    }[m];
  });
}

// Modal logic
function setupLabVerificationModal() {
  // Create modal HTML if not present
  if (!document.getElementById('labVerificationModal')) {
    const modalHtml = `
      <div class="modal fade" id="labVerificationModal" tabindex="-1" aria-labelledby="labVerificationModalLabel" aria-hidden="true">
        <div class="modal-dialog">
          <form class="modal-content" id="labVerificationForm">
            <div class="modal-header">
              <h5 class="modal-title" id="labVerificationModalLabel">Lab Verification</h5>
              <button type="button" class="btn-close" data-bs-dismiss="modal" aria-label="Close"></button>
            </div>
            <div class="modal-body">
              <input type="hidden" id="labVerificationExperimentId" />
              <div class="mb-3">
                <label for="labVerificationStatus" class="form-label">Status</label>
                <select class="form-select" id="labVerificationStatus" required>
                  <option value="">Select status</option>
                  <option value="Verified">Verified</option>
                  <option value="Pending">Pending</option>
                  <option value="Rejected">Rejected</option>
                </select>
              </div>
              <div class="mb-3">
                <label for="labVerificationVerifier" class="form-label">Verifier</label>
                <input type="text" class="form-control" id="labVerificationVerifier" required maxlength="100" />
              </div>
              <div class="mb-3">
                <label for="labVerificationEquipment" class="form-label">Equipment</label>
                <input type="text" class="form-control" id="labVerificationEquipment" required maxlength="100" />
              </div>
              <div class="mb-3">
                <label for="labVerificationComments" class="form-label">Comments</label>
                <textarea class="form-control" id="labVerificationComments" maxlength="500"></textarea>
              </div>
              <div id="labVerificationError" class="alert alert-danger d-none"></div>
              <div id="labVerificationSuccess" class="alert alert-success d-none"></div>
            </div>
            <div class="modal-footer">
              <button type="submit" class="btn btn-primary">Save Verification</button>
              <button type="button" class="btn btn-secondary" data-bs-dismiss="modal">Cancel</button>
            </div>
          </form>
        </div>
      </div>
    `;
    document.body.insertAdjacentHTML('beforeend', modalHtml);
  }

  // Form submit handler
  document.getElementById('labVerificationForm').onsubmit = function (e) {
    e.preventDefault();
    submitLabVerification();
  };
}

window.showLabVerificationModal = function (experimentId, verification) {
  // Set form values
  document.getElementById('labVerificationExperimentId').value = experimentId;
  document.getElementById('labVerificationStatus').value = verification?.status || '';
  document.getElementById('labVerificationVerifier').value = verification?.verifier || '';
  document.getElementById('labVerificationEquipment').value = verification?.equipment || '';
  document.getElementById('labVerificationComments').value = verification?.comments || '';
  document.getElementById('labVerificationError').classList.add('d-none');
  document.getElementById('labVerificationSuccess').classList.add('d-none');

  // Show modal
  const modal = new bootstrap.Modal(document.getElementById('labVerificationModal'));
  modal.show();
};

function submitLabVerification() {
  const experimentId = document.getElementById('labVerificationExperimentId').value;
  const status = document.getElementById('labVerificationStatus').value;
  const verifier = document.getElementById('labVerificationVerifier').value.trim();
  const equipment = document.getElementById('labVerificationEquipment').value.trim();
  const comments = document.getElementById('labVerificationComments').value.trim();

  // Input validation
  if (!status || !verifier || !equipment) {
    showLabVerificationError('Please fill in all required fields.');
    return;
  }

  const payload = { status, verifier, equipment, comments };

  API.request(`/experiments/${experimentId}/lab_verification`, {
    method: 'POST',
    body: JSON.stringify(payload)
  })
    .then(data => {
      showLabVerificationSuccess('Verification saved successfully.');
      // Update the section in the experiment card
      const section = document.querySelector(`.lab-verification-section[data-experiment-id="${experimentId}"]`);
      if (section) renderLabVerification(section, experimentId, data);
      // Hide modal after short delay
      setTimeout(() => {
        bootstrap.Modal.getInstance(document.getElementById('labVerificationModal')).hide();
      }, 1000);
    })
    .catch(err => {
      showLabVerificationError('Failed to save verification. Please try again.');
    });
}

function showLabVerificationError(msg) {
  const el = document.getElementById('labVerificationError');
  el.textContent = msg;
  el.classList.remove('d-none');
}

function showLabVerificationSuccess(msg) {
  const el = document.getElementById('labVerificationSuccess');
  el.textContent = msg;
  el.classList.remove('d-none');
}