/**
 * Unit tests for lab verification UI functionality
 * 
 * This file contains tests for the lab verification UI components,
 * including form validation, API interactions, and DOM manipulation.
 */

// Mock dependencies
const { JSDOM } = require('jsdom');
const sinon = require('sinon');
const { expect } = require('chai');

describe('Lab Verification UI', () => {
  let window, document, API, bootstrap;
  
  // Set up DOM environment before each test
  beforeEach(() => {
    // Create a DOM environment
    const dom = new JSDOM(`
      <!DOCTYPE html>
      <html>
        <body>
          <div class="lab-verification-section" data-experiment-id="test-experiment-id"></div>
          <div id="labVerificationModal">
            <form id="labVerificationForm">
              <input type="hidden" id="labVerificationExperimentId" value="" />
              <select id="labVerificationStatus">
                <option value="">Select status</option>
                <option value="Verified">Verified</option>
                <option value="Pending">Pending</option>
                <option value="Rejected">Rejected</option>
              </select>
              <input type="text" id="labVerificationVerifier" value="" />
              <input type="text" id="labVerificationEquipment" value="" />
              <textarea id="labVerificationComments"></textarea>
              <div id="labVerificationError" class="d-none"></div>
              <div id="labVerificationSuccess" class="d-none"></div>
            </form>
          </div>
        </body>
      </html>
    `, { url: 'http://localhost/' });
    
    // Set up global variables
    window = dom.window;
    document = window.document;
    global.window = window;
    global.document = document;
    
    // Mock Bootstrap
    bootstrap = {
      Modal: class {
        constructor() {
          this.shown = false;
        }
        
        show() {
          this.shown = true;
        }
        
        hide() {
          this.shown = false;
        }
        
        static getInstance() {
          return {
            hide: sinon.spy()
          };
        }
      }
    };
    global.bootstrap = bootstrap;
    
    // Mock API
    API = {
      request: sinon.stub()
    };
    global.API = API;
    
    // Load the lab verification script
    require('../../../static/js/lab-verification');
  });
  
  afterEach(() => {
    // Clean up
    sinon.restore();
    delete global.window;
    delete global.document;
    delete global.bootstrap;
    delete global.API;
  });
  
  describe('loadLabVerification', () => {
    it('should load verification data successfully', async () => {
      const section = document.querySelector('.lab-verification-section');
      const experimentId = 'test-experiment-id';
      const verificationData = {
        id: 'test-verification-id',
        status: 'Verified',
        verifier: 'Test Verifier',
        equipment: 'Test Equipment',
        comments: 'Test Comments'
      };
      
      // Mock API response
      API.request.resolves(verificationData);
      
      // Call the function
      await window.loadLabVerification(section, experimentId);
      
      // Verify API was called correctly
      expect(API.request.calledOnce).to.be.true;
      expect(API.request.firstCall.args[0]).to.equal(`/experiments/${experimentId}/lab_verification`);
      
      // Verify the section was updated with verification data
      expect(section.innerHTML).to.include('Verified');
      expect(section.innerHTML).to.include('Test Verifier');
      expect(section.innerHTML).to.include('Test Equipment');
      expect(section.innerHTML).to.include('Test Comments');
    });
    
    it('should handle verification not found', async () => {
      const section = document.querySelector('.lab-verification-section');
      const experimentId = 'test-experiment-id';
      
      // Mock API response (error)
      API.request.rejects(new Error('Not found'));
      
      // Call the function
      await window.loadLabVerification(section, experimentId);
      
      // Verify API was called correctly
      expect(API.request.calledOnce).to.be.true;
      
      // Verify the section shows "No verification" message
      expect(section.innerHTML).to.include('No lab verification recorded');
      expect(section.innerHTML).to.include('Record Verification');
    });
  });
  
  describe('renderLabVerification', () => {
    it('should render verification data correctly', () => {
      const section = document.querySelector('.lab-verification-section');
      const experimentId = 'test-experiment-id';
      const verification = {
        status: 'Verified',
        verifier: 'Test Verifier',
        equipment: 'Test Equipment',
        comments: 'Test Comments'
      };
      
      // Call the function
      window.renderLabVerification(section, experimentId, verification);
      
      // Verify the section was updated with verification data
      expect(section.innerHTML).to.include('Verified');
      expect(section.innerHTML).to.include('Test Verifier');
      expect(section.innerHTML).to.include('Test Equipment');
      expect(section.innerHTML).to.include('Test Comments');
      expect(section.innerHTML).to.include('Update Verification');
    });
    
    it('should handle different verification statuses', () => {
      const section = document.querySelector('.lab-verification-section');
      const experimentId = 'test-experiment-id';
      
      // Test with 'Pending' status
      const pendingVerification = {
        status: 'Pending',
        verifier: 'Test Verifier',
        equipment: 'Test Equipment'
      };
      
      window.renderLabVerification(section, experimentId, pendingVerification);
      expect(section.innerHTML).to.include('Pending');
      expect(section.innerHTML).to.include('bg-warning');
      
      // Test with 'Rejected' status
      const rejectedVerification = {
        status: 'Rejected',
        verifier: 'Test Verifier',
        equipment: 'Test Equipment'
      };
      
      window.renderLabVerification(section, experimentId, rejectedVerification);
      expect(section.innerHTML).to.include('Rejected');
      expect(section.innerHTML).to.include('bg-danger');
    });
  });
  
  describe('showLabVerificationModal', () => {
    it('should populate modal with verification data', () => {
      const experimentId = 'test-experiment-id';
      const verification = {
        status: 'Verified',
        verifier: 'Test Verifier',
        equipment: 'Test Equipment',
        comments: 'Test Comments'
      };
      
      // Call the function
      window.showLabVerificationModal(experimentId, verification);
      
      // Verify form fields were populated
      expect(document.getElementById('labVerificationExperimentId').value).to.equal(experimentId);
      expect(document.getElementById('labVerificationStatus').value).to.equal(verification.status);
      expect(document.getElementById('labVerificationVerifier').value).to.equal(verification.verifier);
      expect(document.getElementById('labVerificationEquipment').value).to.equal(verification.equipment);
      expect(document.getElementById('labVerificationComments').value).to.equal(verification.comments);
    });
    
    it('should handle empty verification data', () => {
      const experimentId = 'test-experiment-id';
      
      // Call the function with no verification data
      window.showLabVerificationModal(experimentId);
      
      // Verify form fields were populated with defaults
      expect(document.getElementById('labVerificationExperimentId').value).to.equal(experimentId);
      expect(document.getElementById('labVerificationStatus').value).to.equal('');
      expect(document.getElementById('labVerificationVerifier').value).to.equal('');
      expect(document.getElementById('labVerificationEquipment').value).to.equal('');
      expect(document.getElementById('labVerificationComments').value).to.equal('');
    });
  });
  
  describe('submitLabVerification', () => {
    it('should submit verification data successfully', async () => {
      // Set up form data
      document.getElementById('labVerificationExperimentId').value = 'test-experiment-id';
      document.getElementById('labVerificationStatus').value = 'Verified';
      document.getElementById('labVerificationVerifier').value = 'Test Verifier';
      document.getElementById('labVerificationEquipment').value = 'Test Equipment';
      document.getElementById('labVerificationComments').value = 'Test Comments';
      
      // Mock API response
      const responseData = {
        id: 'test-verification-id',
        status: 'Verified',
        verifier: 'Test Verifier',
        equipment: 'Test Equipment',
        comments: 'Test Comments'
      };
      API.request.resolves(responseData);
      
      // Call the function
      await window.submitLabVerification();
      
      // Verify API was called correctly
      expect(API.request.calledOnce).to.be.true;
      expect(API.request.firstCall.args[0]).to.equal('/experiments/test-experiment-id/lab_verification');
      expect(API.request.firstCall.args[1].method).to.equal('POST');
      
      // Verify payload
      const payload = JSON.parse(API.request.firstCall.args[1].body);
      expect(payload.status).to.equal('Verified');
      expect(payload.verifier).to.equal('Test Verifier');
      expect(payload.equipment).to.equal('Test Equipment');
      expect(payload.comments).to.equal('Test Comments');
      
      // Verify success message was shown
      expect(document.getElementById('labVerificationSuccess').classList.contains('d-none')).to.be.false;
    });
    
    it('should validate required fields', async () => {
      // Set up form data with missing fields
      document.getElementById('labVerificationExperimentId').value = 'test-experiment-id';
      document.getElementById('labVerificationStatus').value = 'Verified';
      document.getElementById('labVerificationVerifier').value = ''; // Missing verifier
      document.getElementById('labVerificationEquipment').value = 'Test Equipment';
      
      // Call the function
      await window.submitLabVerification();
      
      // Verify API was not called
      expect(API.request.called).to.be.false;
      
      // Verify error message was shown
      expect(document.getElementById('labVerificationError').classList.contains('d-none')).to.be.false;
      expect(document.getElementById('labVerificationError').textContent).to.include('required fields');
    });
    
    it('should handle API errors', async () => {
      // Set up form data
      document.getElementById('labVerificationExperimentId').value = 'test-experiment-id';
      document.getElementById('labVerificationStatus').value = 'Verified';
      document.getElementById('labVerificationVerifier').value = 'Test Verifier';
      document.getElementById('labVerificationEquipment').value = 'Test Equipment';
      
      // Mock API error
      API.request.rejects(new Error('API Error'));
      
      // Call the function
      await window.submitLabVerification();
      
      // Verify API was called
      expect(API.request.calledOnce).to.be.true;
      
      // Verify error message was shown
      expect(document.getElementById('labVerificationError').classList.contains('d-none')).to.be.false;
      expect(document.getElementById('labVerificationError').textContent).to.include('Failed to save verification');
    });
  });
});