#!/usr/bin/env python3
"""
CryoProtect Analyzer - PubChem Cryoprotectant Data Importer (Supabase Version with Enhanced Progress Tracking and MCP Integration)

This script retrieves cryoprotectant data from PubChem, filters molecules based on
predefined criteria, scores them, and stores the results in a Supabase database using MCP.

Features:
- Batch processing with parameterizable batch size
- Robust checkpointing for resumable operation
- CLI arguments for batch size, checkpoint path, resume/reset
- Efficient Supabase bulk inserts for molecules and properties via MCP
- Enhanced real-time progress tracking with web dashboard
- Detailed statistics (success/failure, time metrics, etc.)
- Error/skipped CID logging for review
- Fetches property type metadata once per run

Prerequisites:
- Python 3.6+ installed
- Supabase project with the CryoProtect schema applied
- supabase-py package installed (pip install supabase)
- python-dotenv package installed (pip install python-dotenv)
- flask package installed (pip install flask)
"""

import os
import time
import requests
import json
import argparse
import threading
import webbrowser
from datetime import datetime, timedelta
import logging
from dotenv import load_dotenv
import socket
import atexit
import signal
import sys
from pathlib import Path
from flask import Flask, render_template_string, jsonify

# Import MCP helper
from use_mcp_tool import execute_sql, get_project_id

# Load environment variables from .env file
load_dotenv()

# Ensure logs directory exists
Path("logs").mkdir(exist_ok=True)

# Set up logging
LOG_FILE = "logs/cryoprotectant_analysis.log"
SKIPPED_CID_LOG = "logs/skipped_cids.log"
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# Sunday-optimized parameters
# Parse the PUBCHEM_API_DELAY value, handling potential comments in the .env file
api_delay_str = os.getenv("PUBCHEM_API_DELAY", "0.15")
if "#" in api_delay_str:
    api_delay_str = api_delay_str.split("#")[0].strip()
PUBCHEM_API_DELAY = float(api_delay_str)  # Sunday-optimized delay
DEFAULT_BATCH_SIZE = int(os.getenv("PUBCHEM_BATCH_SIZE", "100"))   # Sunday-optimized batch size

# CID File for Production
# IMPORTANT: For production imports, this file must be generated using generate_production_cid_list.py.
# It should contain the complete, scientifically accurate set of cryoprotectant CIDs and synonyms.
CID_FILE = "CID-Synonym-curated"  # Production-grade CID list for full database population

# Scoring Weights (Total = 200)
WEIGHTS = {
    "hydrogen_bonding": 50,
    "solubility_polarity": 40,
    "membrane_permeability": 40,
    "toxicity_biocompatibility": 30,
    "protein_stabilization": 20,
    "stability_reactivity": 10,
    "environmental_safety": 10
}

# Stage 1: Core Filtering Criteria
CORE_CRITERIA = {
    "logP_range": (-5, 5),  # Relaxed for testing
    "mw_range": (0, 1000),  # Relaxed for testing
    "TPSA_range": (0, 200), # Relaxed for testing
    "functional_groups": [] # No functional group requirement for testing
}

# Progress tracking dashboard HTML template
DASHBOARD_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>CryoProtect Import Progress</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0-alpha1/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {
            padding: 20px;
            background-color: #f8f9fa;
        }
        .card {
            margin-bottom: 20px;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }
        .card-header {
            font-weight: bold;
            background-color: #e9ecef;
        }
        .progress {
            height: 25px;
        }
        .progress-bar {
            font-size: 14px;
            font-weight: bold;
        }
        .status-running {
            color: #0d6efd;
        }
        .status-paused {
            color: #fd7e14;
        }
        .status-completed {
            color: #198754;
        }
        .status-error {
            color: #dc3545;
        }
        .refresh-text {
            font-size: 12px;
            color: #6c757d;
        }
        .stats-value {
            font-weight: bold;
            font-size: 18px;
        }
        .stats-label {
            font-size: 14px;
            color: #6c757d;
        }
        #log-container {
            max-height: 200px;
            overflow-y: auto;
            background-color: #212529;
            color: #f8f9fa;
            padding: 10px;
            border-radius: 5px;
            font-family: monospace;
        }
        .log-entry {
            margin: 0;
            padding: 2px 0;
        }
        .log-info {
            color: #0dcaf0;
        }
        .log-warning {
            color: #ffc107;
        }
        .log-error {
            color: #dc3545;
        }
        .log-success {
            color: #20c997;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1 class="mb-4">CryoProtect PubChem Import Progress</h1>
        
        <div class="row">
            <div class="col-md-8">
                <div class="card">
                    <div class="card-header d-flex justify-content-between align-items-center">
                        Overall Progress
                        <span id="status" class="status-running">Running</span>
                    </div>
                    <div class="card-body">
                        <div class="progress mb-3">
                            <div id="progress-bar" class="progress-bar progress-bar-striped progress-bar-animated" 
                                 role="progressbar" style="width: 0%;" 
                                 aria-valuenow="0" aria-valuemin="0" aria-valuemax="100">0%</div>
                        </div>
                        <div class="d-flex justify-content-between">
                            <div><span id="processed-count">0</span>/<span id="total-count">0</span> compounds processed</div>
                            <div class="refresh-text">Auto-refreshes every 2 seconds</div>
                        </div>
                    </div>
                </div>
                
                <div class="card">
                    <div class="card-header">Current Batch Information</div>
                    <div class="card-body">
                        <div class="row">
                            <div class="col-md-4 text-center">
                                <div class="stats-value" id="current-batch">0</div>
                                <div class="stats-label">Current Batch</div>
                            </div>
                            <div class="col-md-4 text-center">
                                <div class="stats-value" id="total-batches">0</div>
                                <div class="stats-label">Total Batches</div>
                            </div>
                            <div class="col-md-4 text-center">
                                <div class="stats-value" id="compounds-in-batch">0</div>
                                <div class="stats-label">Compounds in Batch</div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
            
            <div class="col-md-4">
                <div class="card">
                    <div class="card-header">Time Metrics</div>
                    <div class="card-body">
                        <div class="row mb-2">
                            <div class="col-6 stats-label">Elapsed Time:</div>
                            <div class="col-6 stats-value" id="elapsed-time">00:00:00</div>
                        </div>
                        <div class="row mb-2">
                            <div class="col-6 stats-label">Estimated Remaining:</div>
                            <div class="col-6 stats-value" id="eta">00:00:00</div>
                        </div>
                        <div class="row">
                            <div class="col-6 stats-label">Avg. Time per Batch:</div>
                            <div class="col-6 stats-value" id="avg-time-batch">00:00:00</div>
                        </div>
                    </div>
                </div>
                
                <div class="card">
                    <div class="card-header">Success/Failure Statistics</div>
                    <div class="card-body">
                        <div class="row mb-2">
                            <div class="col-8 stats-label">Successful Imports:</div>
                            <div class="col-4 stats-value" id="successful-imports">0</div>
                        </div>
                        <div class="row mb-2">
                            <div class="col-8 stats-label">Skipped Compounds:</div>
                            <div class="col-4 stats-value" id="skipped-compounds">0</div>
                        </div>
                        <div class="row">
                            <div class="col-8 stats-label">Errors:</div>
                            <div class="col-4 stats-value" id="errors">0</div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="card mt-3">
            <div class="card-header">Recent Log Messages</div>
            <div class="card-body">
                <div id="log-container">
                    <p class="log-entry log-info">Waiting for log messages...</p>
                </div>
            </div>
        </div>
    </div>
    
    <script>
        // Function to fetch progress data and update the dashboard
        function updateDashboard() {
            fetch('/progress')
                .then(response => response.json())
                .then(data => {
                    // Update progress bar
                    const progressPercent = data.progress_percentage;
                    document.getElementById('progress-bar').style.width = progressPercent + '%';
                    document.getElementById('progress-bar').setAttribute('aria-valuenow', progressPercent);
                    document.getElementById('progress-bar').textContent = progressPercent + '%';
                    
                    // Update counts
                    document.getElementById('processed-count').textContent = data.total_processed;
                    document.getElementById('total-count').textContent = data.total_compounds;
                    
                    // Update batch information
                    document.getElementById('current-batch').textContent = data.current_batch;
                    document.getElementById('total-batches').textContent = data.total_batches;
                    document.getElementById('compounds-in-batch').textContent = data.compounds_in_current_batch;
                    
                    // Update time metrics
                    document.getElementById('elapsed-time').textContent = data.elapsed_time;
                    document.getElementById('eta').textContent = data.estimated_time_remaining;
                    document.getElementById('avg-time-batch').textContent = data.avg_time_per_batch;
                    
                    // Update success/failure statistics
                    document.getElementById('successful-imports').textContent = data.total_imported;
                    document.getElementById('skipped-compounds').textContent = data.total_skipped;
                    document.getElementById('errors').textContent = data.total_errors;
                    
                    // Update status
                    const statusElement = document.getElementById('status');
                    statusElement.textContent = data.status;
                    statusElement.className = 'status-' + data.status.toLowerCase();
                    
                    // Update log messages
                    const logContainer = document.getElementById('log-container');
                    logContainer.innerHTML = '';
                    data.recent_logs.forEach(log => {
                        const logEntry = document.createElement('p');
                        logEntry.className = 'log-entry';
                        
                        // Add appropriate class based on log level
                        if (log.includes('[ERROR]')) {
                            logEntry.classList.add('log-error');
                        } else if (log.includes('[WARNING]')) {
                            logEntry.classList.add('log-warning');
                        } else if (log.includes('SUCCESS')) {
                            logEntry.classList.add('log-success');
                        } else {
                            logEntry.classList.add('log-info');
                        }
                        
                        logEntry.textContent = log;
                        logContainer.appendChild(logEntry);
                    });
                    
                    // Scroll to bottom of log container
                    logContainer.scrollTop = logContainer.scrollHeight;
                })
                .catch(error => {
                    console.error('Error fetching progress data:', error);
                });
        }
        
        // Update dashboard immediately and then every 2 seconds
        updateDashboard();
        setInterval(updateDashboard, 2000);
    </script>
</body>
</html>
"""
class ProgressTracker:
    """Class to track and manage progress statistics for the PubChem import process."""
    
    def __init__(self, total_compounds, batch_size, checkpoint_path):
        self.total_compounds = total_compounds
        self.batch_size = batch_size
        self.total_batches = (total_compounds + batch_size - 1) // batch_size
        self.checkpoint_path = checkpoint_path
        
        # Progress statistics
        self.start_time = time.time()
        self.current_batch = 0
        self.total_processed = 0
        self.total_imported = 0
        self.total_skipped = 0
        self.total_errors = 0
        self.compounds_in_current_batch = 0
        self.batch_start_time = time.time()
        self.batch_times = []
        self.status = "Running"
        
        # Recent log messages (circular buffer)
        self.max_log_entries = 50
        self.recent_logs = []
        
        # Load from checkpoint if exists
        self._load_from_checkpoint()
        
        # Register signal handlers for graceful shutdown
        signal.signal(signal.SIGINT, self._handle_exit)
        signal.signal(signal.SIGTERM, self._handle_exit)
        
        # Register exit handler
        atexit.register(self._handle_exit)
    
    def _handle_exit(self, *args):
        """Handle script exit (save checkpoint and set status)."""
        if self.status == "Running":
            self.status = "Paused"
            self.save_checkpoint()
        sys.exit(0)
    
    def _load_from_checkpoint(self):
        """Load progress data from checkpoint file if it exists."""
        if os.path.exists(self.checkpoint_path):
            try:
                with open(self.checkpoint_path, "r") as f:
                    checkpoint = json.load(f)
                
                # Basic checkpoint data
                self.current_batch = checkpoint.get("last_completed_batch", 0) + 1
                self.total_processed = checkpoint.get("total_processed", 0)
                self.total_imported = checkpoint.get("total_imported", 0)
                
                # Enhanced progress data
                self.total_skipped = checkpoint.get("total_skipped", 0)
                self.total_errors = checkpoint.get("total_errors", 0)
                self.batch_times = checkpoint.get("batch_times", [])
                
                # Adjust start time to maintain accurate elapsed time
                if "elapsed_seconds" in checkpoint:
                    self.start_time = time.time() - checkpoint["elapsed_seconds"]
                
                logger.info(f"Loaded progress from checkpoint: {self.total_processed}/{self.total_compounds} compounds processed")
            except Exception as e:
                logger.error(f"Error loading checkpoint: {str(e)}")
    
    def save_checkpoint(self):
        """Save progress data to checkpoint file."""
        checkpoint_data = {
            "last_completed_batch": self.current_batch - 1 if self.current_batch > 0 else 0,
            "total_processed": self.total_processed,
            "total_imported": self.total_imported,
            "total_skipped": self.total_skipped,
            "total_errors": self.total_errors,
            "batch_times": self.batch_times[-20:],  # Keep only the last 20 batch times
            "elapsed_seconds": time.time() - self.start_time,
            "timestamp": datetime.now().isoformat(),
            "status": self.status
        }
        
        with open(self.checkpoint_path, "w") as f:
            json.dump(checkpoint_data, f)
    
    def start_batch(self, batch_num, batch_size):
        """Record the start of a new batch."""
        self.current_batch = batch_num + 1  # 1-indexed for display
        self.compounds_in_current_batch = batch_size
        self.batch_start_time = time.time()
        self.add_log_message(f"Starting batch {self.current_batch}/{self.total_batches} with {batch_size} compounds")
    
    def end_batch(self, processed, imported, skipped, errors):
        """Record the end of a batch."""
        batch_time = time.time() - self.batch_start_time
        self.batch_times.append(batch_time)
        
        self.total_processed += processed
        self.total_imported += imported
        self.total_skipped += skipped
        self.total_errors += errors
        
        # Calculate progress metrics
        progress_percent = round((self.total_processed / self.total_compounds) * 100, 1)
        eta = self.estimate_time_remaining()
        
        self.add_log_message(
            f"Completed batch {self.current_batch}/{self.total_batches}: "
            f"{processed} processed, {imported} imported, {skipped} skipped, {errors} errors. "
            f"Progress: {progress_percent}%, ETA: {eta}"
        )
        
        # Save checkpoint after each batch
        self.save_checkpoint()
    
    def add_error(self, error_message):
        """Record an error."""
        self.total_errors += 1
        self.add_log_message(f"ERROR: {error_message}")
    
    def add_skipped(self, cid, reason):
        """Record a skipped compound."""
        self.total_skipped += 1
        self.add_log_message(f"SKIPPED: CID {cid} - {reason}")
    
    def add_log_message(self, message):
        """Add a log message to the recent logs buffer."""
        timestamp = datetime.now().strftime("%H:%M:%S")
        log_entry = f"{timestamp} - {message}"
        
        # Add to circular buffer
        self.recent_logs.append(log_entry)
        if len(self.recent_logs) > self.max_log_entries:
            self.recent_logs.pop(0)
    
    def estimate_time_remaining(self):
        """Estimate time remaining based on average batch processing time."""
        if not self.batch_times:
            return "Unknown"
        
        # Use the last 10 batch times for a more accurate recent average
        recent_batch_times = self.batch_times[-10:] if len(self.batch_times) >= 10 else self.batch_times
        avg_batch_time = sum(recent_batch_times) / len(recent_batch_times)
        
        remaining_batches = self.total_batches - self.current_batch
        eta_seconds = avg_batch_time * remaining_batches
        
        return str(timedelta(seconds=int(eta_seconds)))
    
    def get_avg_time_per_batch(self):
        """Get the average time per batch."""
        if not self.batch_times:
            return "00:00:00"
        
        avg_seconds = sum(self.batch_times) / len(self.batch_times)
        return str(timedelta(seconds=int(avg_seconds)))
    
    def get_elapsed_time(self):
        """Get the total elapsed time."""
        elapsed_seconds = time.time() - self.start_time
        return str(timedelta(seconds=int(elapsed_seconds)))
    
    def get_progress_percentage(self):
        """Get the progress percentage."""
        if self.total_compounds == 0:
            return 0
        return round((self.total_processed / self.total_compounds) * 100, 1)
    
    def set_status(self, status):
        """Set the current status of the import process."""
        self.status = status
        self.add_log_message(f"Status changed to: {status}")
        self.save_checkpoint()
    
    def get_progress_data(self):
        """Get all progress data as a dictionary for the dashboard."""
        return {
            "total_compounds": self.total_compounds,
            "total_processed": self.total_processed,
            "total_imported": self.total_imported,
            "total_skipped": self.total_skipped,
            "total_errors": self.total_errors,
            "current_batch": self.current_batch,
            "total_batches": self.total_batches,
            "compounds_in_current_batch": self.compounds_in_current_batch,
            "elapsed_time": self.get_elapsed_time(),
            "estimated_time_remaining": self.estimate_time_remaining(),
            "avg_time_per_batch": self.get_avg_time_per_batch(),
            "progress_percentage": self.get_progress_percentage(),
            "status": self.status,
            "recent_logs": self.recent_logs
        }

class DashboardServer:
    """Flask web server for the progress tracking dashboard."""
    
    def __init__(self, progress_tracker, host='127.0.0.1', port=None):
        self.progress_tracker = progress_tracker
        self.host = host
        self.port = port or self._find_free_port()
        self.app = Flask(__name__)
        self.server_thread = None
        self.url = f"http://{self.host}:{self.port}"
        
        # Set up routes
        self._setup_routes()
    
    def _find_free_port(self):
        """Find a free port to run the dashboard server."""
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.bind(('', 0))
            return s.getsockname()[1]
    
    def _setup_routes(self):
        """Set up Flask routes for the dashboard."""
        @self.app.route('/')
        def index():
            return render_template_string(DASHBOARD_TEMPLATE)
        
        @self.app.route('/progress')
        def progress():
            return jsonify(self.progress_tracker.get_progress_data())
    
    def start(self):
        """Start the dashboard server in a separate thread."""
        def run_server():
            self.app.run(host=self.host, port=self.port, debug=False, use_reloader=False)
        
        self.server_thread = threading.Thread(target=run_server)
        self.server_thread.daemon = True
        self.server_thread.start()
        
        logger.info(f"Dashboard server started at {self.url}")
        return self.url
    
    def open_browser(self):
        """Open the dashboard in a web browser."""
        webbrowser.open(self.url)

def get_cid_list():
    """Retrieve all CIDs from the downloaded PubChem CID list."""
    if not os.path.exists(CID_FILE):
        logger.warning(f"WARNING: CID file '{CID_FILE}' not found. Download it first.")
        return []

    with open(CID_FILE, "r") as file:
        cids = [int(line.strip().split("\t")[0]) for line in file if line.strip().split("\t")[0].isdigit()]

    logger.info(f"SUCCESS: Loaded {len(cids)} CIDs from PubChem's CID list.")
    return cids

def get_molecule_properties(cid, progress_tracker):
    """Fetch molecular properties and names from PubChem."""
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/"
        "MolecularFormula,MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,"
        "IsomericSMILES,InChI,InChIKey,IUPACName,Title/JSON"
    )
    try:
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            properties = data["PropertyTable"]["Properties"][0]
            return {
                "CID": cid,
                "Molecular Formula": properties.get("MolecularFormula"),
                "Molecular Weight": properties.get("MolecularWeight"),
                "LogP": properties.get("XLogP"),
                "TPSA": properties.get("TPSA"),
                "H-Bond Donors": properties.get("HBondDonorCount"),
                "H-Bond Acceptors": properties.get("HBondAcceptorCount"),
                "SMILES": properties.get("IsomericSMILES"),
                "InChI": properties.get("InChI"),
                "InChIKey": properties.get("InChIKey"),
                "IUPACName": properties.get("IUPACName"),
                "Title": properties.get("Title"),
                "PubChem Link": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
            }
        error_msg = f"No molecular properties found for CID {cid}."
        progress_tracker.add_log_message(f"WARNING: {error_msg}")
        return {"CID": cid, "Error": "No data found"}
    except Exception as e:
        error_msg = f"Error fetching properties for CID {cid}: {str(e)}"
        progress_tracker.add_error(error_msg)
        return {"CID": cid, "Error": str(e)}

def get_additional_properties(cid, progress_tracker):
    """Fetch toxicity, stability, and environmental impact from PubChem."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON"
    try:
        response = requests.get(url)
        
        if response.status_code == 200:
            data = response.json()
            props = data["PC_Compounds"][0]["props"]

            extra_data = {
                "Toxicity": None,
                "Stability": None,
                "Environmental Safety": None
            }

            for prop in props:
                label = prop.get("urn", {}).get("label", "")
                value = prop.get("value", {}).get("sval", "")

                if "Toxicity" in label:
                    extra_data["Toxicity"] = value
                elif "Stability" in label:
                    extra_data["Stability"] = value
                elif "Environmental" in label:
                    extra_data["Environmental Safety"] = value

            return extra_data
        
        error_msg = f"No additional properties found for CID {cid}."
        progress_tracker.add_log_message(f"WARNING: {error_msg}")
        return {"Error": "No additional data found"}
    except Exception as e:
        error_msg = f"Error fetching additional properties for CID {cid}: {str(e)}"
        progress_tracker.add_error(error_msg)
        return {"Error": str(e)}

def filter_molecule(molecule, progress_tracker):
    """Initial filtering based on core cryoprotectant properties, ensuring numerical values."""
    if "Error" in molecule:
        progress_tracker.add_skipped(molecule['CID'], "Error in molecule data")
        return False

    # Check for required fields for Supabase schema
    if not molecule.get("SMILES"):
        progress_tracker.add_skipped(molecule['CID'], "Missing SMILES")
        return False
    if not molecule.get("InChI"):
        progress_tracker.add_skipped(molecule['CID'], "Missing InChI")
        return False
    if not molecule.get("InChIKey"):
        progress_tracker.add_skipped(molecule['CID'], "Missing InChIKey")
        return False
    if not molecule.get("Molecular Formula"):
        progress_tracker.add_skipped(molecule['CID'], "Missing Molecular Formula")
        return False

    try:
        mw = float(molecule["Molecular Weight"]) if molecule["Molecular Weight"] else None
        logp = float(molecule["LogP"]) if molecule["LogP"] else None
        tpsa = float(molecule["TPSA"]) if molecule["TPSA"] else None
    except (ValueError, TypeError):
        progress_tracker.add_skipped(molecule['CID'], "Invalid numerical values")
        return False

    smiles = molecule["SMILES"]

    if mw is None or not (CORE_CRITERIA["mw_range"][0] <= mw <= CORE_CRITERIA["mw_range"][1]):
        progress_tracker.add_skipped(molecule['CID'], f"Molecular weight {mw} outside range {CORE_CRITERIA['mw_range']}")
        return False
    if logp is None or not (CORE_CRITERIA["logP_range"][0] <= logp <= CORE_CRITERIA["logP_range"][1]):
        progress_tracker.add_skipped(molecule['CID'], f"LogP {logp} outside range {CORE_CRITERIA['logP_range']}")
        return False
    if tpsa is None or not (CORE_CRITERIA["TPSA_range"][0] <= tpsa <= CORE_CRITERIA["TPSA_range"][1]):
        progress_tracker.add_skipped(molecule['CID'], f"TPSA {tpsa} outside range {CORE_CRITERIA['TPSA_range']}")
        return False
    if CORE_CRITERIA["functional_groups"]:
        if not any(group in smiles for group in CORE_CRITERIA["functional_groups"]):
            progress_tracker.add_skipped(molecule['CID'], "Missing required functional groups")
            return False

    return True

def score_molecule(molecule, extra_properties):
    """Compute final score out of 200 based on all properties."""
    score = 0

    # Hydrogen Bonding
    score += WEIGHTS["hydrogen_bonding"]

    # Solubility & Permeability
    score += WEIGHTS["solubility_polarity"]
    score += WEIGHTS["membrane_permeability"]

    # Toxicity & Stability
    if extra_properties.get("Toxicity"):
        score += WEIGHTS["toxicity_biocompatibility"]
    if extra_properties.get("Stability"):
        score += WEIGHTS["stability_reactivity"]

    # Environmental Safety
    if extra_properties.get("Environmental Safety"):
        score += WEIGHTS["environmental_safety"]

def get_project_id():
    """
    Get the Supabase project ID from environment variables.
    
    Returns:
        str: Project ID or None if not found
    """
    try:
        # Get project ID from environment variable
        project_id = os.getenv("SUPABASE_PROJECT_ID")
        
        if project_id:
            logger.info(f"Using Supabase project ID from environment: {project_id}")
            return project_id
        
        # If not found in environment, extract from URL
        supabase_url = os.getenv("SUPABASE_URL", "")
        if supabase_url:
            # Extract project ID from URL (format: https://[project_id].supabase.co)
            parts = supabase_url.split(".")
            if len(parts) >= 3 and "supabase" in parts[1]:
                project_id = parts[0].replace("https://", "")
                logger.info(f"Extracted project ID from URL: {project_id}")
                return project_id
        
        logger.warning("No project ID found in environment variables")
        return None
    
    except Exception as e:
        logger.error(f"Error getting project ID: {str(e)}")
        return None
    return score

def insert_molecule_mcp(molecule, extra_properties, project_id, progress_tracker):
    """Insert a molecule into the database using MCP."""
    try:
        # Calculate the cryoprotectant score
        cryo_score = score_molecule(molecule, extra_properties)
        
        # Prepare the SQL query for insertion
        sql = f"""
        INSERT INTO molecules (
            cid, 
            name, 
            smiles, 
            molecular_weight, 
            formula, 
            logp, 
            tpsa, 
            hbond_donors, 
            hbond_acceptors, 
            inchi, 
            inchi_key, 
            iupac_name, 
            data_source, 
            cryo_score, 
            created_at
        ) VALUES (
            {molecule['CID']}, 
            '{molecule.get('Title', '').replace("'", "''")}', 
            '{molecule.get('SMILES', '').replace("'", "''")}', 
            {molecule.get('Molecular Weight', 'NULL')}, 
            '{molecule.get('Molecular Formula', '').replace("'", "''")}', 
            {molecule.get('LogP', 'NULL')}, 
            {molecule.get('TPSA', 'NULL')}, 
            {molecule.get('H-Bond Donors', 'NULL')}, 
            {molecule.get('H-Bond Acceptors', 'NULL')}, 
            '{molecule.get('InChI', '').replace("'", "''")}', 
            '{molecule.get('InChIKey', '').replace("'", "''")}', 
            '{molecule.get('IUPACName', '').replace("'", "''")}', 
            'PubChem', 
            {cryo_score}, 
            NOW()
        ) ON CONFLICT (cid) DO UPDATE SET
            name = EXCLUDED.name,
            smiles = EXCLUDED.smiles,
            molecular_weight = EXCLUDED.molecular_weight,
            formula = EXCLUDED.formula,
            logp = EXCLUDED.logp,
            tpsa = EXCLUDED.tpsa,
            hbond_donors = EXCLUDED.hbond_donors,
            hbond_acceptors = EXCLUDED.hbond_acceptors,
            inchi = EXCLUDED.inchi,
            inchi_key = EXCLUDED.inchi_key,
            iupac_name = EXCLUDED.iupac_name,
            data_source = EXCLUDED.data_source,
            cryo_score = EXCLUDED.cryo_score,
            updated_at = NOW()
        RETURNING id;
        """
        
        # Execute the SQL query using MCP
        result = execute_sql(sql, project_id)
        
        if isinstance(result, dict) and "error" in result:
            progress_tracker.add_error(f"Database error for CID {molecule['CID']}: {result['error']}")
            return False
        
        progress_tracker.add_log_message(f"SUCCESS: Inserted CID {molecule['CID']} ({molecule.get('Title', 'Unknown')}) with score {cryo_score}")
        return True
    
    except Exception as e:
        progress_tracker.add_error(f"Error inserting CID {molecule['CID']}: {str(e)}")
        return False

def process_batch(batch_cids, batch_num, project_id, progress_tracker):
    """Process a batch of CIDs."""
    processed = 0
    imported = 0
    skipped = 0
    errors = 0
    
    progress_tracker.start_batch(batch_num, len(batch_cids))
    
    for cid in batch_cids:
        try:
            # Fetch basic molecular properties
            molecule = get_molecule_properties(cid, progress_tracker)
            processed += 1
            
            # Apply initial filtering
            if not filter_molecule(molecule, progress_tracker):
                skipped += 1
                continue
            
            # Fetch additional properties for scoring
            extra_properties = get_additional_properties(cid, progress_tracker)
            
            # Insert into database using MCP
            if insert_molecule_mcp(molecule, extra_properties, project_id, progress_tracker):
                imported += 1
            else:
                errors += 1
            
            # Respect PubChem API rate limits
            time.sleep(PUBCHEM_API_DELAY)
            
        except Exception as e:
            progress_tracker.add_error(f"Unexpected error processing CID {cid}: {str(e)}")
            errors += 1
    
    progress_tracker.end_batch(processed, imported, skipped, errors)
    return processed, imported, skipped, errors

def main():
    """Main function to run the PubChem data import."""
    # Define global variable to be modified
    global PUBCHEM_API_DELAY
    
    parser = argparse.ArgumentParser(description="Import PubChem data into Supabase database using MCP.")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help=f"Batch size for processing (default: {DEFAULT_BATCH_SIZE})")
    parser.add_argument("--checkpoint", type=str, default="checkpoints/pubchem_import.json", help="Checkpoint file path")
    parser.add_argument("--reset", action="store_true", help="Reset progress and start from beginning")
    parser.add_argument("--no-dashboard", action="store_true", help="Disable web dashboard")
    parser.add_argument("--api-delay", type=float, default=PUBCHEM_API_DELAY, help=f"Delay between PubChem API calls in seconds (default: {PUBCHEM_API_DELAY})")
    parser.add_argument("--target", type=int, default=5000, help="Target number of compounds to import (default: 5000)")
    args = parser.parse_args()
    
    # Update global API delay
    PUBCHEM_API_DELAY = args.api_delay
    
    # Ensure checkpoints directory exists
    os.makedirs(os.path.dirname(args.checkpoint), exist_ok=True)
    
    # Reset checkpoint if requested
    if args.reset and os.path.exists(args.checkpoint):
        os.remove(args.checkpoint)
        logger.info(f"Reset progress: Removed checkpoint file {args.checkpoint}")
    
    # Get Supabase project ID
    project_id = get_project_id()
    if not project_id:
        logger.error("Failed to get Supabase project ID. Make sure you're logged in to Supabase CLI.")
        sys.exit(1)
    
    logger.info(f"Using Supabase project ID: {project_id}")
    
    # Get CID list
    cids = get_cid_list()
    if not cids:
        logger.error("No CIDs found. Make sure the CID file exists.")
        sys.exit(1)
    
    # Limit to target number if specified
    if args.target and args.target < len(cids):
        logger.info(f"Limiting to {args.target} compounds (out of {len(cids)} available)")
        cids = cids[:args.target]
    
    # Initialize progress tracker
    progress_tracker = ProgressTracker(len(cids), args.batch_size, args.checkpoint)
    
    # Start dashboard server if enabled
    if not args.no_dashboard:
        dashboard = DashboardServer(progress_tracker)
        dashboard_url = dashboard.start()
        logger.info(f"Dashboard available at: {dashboard_url}")
        dashboard.open_browser()
    
    # Process CIDs in batches
    total_processed = 0
    total_imported = 0
    total_skipped = 0
    total_errors = 0
    
    # Calculate batches
    num_batches = (len(cids) + args.batch_size - 1) // args.batch_size
    
    # Skip already processed batches based on checkpoint
    start_batch = progress_tracker.current_batch
    
    logger.info(f"Starting import from batch {start_batch + 1}/{num_batches}")
    
    try:
        for i in range(start_batch, num_batches):
            # Get batch of CIDs
            start_idx = i * args.batch_size
            end_idx = min(start_idx + args.batch_size, len(cids))
            batch_cids = cids[start_idx:end_idx]
            
            # Process batch
            processed, imported, skipped, errors = process_batch(batch_cids, i, project_id, progress_tracker)
            
            # Update totals
            total_processed += processed
            total_imported += imported
            total_skipped += skipped
            total_errors += errors
            
            # Log progress
            logger.info(f"Batch {i+1}/{num_batches} complete. Total: {total_processed} processed, {total_imported} imported, {total_skipped} skipped, {total_errors} errors")
        
        # Mark as completed
        progress_tracker.set_status("Completed")
        logger.info(f"Import completed. Total: {total_processed} processed, {total_imported} imported, {total_skipped} skipped, {total_errors} errors")
    
    except KeyboardInterrupt:
        logger.info("Import interrupted by user.")
        progress_tracker.set_status("Paused")
    
    except Exception as e:
        logger.error(f"Import failed: {str(e)}")
        progress_tracker.set_status("Error")
        progress_tracker.add_error(str(e))
    
    # Save final checkpoint
    progress_tracker.save_checkpoint()
    
    # Keep dashboard running if enabled
    if not args.no_dashboard:
        logger.info(f"Dashboard still available at: {dashboard_url}")
        logger.info("Press Ctrl+C to exit")
        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            logger.info("Exiting...")

if __name__ == "__main__":
    main()