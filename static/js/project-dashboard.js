/**
 * CryoProtect Analyzer - Project Dashboard Module
 * 
 * This module provides functions for managing projects and experiments in the user dashboard
 */

const ProjectDashboard = (function() {
    // Cache for dashboard data
    let dashboardData = null;
    let currentProject = null;
    let currentProjectId = null;
    let projectsOffset = 0;
    const projectsLimit = 10;
    
    /**
     * Initialize the project dashboard
     */
    function init() {
        // Set up event listeners
        document.getElementById('create-project-btn').addEventListener('click', showCreateProjectModal);
        document.getElementById('create-project-submit').addEventListener('click', createProject);
        document.getElementById('load-more-projects').addEventListener('click', loadMoreProjects);
        document.getElementById('project-search-btn').addEventListener('click', searchProjects);
        document.getElementById('project-search').addEventListener('keyup', function(event) {
            if (event.key === 'Enter') {
                searchProjects();
            }
        });
        
        document.getElementById('add-experiment-btn').addEventListener('click', showAddExperimentModal);
        document.getElementById('add-experiment-submit').addEventListener('click', addExperimentToProject);
        document.getElementById('update-project-btn').addEventListener('click', updateProject);
        document.getElementById('delete-project-btn').addEventListener('click', confirmDeleteProject);
        
        // Quick access buttons
        document.getElementById('create-experiment-btn').addEventListener('click', redirectToCreateExperiment);
        document.getElementById('import-data-btn').addEventListener('click', showImportDataModal);
        document.getElementById('export-results-btn').addEventListener('click', showExportResultsModal);
        document.getElementById('share-project-btn').addEventListener('click', showShareProjectModal);
        
        // Load dashboard data
        loadDashboardData();
    }
    
    /**
     * Load dashboard data from API
     */
    function loadDashboardData() {
        // Show loading indicators
        document.getElementById('projects-list').innerHTML = '<tr class="placeholder-row"><td colspan="4" class="text-center">Loading projects...</td></tr>';
        document.getElementById('activity-list').innerHTML = '<li class="list-group-item text-center">Loading activity...</li>';
        document.getElementById('saved-searches-list').innerHTML = '<li class="list-group-item text-center">Loading saved searches...</li>';
        
        // Fetch dashboard data
        API.get('/dashboard/projects')
            .then(data => {
                // Cache the data
                dashboardData = data;
                
                // Update dashboard stats
                updateDashboardStats(data);
                
                // Render projects
                renderProjects(data.projects.recent);
                
                // Render recent activity
                renderRecentActivity(data.recentActivity);
                
                // Render saved searches
                renderSavedSearches(data.savedSearches);
                
                // Render experiment distribution chart
                renderExperimentDistributionChart(data.experiments.by_project);
            })
            .catch(error => {
                console.error('Error loading dashboard data:', error);
                
                // Show error messages
                document.getElementById('projects-list').innerHTML = '<tr class="placeholder-row"><td colspan="4" class="text-center text-danger">Error loading projects</td></tr>';
                document.getElementById('activity-list').innerHTML = '<li class="list-group-item text-center text-danger">Error loading activity</li>';
                document.getElementById('saved-searches-list').innerHTML = '<li class="list-group-item text-center text-danger">Error loading saved searches</li>';
            });
    }
    
    /**
     * Update dashboard statistics
     * 
     * @param {Object} data - Dashboard data
     */
    function updateDashboardStats(data) {
        document.getElementById('total-projects').textContent = data.projects.total;
        document.getElementById('total-experiments').textContent = data.experiments.total;
        document.getElementById('shared-projects').textContent = data.projects.by_status.find(s => s.status === 'Shared')?.count || 0;
        document.getElementById('recent-activity').textContent = data.recentActivity.length;
        document.getElementById('projects-count').textContent = data.projects.total;
    }
    
    /**
     * Render projects table
     * 
     * @param {Array} projects - List of projects
     */
    function renderProjects(projects) {
        const projectsList = document.getElementById('projects-list');
        
        if (!projects || projects.length === 0) {
            projectsList.innerHTML = '<tr class="placeholder-row"><td colspan="4" class="text-center">No projects found</td></tr>';
            return;
        }
        
        let html = '';
        
        projects.forEach(project => {
            const date = new Date(project.updated_at);
            const formattedDate = date.toLocaleDateString('en-US', {
                year: 'numeric',
                month: 'short',
                day: 'numeric'
            });
            
            html += `
                <tr data-project-id="${project.id}">
                    <td>
                        <a href="#" class="project-name-link" data-project-id="${project.id}">
                            ${project.name}
                        </a>
                        ${project.is_public ? '<span class="badge bg-info ms-2">Public</span>' : ''}
                    </td>
                    <td>${project.experiment_count || 0}</td>
                    <td>${formattedDate}</td>
                    <td>
                        <div class="btn-group btn-group-sm">
                            <button type="button" class="btn btn-outline-primary view-project-btn" data-project-id="${project.id}">
                                <i class="bi bi-eye"></i>
                            </button>
                            <button type="button" class="btn btn-outline-secondary edit-project-btn" data-project-id="${project.id}">
                                <i class="bi bi-pencil"></i>
                            </button>
                            <button type="button" class="btn btn-outline-danger delete-project-btn" data-project-id="${project.id}">
                                <i class="bi bi-trash"></i>
                            </button>
                        </div>
                    </td>
                </tr>
            `;
        });
        
        projectsList.innerHTML = html;
        
        // Add event listeners to project actions
        document.querySelectorAll('.project-name-link, .view-project-btn').forEach(btn => {
            btn.addEventListener('click', function(event) {
                event.preventDefault();
                const projectId = this.getAttribute('data-project-id');
                viewProject(projectId);
            });
        });
        
        document.querySelectorAll('.edit-project-btn').forEach(btn => {
            btn.addEventListener('click', function() {
                const projectId = this.getAttribute('data-project-id');
                editProject(projectId);
            });
        });
        
        document.querySelectorAll('.delete-project-btn').forEach(btn => {
            btn.addEventListener('click', function() {
                const projectId = this.getAttribute('data-project-id');
                confirmDeleteProject(projectId);
            });
        });
    }
    
    /**
     * Render recent activity list
     * 
     * @param {Array} activities - List of recent activities
     */
    function renderRecentActivity(activities) {
        const activityList = document.getElementById('activity-list');
        
        if (!activities || activities.length === 0) {
            activityList.innerHTML = '<li class="list-group-item text-center">No recent activity</li>';
            return;
        }
        
        let html = '';
        
        activities.forEach(activity => {
            const date = new Date(activity.timestamp);
            const timeAgo = getTimeAgo(date);
            
            html += `
                <li class="list-group-item">
                    <div class="d-flex w-100 justify-content-between">
                        <h6 class="mb-1">${activity.title}</h6>
                        <small class="text-muted">${timeAgo}</small>
                    </div>
                    <p class="mb-1 small">${activity.description}</p>
                    <small class="text-muted">Project: ${activity.project_name}</small>
                </li>
            `;
        });
        
        activityList.innerHTML = html;
    }
    
    /**
     * Render saved searches list
     * 
     * @param {Array} searches - List of saved searches
     */
    function renderSavedSearches(searches) {
        const searchesList = document.getElementById('saved-searches-list');
        
        if (!searches || searches.length === 0) {
            searchesList.innerHTML = '<li class="list-group-item text-center">No saved searches</li>';
            return;
        }
        
        let html = '';
        
        searches.forEach(search => {
            html += `
                <li class="list-group-item">
                    <a href="#" class="saved-search-link" data-query="${search.query}">
                        ${search.name}
                    </a>
                </li>
            `;
        });
        
        searchesList.innerHTML = html;
        
        // Add event listeners to saved search links
        document.querySelectorAll('.saved-search-link').forEach(link => {
            link.addEventListener('click', function(event) {
                event.preventDefault();
                const query = this.getAttribute('data-query');
                applySavedSearch(query);
            });
        });
    }
    
    /**
     * Render experiment distribution chart
     * 
     * @param {Array} data - Experiment distribution data
     */
    function renderExperimentDistributionChart(data) {
        if (!data || data.length === 0) {
            document.getElementById('experiment-distribution-chart').innerHTML = '<div class="text-center py-5">No data available</div>';
            return;
        }
        
        const labels = data.map(item => item.project_name);
        const values = data.map(item => item.count);
        const colors = data.map((_, index) => COLORS[index % COLORS.length]);
        
        const chartData = {
            labels: labels,
            datasets: [{
                label: 'Experiments',
                data: values,
                backgroundColor: colors,
                borderColor: colors.map(color => color.replace('0.7', '1')),
                borderWidth: 1
            }]
        };
        
        const config = {
            type: 'bar',
            data: chartData,
            options: {
                responsive: true,
                maintainAspectRatio: false,
                plugins: {
                    legend: {
                        display: false
                    },
                    tooltip: {
                        callbacks: {
                            label: function(context) {
                                return `Experiments: ${context.raw}`;
                            }
                        }
                    }
                },
                scales: {
                    y: {
                        beginAtZero: true,
                        title: {
                            display: true,
                            text: 'Number of Experiments'
                        },
                        ticks: {
                            precision: 0
                        }
                    }
                }
            }
        };
        
        // Create chart
        new Chart(document.getElementById('experiment-distribution-chart'), config);
    }
    
    /**
     * Show create project modal
     */
    function showCreateProjectModal() {
        // Reset form
        document.getElementById('create-project-form').reset();
        
        // Show modal
        const modal = new bootstrap.Modal(document.getElementById('create-project-modal'));
        modal.show();
    }
    
    /**
     * Create a new project
     */
    function createProject() {
        const name = document.getElementById('project-name').value.trim();
        const description = document.getElementById('project-description').value.trim();
        const isPublic = document.getElementById('project-public').checked;
        
        if (!name) {
            alert('Project name is required');
            return;
        }
        
        const projectData = {
            name: name,
            description: description,
            is_public: isPublic
        };
        
        API.post('/projects', projectData)
            .then(project => {
                // Hide modal
                bootstrap.Modal.getInstance(document.getElementById('create-project-modal')).hide();
                
                // Reload dashboard data
                loadDashboardData();
                
                // Show success message
                alert('Project created successfully');
            })
            .catch(error => {
                console.error('Error creating project:', error);
                alert('Error creating project: ' + error.message);
            });
    }
    
    /**
     * Load more projects
     */
    function loadMoreProjects() {
        projectsOffset += projectsLimit;
        
        API.get(`/projects?limit=${projectsLimit}&offset=${projectsOffset}`)
            .then(projects => {
                if (projects.length === 0) {
                    document.getElementById('load-more-projects').disabled = true;
                    return;
                }
                
                // Append projects to the table
                const projectsList = document.getElementById('projects-list');
                let html = projectsList.innerHTML;
                
                projects.forEach(project => {
                    const date = new Date(project.updated_at);
                    const formattedDate = date.toLocaleDateString('en-US', {
                        year: 'numeric',
                        month: 'short',
                        day: 'numeric'
                    });
                    
                    html += `
                        <tr data-project-id="${project.id}">
                            <td>
                                <a href="#" class="project-name-link" data-project-id="${project.id}">
                                    ${project.name}
                                </a>
                                ${project.is_public ? '<span class="badge bg-info ms-2">Public</span>' : ''}
                            </td>
                            <td>${project.experiment_count || 0}</td>
                            <td>${formattedDate}</td>
                            <td>
                                <div class="btn-group btn-group-sm">
                                    <button type="button" class="btn btn-outline-primary view-project-btn" data-project-id="${project.id}">
                                        <i class="bi bi-eye"></i>
                                    </button>
                                    <button type="button" class="btn btn-outline-secondary edit-project-btn" data-project-id="${project.id}">
                                        <i class="bi bi-pencil"></i>
                                    </button>
                                    <button type="button" class="btn btn-outline-danger delete-project-btn" data-project-id="${project.id}">
                                        <i class="bi bi-trash"></i>
                                    </button>
                                </div>
                            </td>
                        </tr>
                    `;
                });
                
                projectsList.innerHTML = html;
                
                // Add event listeners to project actions
                document.querySelectorAll('.project-name-link, .view-project-btn').forEach(btn => {
                    btn.addEventListener('click', function(event) {
                        event.preventDefault();
                        const projectId = this.getAttribute('data-project-id');
                        viewProject(projectId);
                    });
                });
                
                document.querySelectorAll('.edit-project-btn').forEach(btn => {
                    btn.addEventListener('click', function() {
                        const projectId = this.getAttribute('data-project-id');
                        editProject(projectId);
                    });
                });
                
                document.querySelectorAll('.delete-project-btn').forEach(btn => {
                    btn.addEventListener('click', function() {
                        const projectId = this.getAttribute('data-project-id');
                        confirmDeleteProject(projectId);
                    });
                });
            })
            .catch(error => {
                console.error('Error loading more projects:', error);
                alert('Error loading more projects: ' + error.message);
            });
    }
    
    /**
     * Search projects
     */
    function searchProjects() {
        const searchTerm = document.getElementById('project-search').value.trim();
        
        if (!searchTerm) {
            // Reset to initial projects
            renderProjects(dashboardData.projects.recent);
            return;
        }
        
        // Filter projects by name (client-side filtering for simplicity)
        const filteredProjects = dashboardData.projects.recent.filter(project => 
            project.name.toLowerCase().includes(searchTerm.toLowerCase()) ||
            (project.description && project.description.toLowerCase().includes(searchTerm.toLowerCase()))
        );
        
        renderProjects(filteredProjects);
    }
    
    /**
     * View project details
     * 
     * @param {string} projectId - Project ID
     */
    function viewProject(projectId) {
        // Fetch project details
        API.get(`/projects/${projectId}`)
            .then(project => {
                currentProject = project;
                currentProjectId = projectId;
                
                // Update project details in modal
                document.getElementById('project-detail-name').textContent = project.name;
                document.getElementById('project-detail-description').textContent = project.description || 'No description';
                
                const createdDate = new Date(project.created_at);
                document.getElementById('project-detail-created').textContent = createdDate.toLocaleDateString('en-US', {
                    year: 'numeric',
                    month: 'short',
                    day: 'numeric'
                });
                
                // Set status badge
                const statusBadge = document.getElementById('project-detail-status');
                statusBadge.textContent = project.is_public ? 'Public' : 'Private';
                statusBadge.className = project.is_public ? 'badge bg-info' : 'badge bg-secondary';
                
                // Load project experiments
                loadProjectExperiments(projectId);
                
                // Load project activity
                loadProjectActivity(projectId);
                
                // Update form in settings tab
                document.getElementById('update-project-name').value = project.name;
                document.getElementById('update-project-description').value = project.description || '';
                document.getElementById('update-project-public').checked = project.is_public;
                
                // Show modal
                const modal = new bootstrap.Modal(document.getElementById('project-details-modal'));
                modal.show();
            })
            .catch(error => {
                console.error('Error fetching project details:', error);
                alert('Error fetching project details: ' + error.message);
            });
    }
    
    /**
     * Load project experiments
     * 
     * @param {string} projectId - Project ID
     */
    function loadProjectExperiments(projectId) {
        const experimentsList = document.getElementById('project-experiments-list');
        experimentsList.innerHTML = '<tr class="placeholder-row"><td colspan="5" class="text-center">Loading experiments...</td></tr>';
        
        API.get(`/projects/${projectId}/experiments`)
            .then(experiments => {
                if (!experiments || experiments.length === 0) {
                    experimentsList.innerHTML = '<tr class="placeholder-row"><td colspan="5" class="text-center">No experiments found</td></tr>';
                    return;
                }
                
                let html = '';
                
                experiments.forEach(experiment => {
                    const date = new Date(experiment.date_performed);
                    const formattedDate = date.toLocaleDateString('en-US', {
                        year: 'numeric',
                        month: 'short',
                        day: 'numeric'
                    });
                    
                    // Format value based on type
                    let value = '';
                    if (experiment.numeric_value !== null) {
                        value = experiment.numeric_value;
                    } else if (experiment.text_value !== null) {
                        value = experiment.text_value;
                    } else if (experiment.boolean_value !== null) {
                        value = experiment.boolean_value ? 'Yes' : 'No';
                    }
                    
                    html += `
                        <tr data-experiment-id="${experiment.id}">
                            <td>${experiment.name || 'Unnamed Experiment'}</td>
                            <td>${formattedDate}</td>
                            <td>${experiment.property_name}</td>
                            <td>${value}</td>
                            <td>
                                <div class="btn-group btn-group-sm">
                                    <button type="button" class="btn btn-outline-primary view-experiment-btn" data-experiment-id="${experiment.id}">
                                        <i class="bi bi-eye"></i>
                                    </button>
                                    <button type="button" class="btn btn-outline-danger remove-experiment-btn" data-experiment-id="${experiment.id}">
                                        <i class="bi bi-x-circle"></i>
                                    </button>
                                </div>
                            </td>
                        </tr>
                    `;
                });
                
                experimentsList.innerHTML = html;
                
                // Add event listeners to experiment actions
                document.querySelectorAll('.view-experiment-btn').forEach(btn => {
                    btn.addEventListener('click', function() {
                        const experimentId = this.getAttribute('data-experiment-id');
                        viewExperiment(experimentId);
                    });
                });
                
                document.querySelectorAll('.remove-experiment-btn').forEach(btn => {
                    btn.addEventListener('click', function() {
                        const experimentId = this.getAttribute('data-experiment-id');
                        confirmRemoveExperiment(projectId, experimentId);
                    });
                });
            })
            .catch(error => {
                console.error('Error loading project experiments:', error);
                experimentsList.innerHTML = '<tr class="placeholder-row"><td colspan="5" class="text-center text-danger">Error loading experiments</td></tr>';
            });
    }
    
    /**
     * Load project activity
     * 
     * @param {string} projectId - Project ID
     */
    function loadProjectActivity(projectId) {
        const activityList = document.getElementById('project-activity-list');
        activityList.innerHTML = '<li class="list-group-item text-center">Loading activity...</li>';
        
        API.get(`/projects/${projectId}/activity`)
            .then(activities => {
                if (!activities || activities.length === 0) {
                    activityList.innerHTML = '<li class="list-group-item text-center">No activity found</li>';
                    return;
                }
                
                let html = '';
                
                activities.forEach(activity => {
                    const date = new Date(activity.timestamp);
                    const timeAgo = getTimeAgo(date);
                    
                    html += `
                        <li class="list-group-item">
                            <div class="d-flex w-100 justify-content-between">
                                <h6 class="mb-1">${activity.title}</h6>
                                <small class="text-muted">${timeAgo}</small>
                            </div>
                            <p class="mb-1 small">${activity.description}</p>
                        </li>
                    `;
                });
                
                activityList.innerHTML = html;
            })
            .catch(error => {
                console.error('Error loading project activity:', error);
                activityList.innerHTML = '<li class="list-group-item text-center text-danger">Error loading activity</li>';
            });
    }
    
    /**
     * Get time ago string from date
     * 
     * @param {Date} date - Date object
     * @returns {string} Time ago string
     */
    function getTimeAgo(date) {
        const now = new Date();
        const diffMs = now - date;
        const diffSec = Math.floor(diffMs / 1000);
        const diffMin = Math.floor(diffSec / 60);
        const diffHour = Math.floor(diffMin / 60);
        const diffDay = Math.floor(diffHour / 24);
        
        if (diffSec < 60) {
            return 'just now';
        } else if (diffMin < 60) {
            return `${diffMin} minute${diffMin > 1 ? 's' : ''} ago`;
        } else if (diffHour < 24) {
            return `${diffHour} hour${diffHour > 1 ? 's' : ''} ago`;
        } else if (diffDay < 30) {
            return `${diffDay} day${diffDay > 1 ? 's' : ''} ago`;
        } else {
            return date.toLocaleDateString('en-US', {
                year: 'numeric',
                month: 'short',
                day: 'numeric'
            });
        }
    }
    
    /**
     * Show add experiment modal
     */
    function showAddExperimentModal() {
        if (!currentProjectId) {
            alert('No project selected');
            return;
        }
        
        // Reset form
        document.getElementById('add-experiment-form').reset();
        
        // Load available experiments
        loadAvailableExperiments();
        
        // Show modal
        const modal = new bootstrap.Modal(document.getElementById('add-experiment-modal'));
        modal.show();
    }
    
    /**
     * Load available experiments for adding to project
     */
    function loadAvailableExperiments() {
        const experimentSelect = document.getElementById('experiment-select');
        experimentSelect.innerHTML = '<option value="" selected disabled>Loading experiments...</option>';
        
        API.get('/experiments')
            .then(experiments => {
                if (!experiments || experiments.length === 0) {
                    experimentSelect.innerHTML = '<option value="" selected disabled>No experiments available</option>';
                    return;
                }
                
                let html = '<option value="" selected disabled>Choose an experiment...</option>';
                
                experiments.forEach(experiment => {
                    html += `<option value="${experiment.id}">${experiment.property_name} (${experiment.date_performed})</option>`;
                });
                
                experimentSelect.innerHTML = html;
            })
            .catch(error => {
                console.error('Error loading available experiments:', error);
                experimentSelect.innerHTML = '<option value="" selected disabled>Error loading experiments</option>';
            });
    }
    
    /**
     * Add experiment to project
     */
    function addExperimentToProject() {
        if (!currentProjectId) {
            alert('No project selected');
            return;
        }
        
        const experimentId = document.getElementById('experiment-select').value;
        const notes = document.getElementById('experiment-notes').value.trim();
        
        if (!experimentId) {
            alert('Please select an experiment');
            return;
        }
        
        const data = {
            experiment_id: experimentId,
            notes: notes
        };
        
        API.post(`/projects/${currentProjectId}/experiments`, data)
            .then(result => {
                // Hide modal
                bootstrap.Modal.getInstance(document.getElementById('add-experiment-modal')).hide();
                
                // Reload project experiments
                loadProjectExperiments(currentProjectId);
                
                // Show success message
                alert('Experiment added to project successfully');
            })
            .catch(error => {
                console.error('Error adding experiment to project:', error);
                alert('Error adding experiment to project: ' + error.message);
            });
    }
    
    /**
     * Redirect to create experiment page
     */
    function redirectToCreateExperiment() {
        window.location.href = '/experiments/new';
    }
    
    /**
     * Show import data modal
     */
    function showImportDataModal() {
        alert('Import data functionality not implemented yet');
    }
    
    /**
     * Show export results modal
     */
    function showExportResultsModal() {
        alert('Export results functionality not implemented yet');
    }
    
    /**
     * Show share project modal
     */
    function showShareProjectModal() {
        alert('Share project functionality not implemented yet');
    }
    
    /**
     * Apply saved search
     * 
     * @param {string} query - Search query
     */
    function applySavedSearch(query) {
        document.getElementById('project-search').value = query;
        searchProjects();
    }
    
    // Return public API
    return {
        init: init
    };
})();

// Initialize the project dashboard when the DOM is loaded
document.addEventListener('DOMContentLoaded', function() {
    ProjectDashboard.init();
});
