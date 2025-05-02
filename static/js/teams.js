/**
 * CryoProtect Analyzer - Teams JavaScript
 * 
 * This file contains the JavaScript code for the teams page, including
 * team management, member management, resource sharing, and notifications.
 */

document.addEventListener('DOMContentLoaded', function() {
    // Initialize API client
    const api = new ApiClient(supabaseUrl, supabaseKey);
    
    // DOM elements
    const teamsContainer = document.getElementById('teams-container');
    const createTeamButton = document.getElementById('create-team-button');
    const createTeamModal = document.getElementById('create-team-modal');
    const createTeamForm = document.getElementById('create-team-form');
    const addMemberModal = document.getElementById('add-member-modal');
    const addMemberForm = document.getElementById('add-member-form');
    const shareResourceModal = document.getElementById('share-resource-modal');
    const shareResourceForm = document.getElementById('share-resource-form');
    const teamDetailsModal = document.getElementById('team-details-modal');
    const notificationsButton = document.getElementById('notifications-button');
    const notificationsModal = document.getElementById('notifications-modal');
    const notificationsList = document.getElementById('notifications-list');
    const notificationBadge = document.getElementById('notification-badge');
    const markAllReadButton = document.getElementById('mark-all-read-button');
    
    // Current state
    let currentTeamId = null;
    let unreadNotificationsCount = 0;
    
    // Initialize
    loadTeams();
    loadNotifications();
    setInterval(loadNotifications, 60000); // Refresh notifications every minute
    
    // Event listeners
    createTeamButton.addEventListener('click', () => openModal(createTeamModal));
    createTeamForm.addEventListener('submit', handleCreateTeam);
    addMemberForm.addEventListener('submit', handleAddMember);
    shareResourceForm.addEventListener('submit', handleShareResource);
    notificationsButton.addEventListener('click', () => openModal(notificationsModal));
    markAllReadButton.addEventListener('click', handleMarkAllNotificationsRead);
    
    // Close modals when clicking on close button or cancel button
    document.querySelectorAll('.close-modal, .cancel-button').forEach(element => {
        element.addEventListener('click', () => {
            document.querySelectorAll('.modal').forEach(modal => {
                modal.style.display = 'none';
            });
        });
    });
    
    // Tab switching in team details modal
    document.querySelectorAll('.tab-button').forEach(button => {
        button.addEventListener('click', () => {
            const tabId = button.getAttribute('data-tab');
            
            // Update active tab button
            document.querySelectorAll('.tab-button').forEach(btn => {
                btn.classList.remove('active');
            });
            button.classList.add('active');
            
            // Show selected tab content
            document.querySelectorAll('.tab-content').forEach(content => {
                content.classList.remove('active');
            });
            document.getElementById(`${tabId}-tab`).classList.add('active');
            
            // Load tab-specific data
            if (tabId === 'members' && currentTeamId) {
                loadTeamMembers(currentTeamId);
            } else if (tabId === 'resources' && currentTeamId) {
                loadSharedResources(currentTeamId);
            } else if (tabId === 'activity' && currentTeamId) {
                loadTeamActivity(currentTeamId);
            }
        });
    });
    
    /**
     * Load all teams for the current user
     */
    function loadTeams() {
        teamsContainer.innerHTML = '<div class="loading">Loading teams...</div>';
        
        api.get('/api/v1/teams')
            .then(response => {
                if (response.length === 0) {
                    teamsContainer.innerHTML = '<p>You are not a member of any teams. Create a new team to get started.</p>';
                    return;
                }
                
                teamsContainer.innerHTML = '';
                response.forEach(team => {
                    const teamCard = createTeamCard(team);
                    teamsContainer.appendChild(teamCard);
                });
            })
            .catch(error => {
                console.error('Error loading teams:', error);
                teamsContainer.innerHTML = `<p>Error loading teams: ${error.message}</p>`;
            });
    }
    
    /**
     * Create a team card element
     */
    function createTeamCard(team) {
        const card = document.createElement('div');
        card.className = 'team-card';
        card.innerHTML = `
            <div class="team-header">
                <h3>${team.team_name}</h3>
                <div class="team-actions">
                    <button class="view-team-button secondary-button" data-team-id="${team.team_id}">View Details</button>
                </div>
            </div>
            <p>${team.team_description || 'No description'}</p>
            <div class="team-meta">
                <span>Your role: <strong>${team.role}</strong></span>
                <span>Members: ${team.member_count}</span>
            </div>
        `;
        
        // Add event listener for view details button
        card.querySelector('.view-team-button').addEventListener('click', () => {
            viewTeamDetails(team.team_id);
        });
        
        return card;
    }
    
    /**
     * View team details
     */
    function viewTeamDetails(teamId) {
        currentTeamId = teamId;
        
        api.get(`/api/v1/teams/${teamId}`)
            .then(team => {
                document.getElementById('team-details-name').textContent = team.name;
                document.getElementById('team-details-description').textContent = team.description || 'No description';
                
                // Set team ID for add member and share resource forms
                document.getElementById('team-id-member').value = teamId;
                document.getElementById('team-id-resource').value = teamId;
                
                // Load members tab by default
                loadTeamMembers(teamId);
                
                // Reset tabs
                document.querySelectorAll('.tab-button').forEach(btn => {
                    btn.classList.remove('active');
                });
                document.querySelector('.tab-button[data-tab="members"]').classList.add('active');
                
                document.querySelectorAll('.tab-content').forEach(content => {
                    content.classList.remove('active');
                });
                document.getElementById('members-tab').classList.add('active');
                
                // Add event listeners for add member and share resource buttons
                document.querySelector('.add-member-button').addEventListener('click', () => {
                    openModal(addMemberModal);
                });
                
                document.querySelector('.share-resource-button').addEventListener('click', () => {
                    openModal(shareResourceModal);
                });
                
                // Open the modal
                openModal(teamDetailsModal);
            })
            .catch(error => {
                console.error('Error loading team details:', error);
                alert(`Error loading team details: ${error.message}`);
            });
    }
    
    /**
     * Load team members
     */
    function loadTeamMembers(teamId) {
        const membersList = document.querySelector('.member-list');
        membersList.innerHTML = '<div class="loading">Loading members...</div>';
        
        api.get(`/api/v1/teams/${teamId}/members`)
            .then(members => {
                if (members.length === 0) {
                    membersList.innerHTML = '<p>No members found.</p>';
                    return;
                }
                
                membersList.innerHTML = '';
                members.forEach(member => {
                    const memberItem = document.createElement('div');
                    memberItem.className = 'member-item';
                    memberItem.innerHTML = `
                        <div class="member-info">
                            <div>${member.display_name || member.email}</div>
                            <div class="member-email">${member.email}</div>
                        </div>
                        <div class="member-role role-${member.role}">${member.role}</div>
                        <div class="member-actions">
                            <select class="role-select" data-user-id="${member.user_id}">
                                <option value="admin" ${member.role === 'admin' ? 'selected' : ''}>Admin</option>
                                <option value="editor" ${member.role === 'editor' ? 'selected' : ''}>Editor</option>
                                <option value="viewer" ${member.role === 'viewer' ? 'selected' : ''}>Viewer</option>
                            </select>
                            <button class="remove-member-button" data-user-id="${member.user_id}">Remove</button>
                        </div>
                    `;
                    
                    membersList.appendChild(memberItem);
                });
                
                // Add event listeners for role select and remove buttons
                document.querySelectorAll('.role-select').forEach(select => {
                    select.addEventListener('change', event => {
                        const userId = event.target.getAttribute('data-user-id');
                        const newRole = event.target.value;
                        updateMemberRole(teamId, userId, newRole);
                    });
                });
                
                document.querySelectorAll('.remove-member-button').forEach(button => {
                    button.addEventListener('click', event => {
                        const userId = event.target.getAttribute('data-user-id');
                        removeMember(teamId, userId);
                    });
                });
            })
            .catch(error => {
                console.error('Error loading team members:', error);
                membersList.innerHTML = `<p>Error loading members: ${error.message}</p>`;
            });
    }
    
    /**
     * Load shared resources
     */
    function loadSharedResources(teamId) {
        const resourcesList = document.querySelector('.resource-list');
        resourcesList.innerHTML = '<div class="loading">Loading shared resources...</div>';
        
        api.get(`/api/v1/teams/${teamId}/resources`)
            .then(resources => {
                if (resources.length === 0) {
                    resourcesList.innerHTML = '<p>No shared resources found.</p>';
                    return;
                }
                
                resourcesList.innerHTML = '';
                resources.forEach(resource => {
                    const resourceItem = document.createElement('div');
                    resourceItem.className = 'resource-item';
                    resourceItem.innerHTML = `
                        <div class="resource-info">
                            <div class="resource-type type-${resource.resource_type}">${resource.resource_type}</div>
                            <div>${resource.resource_id}</div>
                        </div>
                        <div class="resource-actions">
                            <button class="unshare-resource-button" data-resource-type="${resource.resource_type}" data-resource-id="${resource.resource_id}">Unshare</button>
                        </div>
                    `;
                    
                    resourcesList.appendChild(resourceItem);
                });
                
                // Add event listeners for unshare buttons
                document.querySelectorAll('.unshare-resource-button').forEach(button => {
                    button.addEventListener('click', event => {
                        const resourceType = event.target.getAttribute('data-resource-type');
                        const resourceId = event.target.getAttribute('data-resource-id');
                        unshareResource(teamId, resourceType, resourceId);
                    });
                });
            })
            .catch(error => {
                console.error('Error loading shared resources:', error);
                resourcesList.innerHTML = `<p>Error loading shared resources: ${error.message}</p>`;
            });
    }
    
    /**
     * Load team activity
     */
    function loadTeamActivity(teamId) {
        const activityFeed = document.querySelector('.activity-feed');
        activityFeed.innerHTML = '<div class="loading">Loading activity...</div>';
        
        api.get(`/api/v1/teams/${teamId}/activity`)
            .then(activities => {
                if (activities.length === 0) {
                    activityFeed.innerHTML = '<p>No activity found.</p>';
                    return;
                }
                
                activityFeed.innerHTML = '';
                activities.forEach(activity => {
                    const activityItem = document.createElement('div');
                    activityItem.className = 'activity-item';
                    
                    // Format the activity message based on the action and resource type
                    let message = '';
                    if (activity.action === 'created' && activity.resource_type === 'team') {
                        message = `Team was created`;
                    } else if (activity.action === 'added_member') {
                        message = `Added a new member to the team`;
                    } else if (activity.action === 'updated_member_role') {
                        message = `Updated a member's role`;
                    } else if (activity.action === 'removed_member') {
                        message = `Removed a member from the team`;
                    } else if (activity.action === 'shared_resource') {
                        message = `Shared a ${activity.resource_type} with the team`;
                    } else if (activity.action === 'unshared_resource') {
                        message = `Unshared a ${activity.resource_type} from the team`;
                    } else {
                        message = `${activity.action} ${activity.resource_type}`;
                    }
                    
                    // Format the timestamp
                    const timestamp = new Date(activity.created_at);
                    const formattedTime = timestamp.toLocaleString();
                    
                    activityItem.innerHTML = `
                        <div class="activity-user">${activity.user_email || 'Unknown user'}</div>
                        <div class="activity-message">${message}</div>
                        <div class="activity-time">${formattedTime}</div>
                    `;
                    
                    activityFeed.appendChild(activityItem);
                });
            })
            .catch(error => {
                console.error('Error loading team activity:', error);
                activityFeed.innerHTML = `<p>Error loading activity: ${error.message}</p>`;
            });
    }
    
    /**
     * Load notifications
     */
    function loadNotifications() {
        api.get('/api/v1/notifications?unread_only=true')
            .then(notifications => {
                unreadNotificationsCount = notifications.length;
                notificationBadge.textContent = unreadNotificationsCount;
                
                if (unreadNotificationsCount > 0) {
                    notificationBadge.style.display = 'inline';
                } else {
                    notificationBadge.style.display = 'none';
                }
            })
            .catch(error => {
                console.error('Error loading notifications:', error);
            });
    }
    
    /**
     * Display notifications in the modal
     */
    function displayNotifications() {
        notificationsList.innerHTML = '<div class="loading">Loading notifications...</div>';
        
        api.get('/api/v1/notifications')
            .then(notifications => {
                if (notifications.length === 0) {
                    notificationsList.innerHTML = '<p>No notifications found.</p>';
                    return;
                }
                
                notificationsList.innerHTML = '';
                notifications.forEach(notification => {
                    const notificationItem = document.createElement('div');
                    notificationItem.className = `notification-item ${notification.is_read ? 'read' : 'unread'}`;
                    
                    // Format the timestamp
                    const timestamp = new Date(notification.created_at);
                    const formattedTime = timestamp.toLocaleString();
                    
                    notificationItem.innerHTML = `
                        <div class="notification-header">
                            <div class="notification-title">${notification.title}</div>
                            <div class="notification-time">${formattedTime}</div>
                        </div>
                        <div class="notification-message">${notification.message}</div>
                        <div class="notification-actions">
                            ${!notification.is_read ? `<button class="mark-read-button" data-notification-id="${notification.id}">Mark as Read</button>` : ''}
                            <button class="delete-notification-button" data-notification-id="${notification.id}">Delete</button>
                        </div>
                    `;
                    
                    notificationsList.appendChild(notificationItem);
                });
                
                // Add event listeners for mark as read and delete buttons
                document.querySelectorAll('.mark-read-button').forEach(button => {
                    button.addEventListener('click', event => {
                        const notificationId = event.target.getAttribute('data-notification-id');
                        markNotificationAsRead(notificationId);
                    });
                });
                
                document.querySelectorAll('.delete-notification-button').forEach(button => {
                    button.addEventListener('click', event => {
                        const notificationId = event.target.getAttribute('data-notification-id');
                        deleteNotification(notificationId);
                    });
                });
            })
            .catch(error => {
                console.error('Error loading notifications:', error);
                notificationsList.innerHTML = `<p>Error loading notifications: ${error.message}</p>`;
            });
    }
    
    /**
     * Handle create team form submission
     */
    function handleCreateTeam(event) {
        event.preventDefault();
        
        const formData = new FormData(createTeamForm);
        const teamData = {
            name: formData.get('name'),
            description: formData.get('description')
        };
        
        api.post('/api/v1/teams', teamData)
            .then(response => {
                createTeamForm.reset();
                closeModal(createTeamModal);
                loadTeams();
            })
            .catch(error => {
                console.error('Error creating team:', error);
                alert(`Error creating team: ${error.message}`);
            });
    }
    
    /**
     * Handle add member form submission
     */
    function handleAddMember(event) {
        event.preventDefault();
        
        const formData = new FormData(addMemberForm);
        const teamId = formData.get('team_id');
        const email = formData.get('email');
        const role = formData.get('role');
        
        // First, we need to get the user ID from the email
        api.get(`/api/v1/users?email=${encodeURIComponent(email)}`)
            .then(users => {
                if (users.length === 0) {
                    throw new Error('User not found with that email');
                }
                
                const userId = users[0].id;
                
                // Now add the member to the team
                return api.post(`/api/v1/teams/${teamId}/members`, {
                    user_id: userId,
                    role: role
                });
            })
            .then(response => {
                addMemberForm.reset();
                closeModal(addMemberModal);
                loadTeamMembers(teamId);
            })
            .catch(error => {
                console.error('Error adding member:', error);
                alert(`Error adding member: ${error.message}`);
            });
    }
    
    /**
     * Handle share resource form submission
     */
    function handleShareResource(event) {
        event.preventDefault();
        
        const formData = new FormData(shareResourceForm);
        const teamId = formData.get('team_id');
        const resourceType = formData.get('resource_type');
        const resourceId = formData.get('resource_id');
        
        api.post(`/api/v1/teams/${teamId}/resources`, {
            resource_type: resourceType,
            resource_id: resourceId
        })
            .then(response => {
                shareResourceForm.reset();
                closeModal(shareResourceModal);
                loadSharedResources(teamId);
            })
            .catch(error => {
                console.error('Error sharing resource:', error);
                alert(`Error sharing resource: ${error.message}`);
            });
    }
    
    /**
     * Update a team member's role
     */
    function updateMemberRole(teamId, userId, newRole) {
        api.put(`/api/v1/teams/${teamId}/members/${userId}`, {
            role: newRole
        })
            .then(response => {
                loadTeamMembers(teamId);
            })
            .catch(error => {
                console.error('Error updating member role:', error);
                alert(`Error updating member role: ${error.message}`);
                loadTeamMembers(teamId); // Reload to reset the select
            });
    }
    
    /**
     * Remove a member from a team
     */
    function removeMember(teamId, userId) {
        if (!confirm('Are you sure you want to remove this member from the team?')) {
            return;
        }
        
        api.delete(`/api/v1/teams/${teamId}/members/${userId}`)
            .then(response => {
                loadTeamMembers(teamId);
            })
            .catch(error => {
                console.error('Error removing member:', error);
                alert(`Error removing member: ${error.message}`);
            });
    }
    
    /**
     * Unshare a resource from a team
     */
    function unshareResource(teamId, resourceType, resourceId) {
        if (!confirm('Are you sure you want to unshare this resource from the team?')) {
            return;
        }
        
        api.delete(`/api/v1/teams/${teamId}/resources/${resourceType}/${resourceId}`)
            .then(response => {
                loadSharedResources(teamId);
            })
            .catch(error => {
                console.error('Error unsharing resource:', error);
                alert(`Error unsharing resource: ${error.message}`);
            });
    }
    
    /**
     * Mark a notification as read
     */
    function markNotificationAsRead(notificationId) {
        api.put(`/api/v1/notifications/${notificationId}`)
            .then(response => {
                displayNotifications();
                loadNotifications();
            })
            .catch(error => {
                console.error('Error marking notification as read:', error);
                alert(`Error marking notification as read: ${error.message}`);
            });
    }
    
    /**
     * Delete a notification
     */
    function deleteNotification(notificationId) {
        api.delete(`/api/v1/notifications/${notificationId}`)
            .then(response => {
                displayNotifications();
                loadNotifications();
            })
            .catch(error => {
                console.error('Error deleting notification:', error);
                alert(`Error deleting notification: ${error.message}`);
            });
    }
    
    /**
     * Handle mark all notifications as read
     */
    function handleMarkAllNotificationsRead() {
        api.put('/api/v1/notifications')
            .then(response => {
                displayNotifications();
                loadNotifications();
            })
            .catch(error => {
                console.error('Error marking all notifications as read:', error);
                alert(`Error marking all notifications as read: ${error.message}`);
            });
    }
    
    /**
     * Open a modal
     */
    function openModal(modal) {
        // Close all other modals first
        document.querySelectorAll('.modal').forEach(m => {
            m.style.display = 'none';
        });
        
        // Open the requested modal
        modal.style.display = 'block';
        
        // If it's the notifications modal, load the notifications
        if (modal === notificationsModal) {
            displayNotifications();
        }
    }
    
    /**
     * Close a modal
     */
    function closeModal(modal) {
        modal.style.display = 'none';
    }
});