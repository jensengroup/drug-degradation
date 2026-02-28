/**
 * Interactive Pathway Visualization Module (Phase 4)
 *
 * This module provides interactive reaction profile viewing with transition state
 * markers and energy barrier annotations, matching the style from paths.ipynb
 */

class PathwayProfileViewer {
    constructor(networkData) {
        this.networkData = networkData;
        this.pathwayAnalyzer = window.pathAnalyzer;
        this.currentProfile = null;
        this.modalId = 'pathwayProfileModal';

        this.initializeViewer();
    }

    /**
     * Initialize the pathway profile viewer UI as a modal
     */
    initializeViewer() {
        // Create modal HTML if it doesn't exist
        if (!document.getElementById(this.modalId)) {
            const modalHTML = `
                <div id="${this.modalId}" class="modal pathway-modal" style="display: none;">
                    <div class="modal-content pathway-modal-content">
                        <div class="modal-header">
                            <h3>Pathway Profile Viewer</h3>
                            <span class="close" onclick="hidePathwayViewer()">&times;</span>
                        </div>
                        <div class="pathway-controls">
                            <div class="control-row">
                                <label for="node-selector">Target Node:</label>
                                <select id="node-selector">
                                    <option value="">Select a node...</option>
                                </select>
                                <label for="pathway-selector">Pathway:</label>
                                <select id="pathway-selector">
                                    <option value="">Select pathway...</option>
                                </select>
                            </div>
                            <div class="control-row">
                                <button id="show-profile-btn" disabled>Show Profile</button>
                                <button id="compare-pathways-btn" disabled>Compare Pathways</button>
                                <button id="export-profile-btn" disabled>Export Data</button>
                            </div>
                        </div>
                        <div class="pathway-display">
                            <div id="profile-plot" class="profile-plot"></div>
                            <div id="profile-info" class="profile-info"></div>
                        </div>
                    </div>
                </div>
            `;
            document.body.insertAdjacentHTML('beforeend', modalHTML);
        }

        this.modal = document.getElementById(this.modalId);
        this.setupEventHandlers();
        this.addStyles();
        this.populateNodeSelector();
    }

    /**
     * Add CSS styles for the pathway viewer
     */
    addStyles() {
        if (document.getElementById('pathway-viewer-styles')) {
            return; // Styles already added
        }

        const style = document.createElement('style');
        style.id = 'pathway-viewer-styles';
        style.textContent = `
            /* Modal Styles */
            .pathway-modal {
                position: fixed;
                z-index: 1000;
                left: 0;
                top: 0;
                width: 100%;
                height: 100%;
                background-color: rgba(0,0,0,0.5);
            }

            .pathway-modal-content {
                background-color: #fefefe;
                margin: 2% auto;
                padding: 20px;
                border: 1px solid #888;
                border-radius: 8px;
                width: 90%;
                max-width: 1200px;
                max-height: 90vh;
                overflow-y: auto;
            }

            .modal-header {
                display: flex;
                justify-content: space-between;
                align-items: center;
                border-bottom: 1px solid #ddd;
                padding-bottom: 10px;
                margin-bottom: 20px;
            }

            .modal-header h3 {
                margin: 0;
                color: #333;
            }

            .close {
                color: #aaa;
                font-size: 28px;
                font-weight: bold;
                cursor: pointer;
            }

            .close:hover,
            .close:focus {
                color: black;
                text-decoration: none;
            }

            .pathway-viewer {
                border: 1px solid #ddd;
                border-radius: 8px;
                padding: 15px;
                margin: 10px 0;
                background: #fafafa;
            }

            .pathway-controls h3 {
                margin: 0 0 15px 0;
                color: #333;
            }

            .control-row {
                display: flex;
                align-items: center;
                gap: 10px;
                margin-bottom: 10px;
            }

            .control-row label {
                font-weight: bold;
                min-width: 80px;
            }

            .control-row select {
                padding: 5px;
                border: 1px solid #ccc;
                border-radius: 4px;
                min-width: 150px;
            }

            .control-row button {
                padding: 8px 15px;
                border: none;
                border-radius: 4px;
                cursor: pointer;
                font-weight: bold;
            }

            .control-row button:enabled {
                background: #007bff;
                color: white;
            }

            .control-row button:disabled {
                background: #ccc;
                color: #666;
                cursor: not-allowed;
            }

            .pathway-display {
                margin-top: 20px;
                border-top: 1px solid #ddd;
                padding-top: 20px;
            }

            .profile-plot {
                min-height: 400px;
                border: 1px solid #ddd;
                border-radius: 4px;
                background: white;
                padding: 10px;
                margin-bottom: 15px;
                position: relative;
            }

            .profile-info {
                background: white;
                border: 1px solid #ddd;
                border-radius: 4px;
                padding: 15px;
            }

            .profile-point {
                position: absolute;
                width: 8px;
                height: 8px;
                border-radius: 50%;
                background: #007bff;
                cursor: pointer;
                transform: translate(-50%, -50%);
            }

            .profile-point.transition-state {
                background: #dc3545;
                width: 10px;
                height: 10px;
            }

            .profile-point:hover {
                transform: translate(-50%, -50%) scale(1.2);
            }

            .profile-line {
                position: absolute;
                height: 2px;
                background: #007bff;
                transform: translateY(-50%);
            }

            .profile-label {
                position: absolute;
                font-size: 12px;
                font-weight: bold;
                background: rgba(255, 255, 255, 0.9);
                padding: 2px 6px;
                border-radius: 3px;
                border: 1px solid #ddd;
                transform: translate(-50%, -100%);
                margin-top: -5px;
            }

            .barrier-annotation {
                position: absolute;
                border-left: 2px dashed #ff6b35;
                color: #ff6b35;
                font-size: 11px;
                font-weight: bold;
            }

            .axis {
                position: absolute;
                color: #666;
                font-size: 12px;
            }

            .axis-x {
                bottom: 5px;
                left: 50%;
                transform: translateX(-50%);
            }

            .axis-y {
                left: 5px;
                top: 50%;
                transform: translateY(-50%) rotate(-90deg);
                transform-origin: center;
            }
        `;

        document.head.appendChild(style);
    }

    /**
     * Set up event handlers for the viewer controls
     */
    setupEventHandlers() {
        const nodeSelector = document.getElementById('node-selector');
        const pathwaySelector = document.getElementById('pathway-selector');
        const showProfileBtn = document.getElementById('show-profile-btn');
        const comparePathwaysBtn = document.getElementById('compare-pathways-btn');
        const exportProfileBtn = document.getElementById('export-profile-btn');

        nodeSelector.addEventListener('change', () => {
            this.onNodeSelectionChange();
        });

        pathwaySelector.addEventListener('change', () => {
            this.onPathwaySelectionChange();
        });

        showProfileBtn.addEventListener('click', () => {
            this.showProfile();
            this.updateMainReactionChart(); // Update main chart when showing profile
        });

        comparePathwaysBtn.addEventListener('click', () => {
            this.comparePathways();
        });

        exportProfileBtn.addEventListener('click', () => {
            this.exportProfileData();
        });

        // Add event listener for pathway selection change to update main chart
        pathwaySelector.addEventListener('change', () => {
            this.onPathwaySelectionChange();
            if (pathwaySelector.value) {
                this.updateMainReactionChart(); // Update main chart when pathway changes
            }
        });
    }

    /**
     * Initialize with network data and populate node selector
     * @param {Object} networkData - Network data object
     */
    initialize(networkData) {
        this.networkData = networkData;
        this.populateNodeSelector();
    }

    /**
     * Populate the node selector with available nodes
     */
    populateNodeSelector() {
        const nodeSelector = document.getElementById('node-selector');
        nodeSelector.innerHTML = '<option value="">Select a node...</option>';

        if (!this.networkData || !this.networkData.precomputedPathways) {
            console.warn('No pathway data available in networkData');
            return;
        }

        // Add nodes that have pathways
        const nodeIds = Object.keys(this.networkData.precomputedPathways)
            .map(id => parseInt(id))
            .sort((a, b) => a - b);

        console.log(`Found ${nodeIds.length} nodes with pathways`);

        for (const nodeId of nodeIds) {
            const pathways = this.networkData.precomputedPathways[nodeId];
            if (pathways && pathways.length > 0) {
                const option = document.createElement('option');
                option.value = nodeId;
                option.textContent = `Node ${nodeId} (${pathways.length} pathways)`;
                nodeSelector.appendChild(option);
            }
        }
    }

    /**
     * Handle node selection change
     */
    onNodeSelectionChange() {
        const nodeSelector = document.getElementById('node-selector');
        const pathwaySelector = document.getElementById('pathway-selector');
        const showProfileBtn = document.getElementById('show-profile-btn');

        pathwaySelector.innerHTML = '<option value="">Select pathway...</option>';
        showProfileBtn.disabled = true;

        const selectedNodeId = nodeSelector.value;
        if (!selectedNodeId) {
            return;
        }

        // Populate pathway selector for this node
        const pathways = this.networkData.precomputedPathways[parseInt(selectedNodeId)];
        if (pathways && pathways.length > 0) {
            pathways.forEach((pathway, index) => {
                const option = document.createElement('option');
                option.value = index;
                option.textContent = `Path ${index + 1} (${pathway.path.length} steps, barrier: ${pathway.barrier_height.toFixed(2)} kcal/mol)`;
                pathwaySelector.appendChild(option);
            });
        }
    }

    /**
     * Handle pathway selection change
     */
    onPathwaySelectionChange() {
        const pathwaySelector = document.getElementById('pathway-selector');
        const showProfileBtn = document.getElementById('show-profile-btn');

        showProfileBtn.disabled = !pathwaySelector.value;
    }

    /**
     * Show the selected pathway profile
     */
    showProfile() {
        const nodeSelector = document.getElementById('node-selector');
        const pathwaySelector = document.getElementById('pathway-selector');

        const nodeId = parseInt(nodeSelector.value);
        const pathwayIndex = parseInt(pathwaySelector.value);

        if (!nodeId || pathwayIndex === undefined) {
            return;
        }

        // Get pathway data from networkData
        const pathways = this.networkData.precomputedPathways[nodeId];
        if (!pathways || !pathways[pathwayIndex]) {
            this.showError('No pathway data available for selected node/pathway');
            return;
        }

        const pathway = pathways[pathwayIndex];
        const profile = pathway.profile;

        if (!profile || profile.length === 0) {
            this.showError('No profile data available for selected pathway');
            return;
        }

        this.currentProfile = {
            nodeId,
            pathwayIndex,
            profile,
            pathway
        };

        this.renderProfile(profile);
        this.updateProfileInfo(nodeId, pathwayIndex, profile);
        this.enableExportButton();
    }

    /**
     * Render the reaction profile plot
     * @param {Array} profile - Reaction profile data
     */
    renderProfile(profile) {
        const plotContainer = document.getElementById('profile-plot');
        plotContainer.innerHTML = '';

        if (!profile || profile.length === 0) {
            plotContainer.innerHTML = '<p>No profile data to display</p>';
            return;
        }

        // Calculate plot dimensions and scaling
        const plotWidth = plotContainer.clientWidth - 40;
        const plotHeight = plotContainer.clientHeight - 40;
        const margin = 20;

        // Find energy range
        const energies = profile.map(p => p.energy);
        const minEnergy = Math.min(...energies);
        const maxEnergy = Math.max(...energies);
        const energyRange = Math.max(maxEnergy - minEnergy, 10); // Minimum range of 10 kcal/mol

        // Find coordinate range
        const coords = profile.map(p => p.reactionCoordinate || p.step);
        const minCoord = Math.min(...coords);
        const maxCoord = Math.max(...coords);
        const coordRange = maxCoord - minCoord || 1;

        // Create axes
        plotContainer.innerHTML = `
            <div class="axis axis-x">Reaction Coordinate</div>
            <div class="axis axis-y">Energy (kcal/mol)</div>
        `;

        // Draw profile points and lines
        for (let i = 0; i < profile.length; i++) {
            const point = profile[i];
            const x = margin + ((point.reactionCoordinate || point.step) - minCoord) / coordRange * (plotWidth - 2 * margin);
            const y = margin + (maxEnergy - point.energy + energyRange * 0.1) / (energyRange * 1.2) * (plotHeight - 2 * margin);

            // Create point element
            const pointEl = document.createElement('div');
            pointEl.className = `profile-point ${point.type === 'transition_state' ? 'transition-state' : ''}`;
            pointEl.style.left = `${x}px`;
            pointEl.style.top = `${y}px`;
            pointEl.title = `${point.nodeId}: ${point.energy.toFixed(1)} kcal/mol`;

            plotContainer.appendChild(pointEl);

            // Create label
            const labelEl = document.createElement('div');
            labelEl.className = 'profile-label';
            labelEl.style.left = `${x}px`;
            labelEl.style.top = `${y}px`;
            labelEl.textContent = point.nodeId;

            plotContainer.appendChild(labelEl);

            // Draw line to next point
            if (i < profile.length - 1) {
                const nextPoint = profile[i + 1];
                const nextX = margin + ((nextPoint.reactionCoordinate || nextPoint.step) - minCoord) / coordRange * (plotWidth - 2 * margin);
                const nextY = margin + (maxEnergy - nextPoint.energy + energyRange * 0.1) / (energyRange * 1.2) * (plotHeight - 2 * margin);

                const lineEl = document.createElement('div');
                lineEl.className = 'profile-line';
                lineEl.style.left = `${x}px`;
                lineEl.style.top = `${y}px`;
                lineEl.style.width = `${Math.sqrt((nextX - x) ** 2 + (nextY - y) ** 2)}px`;
                lineEl.style.transformOrigin = '0 50%';
                lineEl.style.transform = `translateY(-50%) rotate(${Math.atan2(nextY - y, nextX - x)}rad)`;

                plotContainer.appendChild(lineEl);
            }
        }

        // Add barrier annotations
        this.addBarrierAnnotations(plotContainer, profile, plotWidth, plotHeight, margin, minEnergy, maxEnergy, energyRange, minCoord, coordRange);
    }

    /**
     * Add barrier height annotations to the plot
     */
    addBarrierAnnotations(container, profile, plotWidth, plotHeight, margin, minEnergy, maxEnergy, energyRange, minCoord, coordRange) {
        // Find the maximum barrier
        let maxBarrierHeight = 0;
        let maxBarrierIndex = -1;

        let runningMin = profile[0].energy;
        for (let i = 0; i < profile.length; i++) {
            const currentEnergy = profile[i].energy;
            const barrier = currentEnergy - runningMin;

            if (barrier > maxBarrierHeight) {
                maxBarrierHeight = barrier;
                maxBarrierIndex = i;
            }

            runningMin = Math.min(runningMin, currentEnergy);
        }

        // Add annotation for maximum barrier
        if (maxBarrierHeight > 0.1 && maxBarrierIndex >= 0) {
            const point = profile[maxBarrierIndex];
            const x = margin + ((point.reactionCoordinate || point.step) - minCoord) / coordRange * (plotWidth - 2 * margin);
            const y = margin + (maxEnergy - point.energy + energyRange * 0.1) / (energyRange * 1.2) * (plotHeight - 2 * margin);

            const annotation = document.createElement('div');
            annotation.className = 'barrier-annotation';
            annotation.style.left = `${x + 10}px`;
            annotation.style.top = `${y - 30}px`;
            annotation.style.height = '25px';
            annotation.innerHTML = `Max Barrier<br>${maxBarrierHeight.toFixed(1)} kcal/mol`;

            container.appendChild(annotation);
        }
    }

    /**
     * Update the profile information panel
     */
    updateProfileInfo(nodeId, pathwayIndex, profile) {
        const infoContainer = document.getElementById('profile-info');

        const pathway = this.pathwayAnalyzer.getAllPathways(nodeId)[pathwayIndex];
        const barrier = this.pathwayAnalyzer.getPathwayBarrier(nodeId, pathwayIndex);
        const length = this.pathwayAnalyzer.getPathwayLength(nodeId, pathwayIndex);

        // Count transition states
        const transitionStates = profile.filter(p => p.type === 'transition_state');

        infoContainer.innerHTML = `
            <h4>Pathway Information</h4>
            <div class="info-grid">
                <div><strong>Target Node:</strong> ${nodeId}</div>
                <div><strong>Pathway:</strong> ${pathwayIndex + 1}</div>
                <div><strong>Path Length:</strong> ${length} steps</div>
                <div><strong>Max Barrier:</strong> ${barrier.toFixed(1)} kcal/mol</div>
                <div><strong>Transition States:</strong> ${transitionStates.length}</div>
                <div><strong>Path:</strong> ${pathway.path.join(' → ')}</div>
            </div>

            ${transitionStates.length > 0 ? `
                <h5>Transition States</h5>
                <ul>
                    ${transitionStates.map(ts => `
                        <li>${ts.nodeId}: ${ts.energy.toFixed(1)} kcal/mol</li>
                    `).join('')}
                </ul>
            ` : ''}
        `;
    }

    /**
     * Enable the export button
     */
    enableExportButton() {
        const exportBtn = document.getElementById('export-profile-btn');
        exportBtn.disabled = false;
    }

    /**
     * Show error message in the plot area
     */
    showError(message) {
        const plotContainer = document.getElementById('profile-plot');
        plotContainer.innerHTML = `<p style="color: red; text-align: center; padding: 20px;">${message}</p>`;
    }

    /**
     * Compare multiple pathways (Phase 4.2)
     */
    comparePathways() {
        // TODO: Implement pathway comparison
        alert('Pathway comparison feature coming soon!');
    }

    /**
     * Export profile data to CSV (Phase 4.3)
     */
    exportProfileData() {
        if (!this.currentProfile) {
            alert('No profile data to export');
            return;
        }

        const { nodeId, pathwayIndex, profile } = this.currentProfile;

        // Create CSV content
        const headers = ['Step', 'NodeId', 'Energy_kcal_mol', 'ReactionCoordinate', 'Type'];
        const rows = profile.map(point => [
            point.step,
            point.nodeId,
            point.energy.toFixed(3),
            (point.reactionCoordinate || point.step).toFixed(3),
            point.type || 'minimum'
        ]);

        const csvContent = [headers, ...rows]
            .map(row => row.join(','))
            .join('\\n');

        // Download the file
        const blob = new Blob([csvContent], { type: 'text/csv' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `pathway_profile_node${nodeId}_path${pathwayIndex + 1}.csv`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    }

    /**
     * Show pathway profile for a specific node (called from network interaction)
     * @param {string} nodeId - The node ID to show pathways for
     */
    showPathwayProfile(nodeId) {
        this.modal.style.display = 'block';

        // Select the node and show first pathway
        const nodeSelector = document.getElementById('node-selector');
        nodeSelector.value = nodeId;
        this.onNodeSelectionChange();

        // If pathways available, select first one and show profile
        const pathwaySelector = document.getElementById('pathway-selector');
        if (pathwaySelector.options.length > 1) {
            pathwaySelector.selectedIndex = 1; // First pathway (index 0 is "Select pathway...")
            this.onPathwaySelectionChange();
            this.showProfile();
            this.updateMainReactionChart(); // Update main chart when node is selected
        }
    }

    /**
     * Update the main reaction profile chart with precomputed pathway data
     * This integrates pathway selection with the main interface's reaction chart
     */
    updateMainReactionChart() {
        // Get currently selected node and pathway
        const nodeSelector = document.getElementById('node-selector');
        const pathwaySelector = document.getElementById('pathway-selector');

        if (!nodeSelector.value || !pathwaySelector.value) {
            return;
        }

        const nodeId = parseInt(nodeSelector.value);
        const pathwayIndex = parseInt(pathwaySelector.value);

        // Get the pathway data
        const pathways = this.networkData.precomputedPathways[nodeId];
        if (!pathways || !pathways[pathwayIndex]) {
            console.warn(`No pathway data found for node ${nodeId}, pathway ${pathwayIndex}`);
            return;
        }

        const pathway = pathways[pathwayIndex];
        const profile = pathway.profile;

        if (!profile || profile.length === 0) {
            console.warn(`No profile data found for node ${nodeId}, pathway ${pathwayIndex}`);
            return;
        }

        // Access the main chart through the network app UI
        if (!window.networkApp || !window.networkApp.ui || !window.networkApp.ui.reactionChart) {
            console.warn('Main reaction chart not available');
            return;
        }

        const chart = window.networkApp.ui.reactionChart;

        // Prepare chart data from the precomputed profile
        const labels = [];
        const energyData = [];
        const pointColors = [];
        const pointRadiuses = [];

        // Process each point in the reaction profile
        for (let i = 0; i < profile.length; i++) {
            const point = profile[i];

            // Use node ID as label, or step number if no node ID
            labels.push(point.node_id ? `Node ${point.node_id}` : `Step ${i + 1}`);

            // Add energy value (already in kcal/mol from precomputed data)
            energyData.push(point.energy);

            // Color transition states differently
            if (point.type === 'transition_state') {
                pointColors.push('#ff6b6b'); // Red for transition states
                pointRadiuses.push(8);       // Larger radius for transition states
            } else {
                pointColors.push('#4dabf7'); // Blue for regular states
                pointRadiuses.push(6);       // Normal radius
            }
        }

        // Update chart data
        chart.data.labels = labels;
        chart.data.datasets[0].data = energyData;
        chart.data.datasets[0].pointBackgroundColor = pointColors;
        chart.data.datasets[0].pointRadius = pointRadiuses;

        // Update chart title and labels
        chart.data.datasets[0].label = `Pathway to Node ${nodeId} (Barrier: ${pathway.barrier_height.toFixed(1)} kcal/mol)`;

        // Update the chart
        chart.update('active');

        // Show the reaction profile chart if it's hidden
        const chartDiv = document.getElementById('reactionProfileChart');
        if (chartDiv && chartDiv.style.display === 'none') {
            chartDiv.style.display = 'block';
        }

        // Add visual indicator that pathway is active
        if (chartDiv) {
            chartDiv.classList.add('pathway-active');

            // Remove the class after a few seconds to return to normal state
            setTimeout(() => {
                chartDiv.classList.remove('pathway-active');
            }, 3000);
        }

        console.log(`✅ Updated main reaction chart with pathway data for node ${nodeId}, pathway ${pathwayIndex}`);
    }

    /**
     * Show reaction profile for a selected path (sequence of nodes)
     * This method is called when a path is selected in the network visualization
     * @param {Array} pathNodes - Array of node IDs representing the selected path
     */
    showReactionProfileForPath(pathNodes) {
        if (!pathNodes || pathNodes.length < 2) {
            console.warn('Invalid path: need at least 2 nodes');
            return;
        }

        console.log(`🎯 Showing reaction profile for path: [${pathNodes.join(' → ')}]`);

        // Find which target node this path leads to
        const targetNode = pathNodes[pathNodes.length - 1];

        // Look for a matching pathway to this target node
        const pathways = this.networkData.precomputedPathways[targetNode];
        if (!pathways || pathways.length === 0) {
            console.warn(`No pathways found for target node ${targetNode}`);
            return;
        }

        // Find the pathway that matches the selected path
        let matchingPathwayIndex = -1;
        for (let i = 0; i < pathways.length; i++) {
            const pathway = pathways[i];
            if (pathway.path && this.pathsMatch(pathway.path, pathNodes)) {
                matchingPathwayIndex = i;
                break;
            }
        }

        if (matchingPathwayIndex === -1) {
            // If no exact match, use the first (best) pathway for the target node
            console.log(`No exact pathway match found, using best pathway for node ${targetNode}`);
            matchingPathwayIndex = 0;
        }

        // Update the main reaction chart with this pathway
        this.showReactionProfileForNodeAndPathway(targetNode, matchingPathwayIndex);
    }

    /**
     * Helper method to check if two paths match
     * @param {Array} path1 - First path array
     * @param {Array} path2 - Second path array
     * @returns {boolean} True if paths match
     */
    pathsMatch(path1, path2) {
        if (path1.length !== path2.length) return false;
        for (let i = 0; i < path1.length; i++) {
            if (path1[i] !== path2[i]) return false;
        }
        return true;
    }

    /**
     * Show reaction profile for a specific node and pathway index
     * @param {number} nodeId - Target node ID
     * @param {number} pathwayIndex - Index of the pathway
     */
    showReactionProfileForNodeAndPathway(nodeId, pathwayIndex) {
        const pathways = this.networkData.precomputedPathways[nodeId];
        if (!pathways || !pathways[pathwayIndex]) {
            console.warn(`No pathway data found for node ${nodeId}, pathway ${pathwayIndex}`);
            return;
        }

        const pathway = pathways[pathwayIndex];
        const profile = pathway.profile;

        if (!profile || profile.length === 0) {
            console.warn(`No profile data found for node ${nodeId}, pathway ${pathwayIndex}`);
            return;
        }

        // Access the main chart through the network app UI
        if (!window.networkApp || !window.networkApp.ui || !window.networkApp.ui.reactionChart) {
            console.warn('Main reaction chart not available');
            return;
        }

        const chart = window.networkApp.ui.reactionChart;

        // Prepare chart data from the precomputed profile
        const labels = [];
        const energyData = [];
        const pointColors = [];
        const pointRadiuses = [];

        // Process each point in the reaction profile
        for (let i = 0; i < profile.length; i++) {
            const point = profile[i];

            // Use node ID as label, or step number if no node ID
            labels.push(point.node_id ? `Node ${point.node_id}` : `Step ${i + 1}`);

            // Add energy value (already in kcal/mol from precomputed data)
            energyData.push(point.energy);

            // Color transition states differently
            if (point.type === 'transition_state') {
                pointColors.push('#ff6b6b'); // Red for transition states
                pointRadiuses.push(8);       // Larger radius for transition states
            } else {
                pointColors.push('#4dabf7'); // Blue for regular states
                pointRadiuses.push(6);       // Normal radius
            }
        }

        // Update chart data
        chart.data.labels = labels;
        chart.data.datasets[0].data = energyData;
        chart.data.datasets[0].pointBackgroundColor = pointColors;
        chart.data.datasets[0].pointRadius = pointRadiuses;

        // Update chart title and labels
        chart.data.datasets[0].label = `Pathway to Node ${nodeId} (${pathway.barrier_height.toFixed(1)} kcal/mol barrier)`;

        // Update the chart
        chart.update('active');

        // Show the reaction profile chart if it's hidden
        const chartDiv = document.getElementById('reactionProfileChart');
        if (chartDiv && chartDiv.style.display === 'none') {
            chartDiv.style.display = 'block';
        }

        // Add visual indicator that pathway is active
        if (chartDiv) {
            chartDiv.classList.add('pathway-active');

            // Remove the class after a few seconds to return to normal state
            setTimeout(() => {
                chartDiv.classList.remove('pathway-active');
            }, 3000);
        }

        console.log(`✅ Updated main reaction chart with pathway to node ${nodeId} (pathway ${pathwayIndex})`);
    }

    /**
     * Hide the pathway viewer modal
     */
    hideViewer() {
        this.modal.style.display = 'none';
    }
}

// Export for use in other modules
window.PathwayProfileViewer = PathwayProfileViewer;
