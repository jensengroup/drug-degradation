/**
 * Main network application entry point
 */
class NetworkApp {
    constructor() {
        this.networkData = null;
        this.filters = null;
        this.ui = null;
        this.stats = null;
        this.network = null;
        this.originalNodeData = [];
        this.originalNetworkEdgeData = [];
        this.edgeData = [];
        this.graphAdjacency = new Map();
    }

    /**
     * Initialize the application with data from server
     */
    initialize(initialData) {
        // Initialize network data
        this.networkData = new NetworkData(initialData);

        // Create the vis.js network
        this.createNetwork();

        // Initial setup and event listeners (will be setup after network is ready)
        this.setupNetworkEvents();
        this.setupControlEvents();
    }

    /**
     * Create the vis.js network
     */
    createNetwork() {
        // Create initial nodes and edges data
        const nodes = new vis.DataSet();
        const edges = new vis.DataSet();

        // Populate with initial molecular data
        this.populateInitialData(nodes, edges);

        // Create network options
        const options = {
            nodes: {
                shape: 'circle',
                size: 16,
                font: {
                    size: 12,
                    color: '#000000'
                },
                borderWidth: 2
            },
            edges: {
                width: 2,
                color: { inherit: 'from' },
                arrows: { to: { enabled: true, scaleFactor: 1 } }
            },
            physics: {
                enabled: true,
                stabilization: { iterations: 100 }
            },
            interaction: {
                multiselect: true,
                selectConnectedEdges: false
            },
            configure: {
                enabled: false
            }
        };

        // Create network
        const container = document.getElementById('mynetwork');
        this.network = new vis.Network(container, { nodes, edges }, options);

        // Make network globally available for backward compatibility
        window.network = this.network;
    }

    /**
     * Populate network with initial molecular data
     */
    populateInitialData(nodes, edges) {
        // Add nodes from molecular data
        const nodeData = [];
        Object.keys(this.networkData.molEnergyMap).forEach(nodeId => {
            const id = parseInt(nodeId);
            const isSubstrate = id === this.networkData.substrateId;

            nodeData.push({
                id: id,
                label: isSubstrate ? 'Substrate' : id.toString(),
                color: isSubstrate ? '#4E4E4EFF' : '#C2C2C2',
                shape: isSubstrate ? 'ellipse' : 'circle',
                font: isSubstrate ? { color: 'white' } : {},
                title: `Node ${id}\n∆${this.networkData.energyType}: ${(this.networkData.relEnergyMap[id] || 0).toFixed(2)} kcal/mol`
            });
        });
        nodes.add(nodeData);

        // Add edges from reaction data (filter out self-edges)
        const edgeData = this.networkData.originalEdgeData
            .filter(edge => edge.begin !== edge.end) // Remove self-edges
            .map((edge, index) => {
                return {
                    id: index,
                    from: edge.begin,
                    to: edge.end,
                    title: this.getEdgeTitle(edge),
                    color: this.getEdgeColor(edge.type),
                    width: this.getEdgeWidth(edge),
                    // Store the original edge data for later access
                    originalData: edge
                };
            });
        edges.add(edgeData);

        // Also filter originalEdgeData to remove self-edges for consistency
        this.networkData.originalEdgeData = this.networkData.originalEdgeData.filter(edge => edge.begin !== edge.end);


    }

    /**
     * Get edge title for tooltip
     */
    getEdgeTitle(edge) {
        // Use originalData if available (when updating tooltips later)
        const edgeData = edge.originalData || edge;

        let title = edgeData.type.charAt(0).toUpperCase() + edgeData.type.slice(1);
        if (edgeData.type.toLowerCase() === 'reaction' && edgeData.count) {
            title += ` (${edgeData.count}x)`;
        }
        if (edgeData.smaller_products && edgeData.smaller_products.length > 0) {
            title += `<br>Eliminated Products: ${edgeData.smaller_products.join(', ')}`;
        }
        if (edgeData.rxn_energy !== undefined && edgeData.rxn_energy !== null) {
            title += `<br>∆${this.networkData.energyType} ${edgeData.rxn_energy.toFixed(2)} kcal/mol`;
        }
        if (edgeData.barrier !== undefined && edgeData.barrier !== null) {
            title += `<br>Barrier: ${edgeData.barrier.toFixed(2)} kcal/mol`;
        }

        return title;
    }

    /**
     * Get edge color based on reaction type
     */
    getEdgeColor(rxnType) {
        const type = rxnType.toLowerCase();
        if (type === 'tautomerization' || type === 'tautomer') return '#55A69A';
        if (type === 'protonation') return '#E2A334';
        if (type === 'deprotonation') return '#3C56E9';
        if (type === 'reaction') return '#808080';
        return '#00000084';
    }

    /**
     * Get edge width based on reaction count
     */
    getEdgeWidth(edge) {
        if (edge.type.toLowerCase() === 'reaction' && edge.count) {
            return Math.min(Math.max(1, edge.count), 20) / 5;
        }
        return 1;
    }

    /**
     * Setup network event handlers
     */
    setupNetworkEvents() {
        // Wait for network to be ready
        this.network.on("afterDrawing", () => {
            if (this.originalNetworkEdgeData.length === 0) {
                this.captureOriginalData();
                this.setupComponents();
                this.initializeUI();
            }
        });

        // Node selection events
        this.network.on("selectNode", (params) => {
            if (params.nodes.length > 0) {
                const nodeId = params.nodes[0];
                const ctrlKey = params.event.srcEvent && (params.event.srcEvent.ctrlKey || params.event.srcEvent.metaKey);
                const shiftKey = params.event.srcEvent && params.event.srcEvent.shiftKey;

                if (shiftKey && this.ui.lastSelectedNode !== null) {
                    this.ui.toggleNodeSelection(nodeId, false, true);
                } else if (ctrlKey) {
                    this.ui.toggleNodeSelection(nodeId, true, false);
                } else {
                    this.ui.toggleNodeSelection(nodeId, false, false);
                }
            }
        });

        // Double-click to select pathway from substrate to target node
        // Hold Shift+Double-click to open pathway profile viewer instead
        this.network.on("doubleClick", (params) => {
            if (params.nodes.length > 0) {
                const nodeId = params.nodes[0];
                const shiftKey = params.event.srcEvent && params.event.srcEvent.shiftKey;

                if (shiftKey) {
                    // Shift+Double-click: Open pathway profile viewer
                    this.showPathwayProfile(nodeId);
                } else {
                    // Regular double-click: Select pathway nodes
                    this.selectPathwayToNode(nodeId);
                }
            }
        });

        this.network.on("deselectNode", (params) => {
            if (params.previousSelection && params.previousSelection.nodes.length > 0) {
                const ctrlKey = params.event && params.event.srcEvent && (params.event.srcEvent.ctrlKey || params.event.srcEvent.metaKey);

                if (!ctrlKey) {
                    this.ui.clearSelection();
                }
            }
        });
    }

    /**
     * Capture original network data
     */
    captureOriginalData() {
        this.originalNetworkEdgeData = this.network.body.data.edges.get();
        this.originalNodeData = this.network.body.data.nodes.get();


        // Build graph adjacency list for pathfinding
        this.buildGraphAdjacency();
    }

    /**
     * Build graph adjacency list
     */
    buildGraphAdjacency() {
        const originalEdges = this.networkData.originalEdgeData
            .filter(edge => edge.begin !== edge.end) // Filter out self-edges
            .map((edge, i) => ({
                from: edge.begin,
                to: edge.end,
                id: i
            }));

        originalEdges.forEach(edge => {
            // Skip self-edges (double-check)
            if (edge.from === edge.to) return;

            if (!this.graphAdjacency.has(edge.from)) {
                this.graphAdjacency.set(edge.from, []);
            }
            if (!this.graphAdjacency.has(edge.to)) {
                this.graphAdjacency.set(edge.to, []);
            }
            this.graphAdjacency.get(edge.from).push(edge.to);
            this.graphAdjacency.get(edge.to).push(edge.from); // Bidirectional for pathfinding
        });

        // Create edge data for filtering (filter out self-edges)
        this.edgeData = this.networkData.originalEdgeData
            .filter(edge => edge.begin !== edge.end) // Filter out self-edges
            .map((edge, i) => ({
                id: i,
                count: edge.count || 1,
                type: edge.type
            }));

        // Make edge data globally available for stats
        window.edgeData = this.edgeData;
    }

    /**
     * Setup application components
     */
    setupComponents() {
        this.filters = new NetworkFilters(this.networkData, this.network);

        // Debug logging for Node 1045 in filters
        console.log('\n🔍 [DEBUG] NetworkFilters Setup - Node 1045 Check:');
        console.log('=================================================');
        if (this.filters.pathwayFilter && this.filters.pathwayFilter.bestPathBarrierMap) {
            const barrier = this.filters.pathwayFilter.bestPathBarrierMap[1045];
            if (barrier !== undefined) {
                console.log('✅ Node 1045 barrier in pathwayFilter:', barrier, 'kcal/mol');
            } else {
                console.log('❌ Node 1045 NOT found in pathwayFilter.bestPathBarrierMap');
            }
        } else {
            console.log('❌ pathwayFilter or bestPathBarrierMap not available');
        }
        console.log('=================================================');

        this.ui = new NetworkUI(this.networkData, this.network);
        this.stats = new NetworkStats(this.network);

        // Initialize network data viewer
        this.dataViewer = new NetworkDataViewer(this);

        // Initialize pathway profile viewer
        this.pathwayViewer = new PathwayProfileViewer(this.networkData);

        // Initialize components with data
        this.filters.initialize(this.originalNodeData, this.originalNetworkEdgeData, this.edgeData, this.graphAdjacency, null);
        this.ui.initialize(this.originalNodeData, this.originalNetworkEdgeData, this.graphAdjacency, null);

        // Make components globally available
        window.networkDataViewer = this.dataViewer;
        window.pathwayViewer = this.pathwayViewer;

    }

    /**
     * Initialize UI after components are ready
     */
    initializeUI() {
        // Set up barrier threshold slider
        const barrierSlider = document.getElementById("barrierThreshold");
        barrierSlider.min = 0;
        barrierSlider.value = CONFIG.DEFAULT_MAX_BARRIER_THRESHOLD || 50;
        document.getElementById("barrierThresholdValue").innerText = CONFIG.DEFAULT_MAX_BARRIER_THRESHOLD || 50;

        // Initialize graph statistics display
        this.stats.updateGraphStats();

        // Apply initial filtering
        this.applyFilters();

        // Update initial tooltips
        setTimeout(() => {
            NetworkUtils.updateNodeTooltipsOnly(this.network, this.networkData);
            NetworkUtils.updateEdgeTooltips(this.network, this.networkData);
        }, 100);
    }

    /**
     * Setup control event handlers
     */
    setupControlEvents() {
        // Threshold sliders
        document.getElementById("threshold").addEventListener("input", (e) => {
            const val = parseInt(e.target.value);
            document.getElementById("thresholdValue").innerText = val;
            this.applyFilters();
        });

        document.getElementById("barrierThreshold").addEventListener("input", (e) => {
            const val = parseFloat(e.target.value);
            document.getElementById("barrierThresholdValue").innerText = val;
            this.applyFilters();
        });

        // Disable all filters checkbox
        document.getElementById("disableAllFilters").addEventListener("change", () => {
            this.applyFilters();
        });

        // Make global functions available
        window.toggleReactionProfile = () => this.ui.toggleReactionProfile();
        window.showGridModal = () => this.ui.showGridModal();
        window.hideGridModal = () => this.ui.hideGridModal();
        window.hideHelpModal = () => this.ui.hideHelpModal();
        window.hideIndicesModal = () => this.ui.hideIndicesModal();
        window.showChartModal = () => this.ui.showChartModal();
    }

    /**
     * Apply all filters and update network
     */
    applyFilters() {
        if (!this.filters) return;

        const result = this.filters.applyAllFilters();
        if (result) {
            this.updateNetworkState(
                result.finalVisibleNodeIds,
                result.finalVisibleEdgeIds,
                result.currentNetworkEdgeData,
                result.currentEdgeData
            );

            // If maxReactionEnergies are available, store them in the networkData
            if (result.maxReactionEnergies) {
                this.networkData.maxReactionEnergies = result.maxReactionEnergies;

                // Update node tooltips to show max reaction energies
                NetworkUtils.updateNodeTooltipsOnly(this.network, this.networkData);
            }

            // If maxBarrierMap is available, store it for potential debugging/display
            if (result.maxBarrierMap) {
                this.networkData.maxBarrierMap = result.maxBarrierMap;
            }

            // If best path results are available, store them in the networkData
            if (result.bestPathResults) {
                console.log('Updating network with best path results:', result.bestPathResults);

                // Store the best path data in the network data object
                this.networkData.bestPathEnergyMap = result.bestPathResults.bestPathEnergyMap || new Map();
                this.networkData.bestPathBarrierMap = result.bestPathResults.bestPathBarrierMap || new Map();
                this.networkData.bestPathMap = result.bestPathResults.bestPathMap || new Map();

                // Log some statistics about the calculated paths
                const totalNodes = result.bestPathResults.bestPathEnergyMap ? result.bestPathResults.bestPathEnergyMap.size : 0;
                console.log(`Best path calculation complete: ${totalNodes} nodes processed`);

                // Update node tooltips to potentially show best path energies (if UI supports it)
                NetworkUtils.updateNodeTooltipsOnly(this.network, this.networkData);
            }

            this.stats.updateGraphStats();
        }
    }

    /**
     * Update network state with filtered data
     */
    updateNetworkState(finalVisibleNodeIds, finalVisibleEdgeIds, currentNetworkEdgeData, currentEdgeData) {


        const nodes = this.network.body.data.nodes;
        const edges = this.network.body.data.edges;

        const currentNodes = nodes.get();
        const currentEdges = edges.get();
        const currentNodeIds = new Set(currentNodes.map(n => n.id));
        const currentEdgeIds = new Set(currentEdges.map(e => e.id));

        // Determine changes needed
        const nodesToRemove = [];
        const nodesToAdd = [];
        const edgesToRemove = [];
        const edgesToAdd = [];

        // Find nodes to remove/add
        currentNodeIds.forEach(nodeId => {
            if (!finalVisibleNodeIds.has(nodeId)) {
                nodesToRemove.push(nodeId);
            }
        });

        finalVisibleNodeIds.forEach(nodeId => {
            if (!currentNodeIds.has(nodeId)) {
                const nodeData = NetworkUtils.createUpdatedNodeData(
                    nodeId,
                    this.originalNodeData,
                    this.networkData,
                    this.ui.selectedNodes
                );
                if (nodeData) {
                    nodesToAdd.push(nodeData);

                }
            }
        });

        // Find edges to remove/add
        currentEdgeIds.forEach(edgeId => {
            if (!finalVisibleEdgeIds.has(edgeId)) {
                edgesToRemove.push(edgeId);
            }
        });

        finalVisibleEdgeIds.forEach(edgeId => {
            if (!currentEdgeIds.has(edgeId)) {
                // Use current edge data (which might be consolidated)
                const edgeInfo = currentNetworkEdgeData.find(e => e.id === edgeId);
                const edgeDataInfo = currentEdgeData.find(e => e.id === edgeId);
                if (edgeInfo && edgeDataInfo) {
                    // Create edge with updated tooltip showing consolidated count
                    const edgeWithTooltip = {
                        ...edgeInfo,
                        title: this.getEdgeTitle(edgeDataInfo),
                        width: this.getEdgeWidth(edgeDataInfo)
                    };
                    edgesToAdd.push(edgeWithTooltip);
                }
            }
        });

        // Apply changes in correct order: edges first, then nodes
        if (edgesToRemove.length > 0) {
            edges.remove(edgesToRemove);

        }

        if (nodesToRemove.length > 0) {
            nodes.remove(nodesToRemove);

        }

        if (nodesToAdd.length > 0) {
            nodes.add(nodesToAdd);

        }

        if (edgesToAdd.length > 0) {
            edges.add(edgesToAdd);

        }

        // Validate result
        NetworkUtils.validateNetworkState(this.network);

        // Update tooltips after network state changes
        setTimeout(() => {
            NetworkUtils.updateNodeTooltipsOnly(this.network, this.networkData);
            NetworkUtils.updateEdgeTooltips(this.network, this.networkData);
        }, 50);
    }

    /**
     * Select pathway nodes from substrate to target node
     */
    selectPathwayToNode(nodeId) {
        const substrateId = parseInt(this.networkData.substrateId);
        const targetNodeId = parseInt(nodeId);

        console.log(`🛤️ Selecting pathway from substrate (${substrateId}) to node ${targetNodeId}`);

        // If clicking on substrate, just select it
        if (targetNodeId === substrateId) {
            this.ui.clearSelection();
            this.ui.toggleNodeSelection(substrateId, false, false);
            return;
        }

        // Get the best path for this node
        let pathwayNodes = [];

        // Try to get pathway from bestPathMap first
        if (this.networkData.bestPathMap && this.networkData.bestPathMap[targetNodeId]) {
            pathwayNodes = this.networkData.bestPathMap[targetNodeId];
            console.log(`✅ Found best path from bestPathMap:`, pathwayNodes);
        }
        // Fallback to precomputedPathways
        else if (this.networkData.precomputedPathways && this.networkData.precomputedPathways[targetNodeId]) {
            const pathways = this.networkData.precomputedPathways[targetNodeId];
            if (pathways.length > 0) {
                // Use the first (best) pathway
                pathwayNodes = pathways[0].path || [];
                console.log(`✅ Found path from precomputedPathways:`, pathwayNodes);
            }
        }

        // If no pathway found, just select the target node
        if (!pathwayNodes || pathwayNodes.length === 0) {
            console.log(`⚠️ No pathway found for node ${targetNodeId}, selecting target node only`);
            this.ui.clearSelection();
            this.ui.toggleNodeSelection(targetNodeId, false, false);
            return;
        }

        // Clear current selection and select pathway nodes in order
        this.ui.clearSelection();

        // Select nodes in pathway order
        pathwayNodes.forEach((pathNodeId, index) => {
            const nodeIdInt = parseInt(pathNodeId);
            if (index === 0) {
                // First node (substrate) - clear previous and start new selection
                this.ui.toggleNodeSelection(nodeIdInt, false, false);
            } else {
                // Subsequent nodes - add to selection maintaining order
                this.ui.toggleNodeSelection(nodeIdInt, true, true);
            }
        });

        console.log(`🎯 Selected pathway with ${pathwayNodes.length} nodes: [${pathwayNodes.join(' → ')}]`);

        // Show chart modal for the pathway if more than one node
        if (pathwayNodes.length > 1) {
            setTimeout(() => {
                this.ui.showChartModal();
            }, 100);
        }
    }

    /**
     * Show pathway profile for a selected node
     */
    showPathwayProfile(nodeId) {
        if (this.pathwayViewer) {
            this.pathwayViewer.showPathwayProfile(nodeId);
        }
    }

}

// Initialize the application when data is available
window.initializeNetworkApp = function (initialData) {
    console.log('\n🔍 [DEBUG] Network App Initialization - Node 1045 Analysis:');
    console.log('========================================================');

    // Check if node 1045 exists in data
    if (initialData.bestPathBarrierMap && initialData.bestPathBarrierMap[1045] !== undefined) {
        console.log('✅ Node 1045 found in bestPathBarrierMap:', initialData.bestPathBarrierMap[1045], 'kcal/mol');
    } else {
        console.log('❌ Node 1045 NOT found in bestPathBarrierMap');
    }

    if (initialData.precomputedPathways && initialData.precomputedPathways[1045]) {
        console.log('✅ Node 1045 found in precomputedPathways');
        console.log('   Number of pathways:', initialData.precomputedPathways[1045].length);
    } else {
        console.log('❌ Node 1045 NOT found in precomputedPathways');
    }

    console.log('========================================================');

    const app = new NetworkApp();
    app.initialize(initialData);

    // Make app globally available for debugging
    window.networkApp = app;
};

// Global functions for HTML onclick handlers
window.hideGridModal = function () {
    if (window.networkApp && window.networkApp.ui) {
        window.networkApp.ui.hideGridModal();
    }
};

window.hideStereoisomerModal = function () {
    console.log("Stereoisomer modal functionality has been removed");
};

window.hideHelpModal = function () {
    if (window.networkApp && window.networkApp.ui) {
        window.networkApp.ui.hideHelpModal();
    }
};

window.toggleReactionProfile = function () {
    if (window.networkApp && window.networkApp.ui) {
        window.networkApp.ui.toggleReactionProfile();
    }
};

window.showChartModal = function () {
    if (window.networkApp && window.networkApp.ui) {
        window.networkApp.ui.showChartModal();
    }
};

window.showPathwayProfile = function (nodeId) {
    if (window.pathwayViewer) {
        window.pathwayViewer.showPathwayProfile(nodeId);
    }
};

window.hidePathwayViewer = function () {
    if (window.pathwayViewer) {
        window.pathwayViewer.hideViewer();
    }
};

// Global function for pathway profile display
window.showPathwayProfile = function (nodeId) {
    if (window.networkApp) {
        window.networkApp.showPathwayProfile(nodeId);
    }
};

// Global function for pathway selection
window.selectPathwayToNode = function (nodeId) {
    if (window.networkApp) {
        window.networkApp.selectPathwayToNode(nodeId);
    }
};

// Global function to test reaction profile display for a given path
window.showReactionProfileForPath = function (pathNodes) {
    if (window.pathwayViewer) {
        console.log(`🎯 Displaying reaction profile for path: [${pathNodes.join(' → ')}]`);
        window.pathwayViewer.showReactionProfileForPath(pathNodes);
    } else {
        console.warn('PathwayViewer not available');
    }
};

// Global function to simulate node selection for testing
window.selectNodesInSequence = function (nodeIds) {
    if (window.networkApp && window.networkApp.ui) {
        const ui = window.networkApp.ui;

        // Clear current selection
        ui.clearSelection();

        // Select nodes in sequence
        nodeIds.forEach((nodeId, index) => {
            ui.selectedNodes.add(nodeId);
            ui.selectedNodesInOrder.push(nodeId);
            ui.updateNodeColor(nodeId, CONFIG.SELECTION_COLOR);
            if (index === nodeIds.length - 1) {
                ui.lastSelectedNode = nodeId;
            }
        });

        ui.updateSelectedEdges();
        ui.updateSelectionDisplay();
        ui.checkAndUpdateReactionProfile();

        console.log(`✅ Selected nodes in sequence: [${nodeIds.join(' → ')}]`);
    } else {
        console.warn('NetworkApp UI not available');
    }
};
