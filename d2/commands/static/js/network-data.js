/**
 * Network data management and stereoisomer functionality
 */
class NetworkData {
    constructor(initialData) {
        this.svgMap = initialData.svgMap;
        this.molEnergyMap = initialData.molEnergyMap;
        this.molChargeMap = initialData.chargeMap;
        this.molWeightMap = initialData.weightMap;
        this.originalEdgeData = initialData.originalEdgeData;
        this.substrateId = initialData.substrateId;
        this.energyType = initialData.energyType;
        this.stereoisomerGroups = initialData.stereoisomerGroups || {};
        this.nodeLabels = initialData.nodeLabels || {};

        // Use pre-calculated relative energies from Python
        this.relEnergyMap = initialData.relEnergyMap || {};

        // Enhanced pathway data from unified ReactionNetwork (Phase 2)
        this.precomputedPathways = initialData.precomputedPathways || {};
        this.pathwayLookup = initialData.pathwayLookup || {};
        this.pathwayProfiles = initialData.pathwayProfiles || {};
        this.pathwayBarriers = initialData.pathwayBarriers || {};
        this.configData = initialData.configData || {};
        this.advancedFeatures = initialData.advancedFeatures || {};

        // Best path energy data (from Python or calculated dynamically in JS)
        this.bestPathEnergyMap = initialData.bestPathEnergyMap || {};
        this.bestPathBarrierMap = initialData.bestPathBarrierMap || {};
        this.bestPathMap = initialData.bestPathMap || {};

        // Initialize enhanced path analyzer with Python pathways (Phase 2.1)
        this.initializePathAnalyzer();
    }

    /**
     * Initialize path analyzer with pre-computed Python pathways
     */
    initializePathAnalyzer() {
        if (window.pathAnalyzer) {
            console.log("🚀 Initializing PathAnalyzer with Python pathway data...");
            window.pathAnalyzer.initializeWithPythonPathways(this);

            // Log pathway data availability
            if (this.advancedFeatures.hasReactionProfiles) {
                console.log(`✅ Loaded ${this.advancedFeatures.pathwayCount} pre-computed pathways`);
                console.log(`✅ ${this.advancedFeatures.reachableNodes} reachable nodes`);
            }

            if (this.advancedFeatures.hasTransitionStates) {
                console.log("✅ Transition state data available");
            }

            if (this.advancedFeatures.hasSteReoisomers) {
                console.log(`✅ ${this.advancedFeatures.totalSteReoisomerGroups} stereoisomer groups`);
            }
        } else {
            console.warn("⚠️  PathAnalyzer not available - loading order issue?");
        }
    }

    /**
     * Get stereoisomer group information
     */
    getStereoisomerGroups() {
        // Group nodes by their stereoisomer group ID
        const groups = {};
        Object.keys(this.stereoisomerGroups).forEach(nodeId => {
            const groupId = this.stereoisomerGroups[nodeId];
            if (!groups[groupId]) {
                groups[groupId] = [];
            }
            groups[groupId].push(parseInt(nodeId));
        });
        return groups;
    }

    /**
     * Get representative node for a stereoisomer group (lowest energy)
     */
    getStereoisomerRepresentative(groupNodes) {
        if (groupNodes.length === 1) {
            return groupNodes[0];
        }

        // Find the node with lowest relative energy, or substrate if in group
        let representative = groupNodes[0];
        let lowestEnergy = this.relEnergyMap[representative] || 0;

        for (const nodeId of groupNodes) {
            if (nodeId === this.substrateId) {
                return nodeId; // Substrate always wins
            }
            const energy = this.relEnergyMap[nodeId] || 0;
            if (energy < lowestEnergy) {
                lowestEnergy = energy;
                representative = nodeId;
            }
        }

        return representative;
    }

    /**
     * Get collapsed node data for stereoisomer groups
     */
    getCollapsedNodeData(originalNodes) {
        const groups = this.getStereoisomerGroups();
        const collapsedNodes = [];
        const processedGroups = new Set();

        originalNodes.forEach(node => {
            const groupId = this.stereoisomerGroups[node.id];

            if (processedGroups.has(groupId)) {
                return; // Already processed this group
            }

            processedGroups.add(groupId);
            const groupNodes = groups[groupId] || [node.id];
            const representative = this.getStereoisomerRepresentative(groupNodes);

            // Create collapsed node
            const collapsedNode = { ...node };
            collapsedNode.id = representative;

            if (groupNodes.length > 1) {
                // Update label to show it's a group
                collapsedNode.label = node.id === this.substrateId ?
                    'Substrate' : `${representative} (+${groupNodes.length - 1})`;
                collapsedNode.title = this.getCollapsedNodeTitle(groupNodes, representative);
                collapsedNode.stereoisomerGroup = groupNodes;
                collapsedNode.isCollapsed = true;
            }

            collapsedNodes.push(collapsedNode);
        });

        return collapsedNodes;
    }

    /**
     * Get title for collapsed stereoisomer node
     */
    getCollapsedNodeTitle(groupNodes, representative) {
        const energies = groupNodes.map(id => ({
            id: id,
            energy: this.relEnergyMap[id] || 0
        })).sort((a, b) => a.energy - b.energy);

        let title = `Stereoisomer Group (${groupNodes.length} isomers)\n`;
        title += `Representative: Node ${representative}\n\n`;
        title += 'All isomers:\n';

        energies.forEach(({ id, energy }) => {
            const marker = id === representative ? '★ ' : '  ';
            title += `${marker}Node ${id}: ∆${this.energyType} ${energy.toFixed(2)} kcal/mol\n`;
        });

        return title;
    }
}

/**
 * Update network data with best path energies based on filtered edges
 * @param {Object} networkData - The network data object
 * @param {Array} filteredEdges - Currently visible edges after filtering
 */
function updateNetworkDataWithBestPaths(networkData, filteredEdges) {
    if (!window.pathAnalyzer || !networkData || !filteredEdges) {
        console.warn('Path analyzer, network data, or filtered edges not available');
        return null;
    }

    // Calculate best path energies based on current filtered edges
    const pathResults = window.pathAnalyzer.updateNetworkWithBestPaths(
        networkData,
        filteredEdges
    );

    return pathResults;
}

// Export the function for use in other modules
window.updateNetworkDataWithBestPaths = updateNetworkDataWithBestPaths;
