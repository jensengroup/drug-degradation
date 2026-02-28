/**
 * Enhanced Network filtering system with Barrier-Based Filtering using Precomputed Profiles
 */

/**
 * NetworkPathwayFilter - Implements barrier-based filtering using precomputed reaction profiles
 * This class uses barrier_height from precomputed pathways for energy-based filtering
 */
class NetworkPathwayFilter {
    constructor(networkData) {
        this.networkData = networkData;
        this.precomputedPathways = networkData.precomputedPathways || {};
        this.originalEdgeData = networkData.originalEdgeData || [];
        this.bestPathBarrierMap = networkData.bestPathBarrierMap || {};
        this.pathwayAnalyzer = window.pathAnalyzer;

        console.log(`🔧 NetworkPathwayFilter initialized`);
        console.log(`   Nodes with precomputed barriers: ${Object.keys(this.bestPathBarrierMap).length}`);
        console.log(`   Nodes with pathway data: ${Object.keys(this.precomputedPathways).length}`);
        console.log(`   Total edges: ${this.originalEdgeData.length}`);

        // Create bestPathBarrierMap if it doesn't exist
        if (Object.keys(this.bestPathBarrierMap).length === 0 && Object.keys(this.precomputedPathways).length > 0) {
            console.log('🔧 Creating bestPathBarrierMap from precomputedPathways...');
            this.bestPathBarrierMap = this.createBestPathBarrierMap();
        }
    }

    /**
     * Create bestPathBarrierMap from precomputedPathways
     * @returns {Object} Map of nodeId -> minimum barrier height
     */
    createBestPathBarrierMap() {
        const barrierMap = {};
        for (const [nodeIdStr, pathways] of Object.entries(this.precomputedPathways)) {
            const nodeId = parseInt(nodeIdStr);
            const barriers = pathways
                .map(p => p.barrier_height)
                .filter(b => b !== null && b !== undefined && !isNaN(b));

            if (barriers.length > 0) {
                barrierMap[nodeId] = Math.min(...barriers);
            }
        }
        console.log(`   Created barrier map for ${Object.keys(barrierMap).length} nodes`);
        return barrierMap;
    }

    /**
     * Filter network using precomputed barriers from bestPathBarrierMap
     * @param {number} maxBarrier - Maximum barrier height in kcal/mol
     * @returns {Object} Valid nodes and edges from nodes below barrier cutoff
     */
    filterByBarrierHeight(maxBarrier) {
        console.log(`🔍 Filtering with barrier cutoff: ${maxBarrier} kcal/mol`);

        const validNodes = new Set();
        const validEdges = new Set();
        const validPathways = {};
        let totalNodes = 0;
        let filteredNodes = 0;

        // Debug for specific node (like 1045)
        const debugNodeId = 1045;

        // Check if debug node exists in bestPathBarrierMap before filtering
        if (this.bestPathBarrierMap[debugNodeId] !== undefined) {
            console.log(`🔍 PRE-FILTER CHECK: Node ${debugNodeId} barrier: ${this.bestPathBarrierMap[debugNodeId]} kcal/mol`);
            console.log(`   Should pass filter (≤${maxBarrier}): ${this.bestPathBarrierMap[debugNodeId] <= maxBarrier}`);
        }

        // Filter nodes based on bestPathBarrierMap
        for (const [nodeIdStr, barrierHeight] of Object.entries(this.bestPathBarrierMap)) {
            const nodeId = parseInt(nodeIdStr);
            totalNodes++;

            // Debug logging for specific node
            if (nodeId === debugNodeId) {
                console.log(`🔍 PROCESSING Node ${debugNodeId}:`);
                console.log(`   - Barrier: ${barrierHeight} kcal/mol`);
                console.log(`   - Cutoff: ${maxBarrier} kcal/mol`);
                console.log(`   - Comparison: ${barrierHeight} <= ${maxBarrier} = ${barrierHeight <= maxBarrier}`);
            }

            // Keep nodes with barriers below cutoff
            if (barrierHeight <= maxBarrier) {
                validNodes.add(nodeId);
                filteredNodes++;

                // Include pathway data if available
                if (this.precomputedPathways[nodeId]) {
                    validPathways[nodeId] = this.precomputedPathways[nodeId];

                    // Collect edges from all pathways for this valid node
                    for (const pathway of this.precomputedPathways[nodeId]) {
                        if (pathway.path && pathway.path.length > 1) {
                            for (let i = 0; i < pathway.path.length - 1; i++) {
                                const edgeKey = `${pathway.path[i]}-${pathway.path[i + 1]}`;
                                validEdges.add(edgeKey);
                            }
                        }
                    }
                }

                if (nodeId === debugNodeId) {
                    console.log(`   ✅ Node ${debugNodeId} KEPT: Added to validNodes set`);
                    console.log(`   - validNodes.has(${debugNodeId}): ${validNodes.has(debugNodeId)}`);
                    console.log(`   - Pathway data available: ${this.precomputedPathways[nodeId] ? 'YES' : 'NO'}`);
                }
            } else if (nodeId === debugNodeId) {
                console.log(`   ❌ Node ${debugNodeId} FILTERED OUT: ${barrierHeight} > ${maxBarrier}`);
            }
        }

        // Collect valid edge objects
        const validEdgeObjects = [];
        for (const edge of this.originalEdgeData) {
            const edgeKey = `${edge.begin}-${edge.end}`;
            if (validEdges.has(edgeKey)) {
                validEdgeObjects.push(edge);
            }
        }

        // Final check for debug node in results
        console.log(`🔍 POST-FILTER CHECK: Node ${debugNodeId} in results?`);
        console.log(`   - In validNodes set: ${validNodes.has(debugNodeId)}`);
        console.log(`   - In validNodes array: ${Array.from(validNodes).includes(debugNodeId)}`);
        console.log(`   - In validPathways: ${validPathways[debugNodeId] ? 'YES' : 'NO'}`);

        console.log(`   Results: ${filteredNodes}/${totalNodes} nodes pass filter`);
        console.log(`   Valid nodes: ${validNodes.size}, Valid edges: ${validEdgeObjects.length}`);

        // Convert to old format for compatibility with existing visualization code
        const filteredNodesCompat = {};
        for (const nodeId of validNodes) {
            const pathways = validPathways[nodeId] || [];
            if (pathways.length > 0) {
                const bestPathway = pathways[0]; // Assume first is best
                filteredNodesCompat[nodeId.toString()] = {
                    bestPathway: bestPathway,
                    allPathways: pathways,
                    barrierHeight: this.bestPathBarrierMap[nodeId] || Number.POSITIVE_INFINITY,
                    pathLength: bestPathway.path_length || (bestPathway.path ? bestPathway.path.length - 1 : 0),
                    reactionPath: bestPathway.path || [],
                    energyProfile: bestPathway.profile || []
                };
            }
        }

        // Final validation log
        console.log(`🔍 FINAL RESULT CHECK: Node ${debugNodeId} in final results?`, Array.from(validNodes).includes(debugNodeId));
        if (Array.from(validNodes).includes(debugNodeId)) {
            console.log(`   ✅ SUCCESS: Node ${debugNodeId} is in the final results!`);
        } else {
            console.log(`   ❌ PROBLEM: Node ${debugNodeId} is missing from final results despite passing filter!`);
        }

        return {
            filteredNodes: filteredNodesCompat,
            statistics: {
                totalNodes,
                passedNodes: filteredNodes,
                filterRate: (filteredNodes / totalNodes * 100).toFixed(1) + '%'
            }
        };
    }

    /**
     * Filter nodes by combined criteria: barrier height and edge count
     * @param {Object} criteria - {maxBarrier, minCount}
     * @returns {Object} Filtered nodes meeting all criteria
     */
    filterByCombinedCriteria(criteria) {
        const { maxBarrier, minCount } = criteria;
        console.log(`🔍 Combined filtering: barrier ≤ ${maxBarrier}, count ≥ ${minCount || 1}`);

        // Step 1: Filter by barrier height first
        const barrierFiltered = this.filterByBarrierHeight(maxBarrier);

        if (!minCount || minCount <= 1) {
            // No edge count filtering needed
            return barrierFiltered;
        }

        // Step 2: Apply edge count filtering
        const validEdgeKeys = new Set();

        for (const edge of this.originalEdgeData) {
            if (this.edgeMeetsCountRequirement(edge, minCount)) {
                validEdgeKeys.add(`${edge.begin}-${edge.end}`);
            }
        }

        // Step 3: Filter barrier-passed nodes by edge count
        const finalFiltered = {};
        let edgeFilterPassed = 0;

        for (const [nodeId, nodeData] of Object.entries(barrierFiltered.filteredNodes)) {
            const pathway = nodeData.bestPathway;

            if (this.pathwayUsesValidEdges(pathway, validEdgeKeys)) {
                finalFiltered[nodeId] = nodeData;
                edgeFilterPassed++;
            }
        }

        console.log(`   Barrier filter: ${barrierFiltered.statistics.passedNodes} nodes`);
        console.log(`   Edge filter: ${edgeFilterPassed} nodes`);

        return {
            filteredNodes: finalFiltered,
            statistics: {
                totalNodes: barrierFiltered.statistics.totalNodes,
                barrierPassed: barrierFiltered.statistics.passedNodes,
                finalPassed: edgeFilterPassed,
                filterRate: (edgeFilterPassed / barrierFiltered.statistics.totalNodes * 100).toFixed(1) + '%'
            }
        };
    }

    /**
     * Legacy filter method for compatibility
     * @param {number} minCount - Minimum count threshold
     * @returns {Object} Filtered pathways and edges
     */
    filterByMinCount(minCount) {
        console.log(`🔍 Filtering by minimum count: ${minCount}`);

        // First filter edges using exact Python logic
        const validEdges = [];
        const validEdgeKeys = new Set();

        if (this.networkData.originalEdgeData) {
            for (const edge of this.networkData.originalEdgeData) {
                if (this.edgeMeetsCountRequirement(edge, minCount)) {
                    validEdges.push(edge);
                    validEdgeKeys.add(`${edge.begin}-${edge.end}`);
                }
            }
        }

        // Filter pathways using Python's exact approach: sort first, then filter, then take top 5
        const filteredPathways = {};

        if (this.pathwayAnalyzer && this.pathwayAnalyzer.precomputedPathways) {
            for (const [nodeIdStr, pathways] of Object.entries(this.pathwayAnalyzer.precomputedPathways)) {
                const nodeId = parseInt(nodeIdStr);

                // Sort pathways using Python's order: (path_length, barrier_height)
                const sortedPathways = this.sortPathwaysPythonStyle(pathways);
                const validPathways = [];

                // Filter sorted pathways (Python approach)
                for (const pathway of sortedPathways) {
                    if (this.pathwayUsesValidEdges(pathway, validEdgeKeys)) {
                        validPathways.push(pathway);

                        // Take only top 5 like Python does
                        if (validPathways.length >= 5) {
                            break;
                        }
                    }
                }

                if (validPathways.length > 0) {
                    filteredPathways[nodeId] = validPathways;
                }
            }
        }

        console.log(`✅ Count filtering: ${validEdges.length} valid edges, ${Object.keys(filteredPathways).length} reachable nodes`);

        return {
            filteredPathways,
            validEdges,
            reachableNodes: Object.keys(filteredPathways).length
        };
    }

    /**
     * Filter pathways by maximum barrier height (Phase 3.1)
     * @param {number} maxBarrier - Maximum barrier height in kcal/mol
     * @returns {Object} Filtered pathways
     */
    filterByMaxBarrier(maxBarrier) {
        console.log(`🔍 Filtering by maximum barrier: ${maxBarrier} kcal/mol`);

        const filteredPathways = {};

        if (this.pathwayAnalyzer && this.pathwayAnalyzer.precomputedPathways) {
            for (const [nodeIdStr, pathways] of Object.entries(this.pathwayAnalyzer.precomputedPathways)) {
                const nodeId = parseInt(nodeIdStr);
                const validPathways = [];

                for (const pathway of pathways) {
                    const barrierHeight = pathway.barrier_height || 0;
                    if (barrierHeight <= maxBarrier) {
                        validPathways.push(pathway);
                    }
                }

                if (validPathways.length > 0) {
                    filteredPathways[nodeId] = validPathways;
                }
            }
        }

        console.log(`✅ Barrier filtering: ${Object.keys(filteredPathways).length} nodes with barriers ≤ ${maxBarrier} kcal/mol`);

        return {
            filteredPathways,
            reachableNodes: Object.keys(filteredPathways).length
        };
    }

    /**
     * Filter pathways by maximum path length (Phase 3.1)
     * @param {number} maxLength - Maximum number of reaction steps
     * @returns {Object} Filtered pathways
     */
    filterByMaxPathLength(maxLength) {
        console.log(`🔍 Filtering by maximum path length: ${maxLength} steps`);

        const filteredPathways = {};

        if (this.pathwayAnalyzer && this.pathwayAnalyzer.precomputedPathways) {
            for (const [nodeIdStr, pathways] of Object.entries(this.pathwayAnalyzer.precomputedPathways)) {
                const nodeId = parseInt(nodeIdStr);
                const validPathways = [];

                for (const pathway of pathways) {
                    const pathLength = pathway.path_length || (pathway.path ? pathway.path.length - 1 : 0);
                    if (pathLength <= maxLength) {
                        validPathways.push(pathway);
                    }
                }

                if (validPathways.length > 0) {
                    filteredPathways[nodeId] = validPathways;
                }
            }
        }

        console.log(`✅ Length filtering: ${Object.keys(filteredPathways).length} nodes with paths ≤ ${maxLength} steps`);

        return {
            filteredPathways,
            reachableNodes: Object.keys(filteredPathways).length
        };
    }

    /**
     * Combined filtering with multiple criteria using Python's exact algorithm
     * @param {Object} criteria - Filtering criteria {minCount, maxBarrier, maxLength}
     * @returns {Object} Filtered pathways meeting all criteria
     */
    filterByCombinedCriteria(criteria) {
        console.log(`🔍 Combined filtering (Python-style):`, criteria);

        // First determine valid edges based on count requirement
        const validEdgeKeys = new Set();

        if (criteria.minCount !== undefined && this.networkData.originalEdgeData) {
            for (const edge of this.networkData.originalEdgeData) {
                if (this.edgeMeetsCountRequirement(edge, criteria.minCount)) {
                    validEdgeKeys.add(`${edge.begin}-${edge.end}`);
                }
            }
        } else {
            // If no count filter, all edges are valid
            if (this.networkData.originalEdgeData) {
                for (const edge of this.networkData.originalEdgeData) {
                    validEdgeKeys.add(`${edge.begin}-${edge.end}`);
                }
            }
        }

        // Filter pathways using Python's exact approach
        const filteredPathways = {};

        if (this.pathwayAnalyzer && this.pathwayAnalyzer.precomputedPathways) {
            for (const [nodeIdStr, pathways] of Object.entries(this.pathwayAnalyzer.precomputedPathways)) {
                const nodeId = parseInt(nodeIdStr);

                // Sort pathways first (Python approach)
                const sortedPathways = this.sortPathwaysPythonStyle(pathways);
                const validPathways = [];

                // Apply all filters on sorted pathways
                for (const pathway of sortedPathways) {
                    // Check count filter (via edge validation)
                    if (!this.pathwayUsesValidEdges(pathway, validEdgeKeys)) {
                        continue;
                    }

                    // Check barrier threshold
                    if (criteria.maxBarrier !== undefined &&
                        (pathway.barrier_height || Number.POSITIVE_INFINITY) > criteria.maxBarrier) {
                        continue;
                    }

                    // Check length threshold
                    if (criteria.maxLength !== undefined) {
                        const pathLength = pathway.path_length || (pathway.path ? pathway.path.length - 1 : 0);
                        if (pathLength > criteria.maxLength) {
                            continue;
                        }
                    }

                    validPathways.push(pathway);

                    // Keep only top 5 paths per node (Python behavior)
                    if (validPathways.length >= 5) {
                        break;
                    }
                }

                if (validPathways.length > 0) {
                    filteredPathways[nodeId] = validPathways;
                }
            }
        }

        console.log(`✅ Combined filtering: ${Object.keys(filteredPathways).length} nodes meet all criteria`);

        return {
            filteredPathways,
            reachableNodes: Object.keys(filteredPathways).length,
            criteria
        };
    }

    /**
     * Combined filtering with multiple criteria (Phase 3.1) - DEPRECATED
     * Use filterByCombinedCriteria for exact Python matching
     * @param {Object} criteria - Filtering criteria {minCount, maxBarrier, maxLength}
     * @returns {Object} Filtered pathways meeting all criteria
     */
    filterByCombinedCriteriaOld(criteria) {
        console.log(`🔍 Combined filtering:`, criteria);

        let filteredPathways = {};

        // Start with all pathways if no specific criteria
        if (this.pathwayAnalyzer && this.pathwayAnalyzer.precomputedPathways) {
            filteredPathways = { ...this.pathwayAnalyzer.precomputedPathways };
        }

        // Apply count filtering first (affects edge availability)
        if (criteria.minCount !== undefined) {
            const countResult = this.filterByMinCount(criteria.minCount);
            filteredPathways = countResult.filteredPathways;
        }

        // Apply barrier filtering
        if (criteria.maxBarrier !== undefined) {
            const validPathways = {};
            for (const [nodeIdStr, pathways] of Object.entries(filteredPathways)) {
                const nodeId = parseInt(nodeIdStr);
                const barrierValidPathways = pathways.filter(pathway =>
                    (pathway.barrier_height || 0) <= criteria.maxBarrier
                );

                if (barrierValidPathways.length > 0) {
                    validPathways[nodeId] = barrierValidPathways;
                }
            }
            filteredPathways = validPathways;
        }

        // Apply length filtering
        if (criteria.maxLength !== undefined) {
            const validPathways = {};
            for (const [nodeIdStr, pathways] of Object.entries(filteredPathways)) {
                const nodeId = parseInt(nodeIdStr);
                const lengthValidPathways = pathways.filter(pathway => {
                    const pathLength = pathway.path_length || (pathway.path ? pathway.path.length - 1 : 0);
                    return pathLength <= criteria.maxLength;
                });

                if (lengthValidPathways.length > 0) {
                    validPathways[nodeId] = lengthValidPathways;
                }
            }
            filteredPathways = validPathways;
        }

        console.log(`✅ Combined filtering: ${Object.keys(filteredPathways).length} nodes meet all criteria`);

        return {
            filteredPathways,
            reachableNodes: Object.keys(filteredPathways).length,
            criteria
        };
    }

    /**
     * Check if an edge meets count requirement (Phase 3.2)
     * Ports exact Python logic: reaction edges check count, others pass
     * @param {Object} edge - Edge object
     * @param {number} minCount - Minimum count threshold
     * @returns {boolean} True if edge meets requirement
     */
    edgeMeetsCountRequirement(edge, minCount) {
        const edgeType = (edge.type || "").toLowerCase();

        // Non-reaction edges always pass the count requirement (exact Python logic)
        if (edgeType !== "reaction") {
            return true;
        }

        // Reaction edges must meet the count threshold
        const edgeCount = edge.count || 1;
        return edgeCount >= minCount;
    }

    /**
     * Sort pathways using Python's exact sorting logic: path_length first, then barrier_height
     * @param {Array} pathways - Array of pathway objects
     * @returns {Array} Sorted pathways
     */
    sortPathwaysPythonStyle(pathways) {
        return [...pathways].sort((a, b) => {
            // Python sorts by (path_length, barrier_height)
            const aLength = a.path_length || 0;
            const bLength = b.path_length || 0;

            if (aLength !== bLength) {
                return aLength - bLength;
            }

            const aBarrier = a.barrier_height || Number.POSITIVE_INFINITY;
            const bBarrier = b.barrier_height || Number.POSITIVE_INFINITY;

            return aBarrier - bBarrier;
        });
    }

    /**
     * Check if a pathway uses only valid edges
     * @param {Object} pathway - Pathway object with path array
     * @param {Set} validEdgeKeys - Set of valid edge keys
     * @returns {boolean} True if all edges in path are valid
     */
    pathwayUsesValidEdges(pathway, validEdgeKeys) {
        if (!pathway.path || pathway.path.length < 2) {
            return false;
        }

        // Check each step in the pathway
        for (let i = 0; i < pathway.path.length - 1; i++) {
            const startId = pathway.path[i];
            const endId = pathway.path[i + 1];
            const edgeKey = `${startId}-${endId}`;

            if (!validEdgeKeys.has(edgeKey)) {
                return false;
            }
        }

        return true;
    }

    /**
     * Get visualization data for filtered nodes
     * @param {Object} filteredResult - Result from filtering methods
     * @returns {Object} Data ready for network visualization
     */
    getVisualizationData(filteredResult) {
        const nodes = [];
        const edges = [];
        const pathwayData = {};

        // Process filtered nodes
        for (const [nodeId, nodeData] of Object.entries(filteredResult.filteredNodes)) {
            const numNodeId = parseInt(nodeId);

            // Add node for visualization
            nodes.push({
                id: numNodeId,
                barrier: nodeData.barrierHeight,
                pathLength: nodeData.pathLength,
                hasPathway: true,
                energyProfile: nodeData.energyProfile
            });

            // Store pathway data
            pathwayData[nodeId] = {
                path: nodeData.reactionPath,
                profile: nodeData.energyProfile,
                barrier: nodeData.barrierHeight
            };

            // Add edges from the pathway
            const path = nodeData.reactionPath;
            for (let i = 0; i < path.length - 1; i++) {
                const sourceId = path[i];
                const targetId = path[i + 1];

                // Find the original edge data
                const originalEdge = this.originalEdgeData.find(e =>
                    e.begin === sourceId && e.end === targetId
                );

                if (originalEdge) {
                    edges.push({
                        source: sourceId,
                        target: targetId,
                        type: originalEdge.type || 'reaction',
                        count: originalEdge.count || 1,
                        inPathway: true
                    });
                }
            }
        }

        return {
            nodes,
            edges,
            pathwayData,
            statistics: filteredResult.statistics
        };
    }
}

/**
 * Legacy NetworkFilters class with enhanced pathway integration (Phase 3.3)
 */
class NetworkFilters {
    constructor(networkData, network) {
        this.networkData = networkData;
        this.network = network;
        this.originalNodeData = [];
        this.originalNetworkEdgeData = [];
        this.edgeData = [];
        this.graphAdjacency = new Map();
        this.stereoisomerManager = null;

        // Phase 3.3: Enhanced pathway filtering integration
        this.pathwayFilter = new NetworkPathwayFilter(networkData);
    }

    /**
     * Initialize filter system with network data (enhanced for Phase 3.3)
     */
    initialize(originalNodeData, originalNetworkEdgeData, edgeData, graphAdjacency, stereoisomerManager = null) {
        this.originalNodeData = originalNodeData;
        this.originalNetworkEdgeData = originalNetworkEdgeData;
        this.edgeData = edgeData;
        this.graphAdjacency = graphAdjacency;
        this.stereoisomerManager = stereoisomerManager;

        // If we have a stereoisomer manager, set up the original edges
        if (this.stereoisomerManager && typeof this.stereoisomerManager.setOriginalEdges === 'function') {
            this.stereoisomerManager.setOriginalEdges(originalNetworkEdgeData, edgeData);
        } else if (this.stereoisomerManager) {
            console.error('StereoisomerManager is missing setOriginalEdges method');
        }

        // Initialize pathway filter with current data
        this.pathwayFilter = new NetworkPathwayFilter(this.networkData);
    }

    /**
     * Enhanced filtering with pathway-based criteria (Phase 3.3)
     * @param {Object} criteria - Combined filtering criteria
     * @returns {Object} Filtered network state
     */
    applyEnhancedFiltering(criteria) {
        console.log("🔄 Applying enhanced pathway-based filtering...", criteria);

        let filteredResults = null;

        // Apply pathway-based filtering if criteria provided
        if (criteria.minCount !== undefined || criteria.maxBarrier !== undefined || criteria.maxLength !== undefined) {
            filteredResults = this.pathwayFilter.filterByCombinedCriteria(criteria);
        }

        // Convert pathway filter results to network filter format
        if (filteredResults) {
            const filteredEdges = this.convertPathwaysToEdges(filteredResults.filteredPathways);
            const visibleNodeIds = this.getNodesFromPathways(filteredResults.filteredPathways);

            return {
                passingNodeIds: new Set(visibleNodeIds),
                failedNodeIds: new Set(),
                filteredEdges: filteredEdges,
                pathwayResults: filteredResults
            };
        }

        // No filtering applied
        return {
            passingNodeIds: new Set(Object.keys(this.networkData.relEnergyMap).map(id => parseInt(id))),
            failedNodeIds: new Set(),
            filteredEdges: this.originalNetworkEdgeData || [],
            pathwayResults: null
        };
    }

    /**
     * Convert filtered pathways to edge list for network display
     * @param {Object} filteredPathways - Filtered pathway data
     * @returns {Array} Array of edges used in filtered pathways
     */
    convertPathwaysToEdges(filteredPathways) {
        const usedEdgeKeys = new Set();
        const usedEdges = [];

        // Collect all edges used in filtered pathways
        for (const [nodeId, pathways] of Object.entries(filteredPathways)) {
            for (const pathway of pathways) {
                if (pathway.path && pathway.path.length > 1) {
                    for (let i = 0; i < pathway.path.length - 1; i++) {
                        const edgeKey = `${pathway.path[i]}-${pathway.path[i + 1]}`;
                        usedEdgeKeys.add(edgeKey);
                    }
                }
            }
        }

        // Find corresponding edge objects
        if (this.networkData.originalEdgeData) {
            for (const edge of this.networkData.originalEdgeData) {
                const edgeKey = `${edge.begin}-${edge.end}`;
                if (usedEdgeKeys.has(edgeKey)) {
                    usedEdges.push(edge);
                }
            }
        }

        return usedEdges;
    }

    /**
     * Get all nodes involved in filtered pathways
     * @param {Object} filteredPathways - Filtered pathway data
     * @returns {Array} Array of node IDs
     */
    getNodesFromPathways(filteredPathways) {
        const nodeIds = new Set();

        // Add substrate
        nodeIds.add(this.networkData.substrateId);

        // Add all nodes from pathways
        for (const [nodeId, pathways] of Object.entries(filteredPathways)) {
            nodeIds.add(parseInt(nodeId));

            for (const pathway of pathways) {
                if (pathway.path) {
                    for (const pathNodeId of pathway.path) {
                        nodeIds.add(pathNodeId);
                    }
                }
            }
        }

        return Array.from(nodeIds);
    }

    /**
     * Filter edges by count threshold
     */
    filterEdgesByCount(threshold) {
        const passingEdgeIds = new Set();
        const failedEdgeIds = new Set();

        this.originalNetworkEdgeData.forEach(edgeInfo => {
            const edgeDataInfo = this.edgeData.find(e => e.id === edgeInfo.id);
            const isReactionEdge = edgeDataInfo && edgeDataInfo.type.toLowerCase() === "reaction";

            if (!isReactionEdge || edgeDataInfo.count >= threshold) {
                passingEdgeIds.add(edgeInfo.id);
            } else {
                failedEdgeIds.add(edgeInfo.id);
            }
        });

        return { passingEdgeIds, failedEdgeIds };
    }

    /**
     * Filter edges by count threshold using custom edge data
     */
    filterEdgesByCountWithEdgeData(threshold, networkEdgeData, edgeData) {
        const passingEdgeIds = new Set();
        const failedEdgeIds = new Set();

        networkEdgeData.forEach(edgeInfo => {
            const edgeDataInfo = edgeData.find(e => e.id === edgeInfo.id);
            const isReactionEdge = edgeDataInfo && edgeDataInfo.type.toLowerCase() === "reaction";

            if (!isReactionEdge || edgeDataInfo.count >= threshold) {
                passingEdgeIds.add(edgeInfo.id);
            } else {
                failedEdgeIds.add(edgeInfo.id);
            }
        });

        return { passingEdgeIds, failedEdgeIds };
    }

    /**
     * Filter edges by maximum barrier threshold
     * Calculates max barrier for each edge as: barrier - min_energy_along_path_to_substrate
     */
    filterEdgesByMaxBarrier(threshold, visibleEdgeIds) {
        const passingEdgeIds = new Set();
        const failedEdgeIds = new Set();
        const maxBarrierMap = {};

        // If no edges provided, return empty sets
        if (!visibleEdgeIds || visibleEdgeIds.size === 0) {
            return { passingEdgeIds, failedEdgeIds, maxBarrierMap };
        }

        // Calculate max barriers for all edges
        const edgeMaxBarriers = this.calculateMaxBarriersForEdges(visibleEdgeIds);

        // Filter edges based on barrier threshold
        for (const edgeId of visibleEdgeIds) {
            const maxBarrier = edgeMaxBarriers[edgeId];

            // Only apply filter to edges that have a barrier value
            if (maxBarrier !== undefined && maxBarrier !== null) {
                maxBarrierMap[edgeId] = maxBarrier;

                if (Math.abs(maxBarrier) <= threshold) {
                    passingEdgeIds.add(edgeId);
                } else {
                    failedEdgeIds.add(edgeId);
                }
            } else {
                // If no barrier data, include the edge (e.g., for conformer edges)
                passingEdgeIds.add(edgeId);
            }
        }

        return { passingEdgeIds, failedEdgeIds, maxBarrierMap };
    }

    /**
     * Filter nodes by maximum barrier using bestPathBarrierMap (our working approach)
     * @param {number} threshold - Maximum barrier threshold in kcal/mol
     * @param {Set} visibleNodeIds - Set of visible node IDs to filter
     * @returns {Object} Object containing passing nodes and barrier map
     */
    filterNodesByMaxBarrier(threshold, visibleNodeIds) {
        console.log(`🔍 [DEBUG] filterNodesByMaxBarrier called with ${visibleNodeIds.size} nodes`);
        console.log(`   Node 1045 in input visibleNodeIds: ${visibleNodeIds.has(1045)}`);
        if (!visibleNodeIds.has(1045)) {
            console.log(`   ❌ Node 1045 was ALREADY FILTERED OUT before barrier filter!`);
        }

        const passingNodeIds = new Set();
        const failedNodeIds = new Set();
        const maxBarrierMap = {};

        // Get bestPathBarrierMap from network data or pathwayFilter
        let bestPathBarrierMap = this.networkData.bestPathBarrierMap;

        // If not available in networkData, try to get it from pathwayFilter
        if (!bestPathBarrierMap && this.pathwayFilter && this.pathwayFilter.bestPathBarrierMap) {
            bestPathBarrierMap = this.pathwayFilter.bestPathBarrierMap;
        }

        // If still not available, create it from precomputedPathways
        if (!bestPathBarrierMap && this.networkData.precomputedPathways) {
            console.log('🔧 Creating bestPathBarrierMap for barrier filtering...');
            bestPathBarrierMap = {};
            for (const [nodeIdStr, pathways] of Object.entries(this.networkData.precomputedPathways)) {
                const nodeId = parseInt(nodeIdStr);
                const barriers = pathways
                    .map(p => p.barrier_height)
                    .filter(b => b !== null && b !== undefined && !isNaN(b));

                if (barriers.length > 0) {
                    bestPathBarrierMap[nodeId] = Math.min(...barriers);
                }
            }
            console.log(`   Created barrier map for ${Object.keys(bestPathBarrierMap).length} nodes`);
        }

        if (!bestPathBarrierMap) {
            console.warn('⚠️ No bestPathBarrierMap available for barrier filtering - including all nodes');
            return { passingNodeIds: visibleNodeIds, failedNodeIds, maxBarrierMap };
        }

        // Debug for specific node (like 1045)
        const debugNodeId = 1045;
        if (bestPathBarrierMap[debugNodeId] !== undefined) {
            console.log(`🔍 Node ${debugNodeId} barrier check: ${bestPathBarrierMap[debugNodeId]} kcal/mol vs ${threshold} kcal/mol`);
            console.log(`   Should pass: ${bestPathBarrierMap[debugNodeId] <= threshold}`);
        }

        // Filter nodes based on their best pathway barriers
        for (const nodeId of visibleNodeIds) {
            const nodeIdNum = parseInt(nodeId);
            const barrierHeight = bestPathBarrierMap[nodeIdNum];

            if (barrierHeight !== undefined && barrierHeight !== null) {
                maxBarrierMap[nodeIdNum] = barrierHeight;

                if (barrierHeight <= threshold) {
                    passingNodeIds.add(nodeIdNum);

                    if (nodeIdNum === debugNodeId) {
                        console.log(`   ✅ Node ${debugNodeId} PASSED barrier filter`);
                    }
                } else {
                    failedNodeIds.add(nodeIdNum);

                    if (nodeIdNum === debugNodeId) {
                        console.log(`   ❌ Node ${debugNodeId} FAILED barrier filter`);
                    }
                }
            } else {
                // If no barrier data, include the node
                passingNodeIds.add(nodeIdNum);
            }
        }

        console.log(`🔍 Barrier filtering: ${passingNodeIds.size}/${visibleNodeIds.size} nodes passed (threshold: ${threshold} kcal/mol)`);

        // Final debug check for specific node
        if (debugNodeId && visibleNodeIds.has(debugNodeId)) {
            console.log(`🔍 Final barrier filter result for node ${debugNodeId}: ${passingNodeIds.has(debugNodeId) ? 'INCLUDED' : 'EXCLUDED'}`);
        }

        return { passingNodeIds, failedNodeIds, maxBarrierMap };
    }

    /**
     * Calculate maximum barriers for all edges
     * Max barrier = barrier - min_energy_along_path_from_substrate_to_source_node
     */
    calculateMaxBarriersForEdges(visibleEdgeIds) {
        const maxBarrierMap = {};
        const substrateId = parseInt(this.networkData.substrateId);
        const relEnergyMap = this.networkData.relEnergyMap;

        // Calculate minimum energy paths from substrate to all nodes
        const minEnergyPaths = this.calculateMinEnergyPathsFromSubstrate(visibleEdgeIds);

        // Process each edge
        for (const edgeId of visibleEdgeIds) {
            const edge = this.originalNetworkEdgeData.find(e => e.id === edgeId);
            const edgeData = this.edgeData.find(e => e.id === edgeId);

            if (!edge || !edgeData) continue;

            // Check for barrier information in multiple locations
            let barrierValue = edge.barrier;
            if (barrierValue === undefined || barrierValue === null) {
                barrierValue = edgeData.barrier;
            }
            if (barrierValue === undefined || barrierValue === null && edge.originalData) {
                barrierValue = edge.originalData.barrier;
            }

            // Only process edges that have barrier information
            if (barrierValue === undefined || barrierValue === null) {
                continue;
            }

            const fromId = parseInt(edge.from);
            const toId = parseInt(edge.to);

            // Get the minimum energy along the path from substrate to the source node
            const minEnergyToSource = minEnergyPaths[fromId];

            // If source node is not reachable from substrate, skip this edge
            if (minEnergyToSource === undefined) {
                continue;
            }

            // Calculate max barrier: barrier - min_energy_along_path_to_source
            const maxBarrier = barrierValue - minEnergyToSource;
            maxBarrierMap[edgeId] = maxBarrier;
        }

        return maxBarrierMap;
    }

    /**
     * Calculate minimum energy along paths from substrate to all reachable nodes
     * Uses best path energies calculated by PathAnalyzer for unified energy calculations
     * Returns a map of nodeId -> minimum_energy_along_path_from_substrate
     */
    calculateMinEnergyPathsFromSubstrate(visibleEdgeIds) {
        const substrateId = parseInt(this.networkData.substrateId);

        // Use best path energies from PathAnalyzer - these already represent
        // the optimal energy to reach each node from the substrate
        const bestPathEnergyMap = this.networkData.bestPathEnergyMap;

        if (!bestPathEnergyMap) {
            throw new Error('bestPathEnergyMap is required but not available. Ensure PathAnalyzer has run first.');
        }

        // Build set of nodes that are connected via visible edges
        const connectedNodes = new Set();
        for (const edgeId of visibleEdgeIds) {
            const edge = this.originalNetworkEdgeData.find(e => e.id === edgeId);
            if (!edge) continue;

            connectedNodes.add(parseInt(edge.from));
            connectedNodes.add(parseInt(edge.to));
        }

        // For nodes with best path data and connected via visible edges, use the best path energy
        // The best path energies already represent the minimum energy barrier needed
        // to reach each node from the substrate via the optimal path
        const minEnergyPaths = {};
        for (const nodeId of connectedNodes) {
            if (bestPathEnergyMap[nodeId] !== undefined) {
                minEnergyPaths[nodeId] = bestPathEnergyMap[nodeId];
            }
        }

        return minEnergyPaths;
    }

    /**
     * Determine if an edge is directed away from substrate
     */
    isEdgeAwayFromSubstrate(fromNode, toNode) {
        const fromEnergy = this.networkData.relEnergyMap[fromNode] || 0;
        const toEnergy = this.networkData.relEnergyMap[toNode] || 0;
        return toEnergy > fromEnergy;
    }

    /**
     * Apply bidirectional filtering rule
     */
    applyBidirectionalEdgeFiltering(passingEdgeIds, failedEdgeIds) {
        const finalPassingEdgeIds = new Set(passingEdgeIds);

        // Build adjacency graph from ALL passing edges
        const forwardAdjacency = new Map(); // substrate -> other nodes

        this.originalNetworkEdgeData.forEach(edgeInfo => {
            if (passingEdgeIds.has(edgeInfo.id)) {
                // Add to forward adjacency
                if (!forwardAdjacency.has(edgeInfo.from)) {
                    forwardAdjacency.set(edgeInfo.from, new Set());
                }
                forwardAdjacency.get(edgeInfo.from).add(edgeInfo.to);
            }
        });

        // Function to check if there's any path from substrate to target node
        const hasPathFromSubstrate = (targetNode) => {
            if (targetNode === this.networkData.substrateId) return true;

            const visited = new Set();
            const queue = [this.networkData.substrateId];
            visited.add(this.networkData.substrateId);

            while (queue.length > 0) {
                const currentNode = queue.shift();

                if (forwardAdjacency.has(currentNode)) {
                    for (const nextNode of forwardAdjacency.get(currentNode)) {
                        if (nextNode === targetNode) return true;

                        if (!visited.has(nextNode)) {
                            visited.add(nextNode);
                            queue.push(nextNode);
                        }
                    }
                }
            }
            return false;
        };

        // Filter out edges that lead TO substrate but have no forward path FROM substrate
        this.originalNetworkEdgeData.forEach(edgeInfo => {
            if (passingEdgeIds.has(edgeInfo.id) && edgeInfo.to === this.networkData.substrateId) {
                // This is an edge leading to the substrate (X -> substrate)
                // Check if there's any path from substrate to X
                if (!hasPathFromSubstrate(edgeInfo.from)) {
                    finalPassingEdgeIds.delete(edgeInfo.id);
                }
            }
        });

        return finalPassingEdgeIds;
    }

    /**
     * Apply bidirectional filtering rule using custom edge data
     */
    applyBidirectionalEdgeFilteringWithEdgeData(passingEdgeIds, failedEdgeIds, networkEdgeData, currentEdgeData) {
        const finalPassingEdgeIds = new Set(passingEdgeIds);

        // Build adjacency graph from ALL passing edges (regardless of count threshold)
        const forwardAdjacency = new Map(); // substrate -> other nodes
        const backwardAdjacency = new Map(); // other nodes -> substrate

        networkEdgeData.forEach(edgeInfo => {
            if (passingEdgeIds.has(edgeInfo.id)) {
                // Add to forward adjacency
                if (!forwardAdjacency.has(edgeInfo.from)) {
                    forwardAdjacency.set(edgeInfo.from, new Set());
                }
                forwardAdjacency.get(edgeInfo.from).add(edgeInfo.to);

                // Add to backward adjacency
                if (!backwardAdjacency.has(edgeInfo.to)) {
                    backwardAdjacency.set(edgeInfo.to, new Set());
                }
                backwardAdjacency.get(edgeInfo.to).add(edgeInfo.from);
            }
        });

        // Function to check if there's any path from substrate to target node
        const hasPathFromSubstrate = (targetNode) => {
            if (targetNode === this.networkData.substrateId) return true;

            const visited = new Set();
            const queue = [this.networkData.substrateId];
            visited.add(this.networkData.substrateId);

            while (queue.length > 0) {
                const currentNode = queue.shift();

                if (forwardAdjacency.has(currentNode)) {
                    for (const nextNode of forwardAdjacency.get(currentNode)) {
                        if (nextNode === targetNode) return true;

                        if (!visited.has(nextNode)) {
                            visited.add(nextNode);
                            queue.push(nextNode);
                        }
                    }
                }
            }
            return false;
        };

        // Filter out edges that lead TO substrate but have no forward path FROM substrate
        networkEdgeData.forEach(edgeInfo => {
            if (passingEdgeIds.has(edgeInfo.id) && edgeInfo.to === this.networkData.substrateId) {
                // This is an edge leading to the substrate (X -> substrate)
                // Check if there's any path from substrate to X
                if (!hasPathFromSubstrate(edgeInfo.from)) {
                    finalPassingEdgeIds.delete(edgeInfo.id);
                }
            }
        });

        return finalPassingEdgeIds;
    }

    /**
     * Filter edges by node visibility
     */
    filterEdgesByNodeVisibility(visibleNodeIds) {
        const visibleEdgeIds = new Set();

        this.originalNetworkEdgeData.forEach(edgeInfo => {
            const fromVisible = visibleNodeIds.has(edgeInfo.from);
            const toVisible = visibleNodeIds.has(edgeInfo.to);

            if (fromVisible && toVisible) {
                visibleEdgeIds.add(edgeInfo.id);
            }
        });

        return visibleEdgeIds;
    }

    /**
     * Filter edges by node visibility using custom edge data
     */
    filterEdgesByNodeVisibilityWithEdgeData(visibleNodeIds, networkEdgeData) {
        const visibleEdgeIds = new Set();

        networkEdgeData.forEach(edgeInfo => {
            const fromVisible = visibleNodeIds.has(edgeInfo.from);
            const toVisible = visibleNodeIds.has(edgeInfo.to);

            if (fromVisible && toVisible) {
                visibleEdgeIds.add(edgeInfo.id);
            }
        });

        return visibleEdgeIds;
    }

    /**
     * Find nodes connected to substrate using BFS
     */
    findConnectedNodes(substrateId, visibleEdges) {
        // Build temporary adjacency list from visible edges only
        const tempAdjacency = new Map();
        visibleEdges.forEach(edge => {
            if (!tempAdjacency.has(edge.from)) {
                tempAdjacency.set(edge.from, []);
            }
            if (!tempAdjacency.has(edge.to)) {
                tempAdjacency.set(edge.to, []);
            }
            tempAdjacency.get(edge.from).push(edge.to);
            tempAdjacency.get(edge.to).push(edge.from); // Bidirectional
        });

        // BFS to find all connected nodes
        const connected = new Set();
        const queue = [substrateId];
        connected.add(substrateId);

        while (queue.length > 0) {
            const currentNode = queue.shift();
            const neighbors = tempAdjacency.get(currentNode) || [];

            for (const neighbor of neighbors) {
                if (!connected.has(neighbor)) {
                    connected.add(neighbor);
                    queue.push(neighbor);
                }
            }
        }

        return connected;
    }

    /**
     * Filter by connectivity
     */
    filterByConnectivity(visibleNodeIds, connectedNodes, showDisconnected) {
        const finalVisibleNodeIds = new Set();

        visibleNodeIds.forEach(nodeId => {
            const isConnected = connectedNodes.has(nodeId);
            const isSubstrate = nodeId == this.networkData.substrateId;

            if (isConnected || (showDisconnected && !isSubstrate)) {
                finalVisibleNodeIds.add(nodeId);
            }
        });

        // Ensure substrate is always included if it passed energy filter
        if (visibleNodeIds.has(this.networkData.substrateId)) {
            finalVisibleNodeIds.add(this.networkData.substrateId);
        }

        return finalVisibleNodeIds;
    }

    /**
     * Get filter parameters from UI
     */
    getFilterParameters() {
        return {
            countThreshold: parseInt(document.getElementById("threshold").value),
            maxBarrierThreshold: parseFloat(document.getElementById("barrierThreshold").value),
            disableAllFilters: document.getElementById('disableAllFilters').checked
        };
    }

    /**
     * Main filtering coordinator
     */
    applyAllFilters() {
        console.log('\n🔍 [DEBUG] Starting applyAllFilters - Node 1045 tracking:');
        console.log('=======================================================');

        if (this.originalNodeData.length === 0 || this.originalNetworkEdgeData.length === 0) {
            console.warn("Original data not yet available");
            return;
        }

        // Get filter parameters
        const filterParams = this.getFilterParameters();

        // Create a set of all node IDs
        const allNodeIds = new Set(this.originalNodeData.map(node => parseInt(node.id)));
        console.log(`📋 Total nodes: ${allNodeIds.size}, Node 1045 in all nodes: ${allNodeIds.has(1045)}`);

        // Set up current edge data
        let currentNetworkEdgeData = this.originalNetworkEdgeData;
        let currentEdgeData = this.edgeData;

        // If "Disable all filters" is checked, show all nodes and edges
        if (filterParams.disableAllFilters) {
            console.log('🚫 All filters disabled - showing everything');
            // Filter out self-edges only
            const selfEdgeFilteredIds = this.filterSelfEdges(
                new Set(currentNetworkEdgeData.map(e => e.id)),
                currentNetworkEdgeData
            );

            // Calculate best path energies for unfiltered network
            let bestPathResults = null;
            if (window.updateNetworkDataWithBestPaths) {
                try {
                    // Create edge list for path analysis (all edges except self-edges)
                    const filteredEdgeList = currentNetworkEdgeData.filter(edge =>
                        selfEdgeFilteredIds.has(edge.id)
                    );

                    // Calculate best paths with unfiltered network
                    bestPathResults = window.updateNetworkDataWithBestPaths(
                        this.networkData,
                        filteredEdgeList
                    );
                } catch (error) {
                    console.warn('Error calculating best path energies:', error);
                }
            }

            return {
                finalVisibleNodeIds: allNodeIds,
                finalVisibleEdgeIds: selfEdgeFilteredIds,
                currentNetworkEdgeData,
                currentEdgeData,
                bestPathResults
            };
        }

        // Handle stereoisomer filtering (if enabled)
        if (this.stereoisomerManager && this.stereoisomerManager.collapsed) {
            console.log('🔄 Applying stereoisomer consolidation...');
            const consolidatedEdges = this.stereoisomerManager.consolidateEdgesForVisibleNodes([...allNodeIds]);
            if (consolidatedEdges) {
                currentNetworkEdgeData = consolidatedEdges.networkEdgeData;
                currentEdgeData = consolidatedEdges.edgeData;
            }
            this.stereoisomerManager.collapseStereoisomers([...allNodeIds]);
        }

        // Apply count filter only (remove energy filter)
        console.log('🔄 Applying count filter...');
        const { passingEdgeIds: countFilteredEdgeIds } = this.filterEdgesByCount(filterParams.countThreshold);

        // Build visible nodes set from count-filtered edges
        const countFilteredEdges = this.originalNetworkEdgeData.filter(edge =>
            countFilteredEdgeIds.has(edge.id)
        );

        const countFilteredNodeIds = new Set();
        countFilteredEdges.forEach(edge => {
            countFilteredNodeIds.add(parseInt(edge.from));
            countFilteredNodeIds.add(parseInt(edge.to));
        });

        // Always include substrate
        countFilteredNodeIds.add(parseInt(this.networkData.substrateId));

        console.log(`   After count filter: ${countFilteredNodeIds.size} nodes, Node 1045: ${countFilteredNodeIds.has(1045)}`);

        // Calculate best path energies BEFORE barrier filtering so they're available
        let bestPathResults = null;
        if (window.updateNetworkDataWithBestPaths) {
            try {
                // Create filtered edge list for initial path analysis
                const initialFilteredEdgeList = currentNetworkEdgeData.filter(edge =>
                    countFilteredEdgeIds.has(edge.id)
                );

                // Calculate best paths with count-filtered network
                bestPathResults = window.updateNetworkDataWithBestPaths(
                    this.networkData,
                    initialFilteredEdgeList
                );
            } catch (error) {
                console.warn('Error calculating initial best path energies:', error);
            }
        }

        // Apply barrier filter using our working node-based approach
        // Use bestPathBarrierMap for filtering nodes, not edges
        console.log('🔄 Applying barrier filter...');
        const { passingNodeIds: barrierFilteredNodeIds, maxBarrierMap } = this.filterNodesByMaxBarrier(
            filterParams.maxBarrierThreshold,
            countFilteredNodeIds  // Use count-filtered nodes instead of energy-filtered
        );
        console.log(`   After barrier filter: ${barrierFilteredNodeIds.size} nodes, Node 1045: ${barrierFilteredNodeIds.has(1045)}`);

        // Apply stereoisomer filtering to barrier-filtered nodes (if enabled)
        console.log('🔄 Applying stereoisomer node filtering...');
        const energyFilteredNodeIds = this.stereoisomerManager
            ? this.stereoisomerManager.filterVisibleNodes(barrierFilteredNodeIds)
            : barrierFilteredNodeIds;
        console.log(`   After stereoisomer filter: ${energyFilteredNodeIds.size} nodes, Node 1045: ${energyFilteredNodeIds.has(1045)}`);

        // Filter edges to only include those between barrier-filtered nodes
        const barrierFilteredEdgeIds = new Set();
        for (const edgeId of countFilteredEdgeIds) {  // Use count-filtered edges
            const edge = currentNetworkEdgeData.find(e => e.id === edgeId);
            if (edge && barrierFilteredNodeIds.has(edge.from) && barrierFilteredNodeIds.has(edge.to)) {
                barrierFilteredEdgeIds.add(edgeId);
            }
        }

        // Apply directed reachability filter as the final step
        console.log('🔄 Applying directed reachability filter...');
        const { reachableNodeIds: directedReachableNodeIds, filteredEdgeIds: directedFilteredEdgeIds } =
            this.filterByDirectedReachability(energyFilteredNodeIds, barrierFilteredEdgeIds);
        console.log(`   After directed reachability: ${directedReachableNodeIds.size} nodes, Node 1045: ${directedReachableNodeIds.has(1045)}`);

        // Final edge filtering to ensure no self-edges
        const finalVisibleEdgeIds = new Set();
        directedFilteredEdgeIds.forEach(edgeId => {
            const edge = currentNetworkEdgeData.find(e => e.id === edgeId);
            if (edge && directedReachableNodeIds.has(edge.from) && directedReachableNodeIds.has(edge.to) && edge.from !== edge.to) {
                finalVisibleEdgeIds.add(edgeId);
            }
        });

        console.log(`🏁 FINAL RESULT: ${directedReachableNodeIds.size} nodes, Node 1045: ${directedReachableNodeIds.has(1045)}`);
        console.log('=======================================================');

        // Update best path energies after all filtering is complete
        if (window.updateNetworkDataWithBestPaths) {
            try {
                // Create filtered edge list for path analysis
                const filteredEdgeList = currentNetworkEdgeData.filter(edge =>
                    finalVisibleEdgeIds.has(edge.id)
                );

                // Recalculate best paths with final filtered network
                bestPathResults = window.updateNetworkDataWithBestPaths(
                    this.networkData,
                    filteredEdgeList
                );
            } catch (error) {
                console.warn('Error calculating final best path energies:', error);
            }
        }

        return {
            finalVisibleNodeIds: directedReachableNodeIds,
            finalVisibleEdgeIds,
            currentNetworkEdgeData,
            currentEdgeData,
            maxReactionEnergies: {},  // Empty object since energy filter was removed
            maxBarrierMap,
            bestPathResults
        };
    }

    /**
     * Filter out self-edges
     */
    filterSelfEdges(edgeIds, networkEdgeData) {
        const filteredEdgeIds = new Set();

        edgeIds.forEach(edgeId => {
            const edge = networkEdgeData.find(e => e.id === edgeId);
            if (edge && edge.from !== edge.to) {
                filteredEdgeIds.add(edgeId);
            }
        });

        return filteredEdgeIds;
    }

    /**
     * Filter nodes to only include those reachable from the substrate via directed paths
     * @param {Set} nodeIds - Set of node IDs to filter
     * @param {Set} edgeIds - Set of edge IDs to use for reachability check
     * @returns {Object} Object containing filtered nodes and edges
     */
    filterByDirectedReachability(nodeIds, edgeIds) {
        const substrateId = parseInt(this.networkData.substrateId);
        console.log(`🔍 [DEBUG] Directed reachability - Substrate: ${substrateId}, Input nodes: ${nodeIds.size}, Node 1045 in input: ${nodeIds.has(1045)}`);

        // Build a directed graph using only the specified edges
        const directedGraph = new Map();

        // Initialize graph for all nodes
        nodeIds.forEach(nodeId => {
            directedGraph.set(parseInt(nodeId), []);
        });

        // Add edges to the graph (in forward direction only)
        let node1045Edges = [];
        this.originalNetworkEdgeData.forEach(edge => {
            if (edgeIds.has(edge.id)) {
                const fromId = parseInt(edge.from);
                const toId = parseInt(edge.to);

                if (directedGraph.has(fromId)) {
                    directedGraph.get(fromId).push(toId);
                }

                // Track edges involving node 1045
                if (fromId === 1045 || toId === 1045) {
                    node1045Edges.push(`${fromId} -> ${toId}`);
                }
            }
        });

        if (nodeIds.has(1045)) {
            console.log(`   Node 1045 edges in filtered graph: [${node1045Edges.join(', ')}]`);
            console.log(`   Node 1045 neighbors: [${(directedGraph.get(1045) || []).join(', ')}]`);
        }

        // Find all nodes reachable from substrate using BFS
        const reachableNodes = new Set([substrateId]);
        const queue = [substrateId];

        while (queue.length > 0) {
            const currentNode = queue.shift();

            // Process all neighbors
            const neighbors = directedGraph.get(currentNode) || [];
            for (const neighborId of neighbors) {
                if (!reachableNodes.has(neighborId)) {
                    reachableNodes.add(neighborId);
                    queue.push(neighborId);

                    // Log if we reach node 1045
                    if (neighborId === 1045) {
                        console.log(`   ✅ Node 1045 reached from ${currentNode}!`);
                    }
                }
            }
        }

        console.log(`   Final reachable nodes: ${reachableNodes.size}, Node 1045 reachable: ${reachableNodes.has(1045)}`);

        if (nodeIds.has(1045) && !reachableNodes.has(1045)) {
            console.log(`   ❌ Node 1045 FILTERED OUT by directed reachability!`);
        }

        while (queue.length > 0) {
            const currentNode = queue.shift();

            // Process all neighbors
            const neighbors = directedGraph.get(currentNode) || [];
            for (const neighborId of neighbors) {
                if (!reachableNodes.has(neighborId)) {
                    reachableNodes.add(neighborId);
                    queue.push(neighborId);
                }
            }
        }

        // Filter edges to only include those between reachable nodes
        const filteredEdgeIds = new Set();
        edgeIds.forEach(edgeId => {
            const edge = this.originalNetworkEdgeData.find(e => e.id === edgeId);
            if (edge && reachableNodes.has(parseInt(edge.from)) && reachableNodes.has(parseInt(edge.to))) {
                filteredEdgeIds.add(edgeId);
            }
        });

        return {
            reachableNodeIds: reachableNodes,
            filteredEdgeIds: filteredEdgeIds
        };
    }
}
