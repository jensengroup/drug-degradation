/**
 * Enhanced Path Analysis Module with Python Integration
 *
 * This module now leverages pre-computed pathways from the Python ReactionNetwork
 * for improved performance and consistency. The BFS algorithms are maintained
 * for backward compatibility and edge filtering scenarios.
 *
 * Key features:
 * - Uses pre-computed Python pathways when available (Phase 2 enhancement)
 * - Falls back to JavaScript BFS for filtered edge scenarios
 * - Exact API compatibility maintained
 * - Enhanced with reaction profile and transition state support
 */

class PathAnalyzer {
    constructor() {
        this.edgeEnergies = new Map();
        this.results = new Map();
        this.precomputedPathways = null;  // Will store Python pathways
        this.pathwayProfiles = null;      // Pre-computed reaction profiles
        this.pathwayBarriers = null;      // Pre-computed barriers
    }

    /**
     * Initialize with pre-computed pathway data from Python
     * @param {Object} networkData - Network data with enhanced pathway information
     */
    initializeWithPythonPathways(networkData) {
        // Store pre-computed pathway data from unified ReactionNetwork
        this.precomputedPathways = networkData.precomputedPathways || {};
        this.pathwayProfiles = networkData.pathwayProfiles || {};
        this.pathwayBarriers = networkData.pathwayBarriers || {};

        console.log(`✅ Loaded ${Object.keys(this.precomputedPathways).length} pre-computed pathways from Python`);

        // Build pathway lookup for fast access
        this.pathwayLookup = this.buildPathwayLookup(networkData);
    }

    /**
     * Build optimized pathway lookup structure
     * @param {Object} networkData - Network data
     * @returns {Object} Optimized lookup structure
     */
    buildPathwayLookup(networkData) {
        const lookup = {};

        if (networkData.pathwayLookup) {
            // Use pre-built lookup from Python if available
            return networkData.pathwayLookup;
        }

        // Build lookup from precomputed pathways
        for (const [nodeId, pathways] of Object.entries(this.precomputedPathways)) {
            if (pathways && pathways.length > 0) {
                lookup[nodeId] = pathways.map((pathway, index) => ({
                    pathIndex: index,
                    path: pathway.path,
                    pathLength: pathway.path_length || pathway.path.length - 1,
                    barrierHeight: pathway.barrier_height
                }));
            }
        }

        return lookup;
    }

    /**
     * Get reaction profile for a specific pathway (Phase 2.2 enhancement)
     * @param {number} nodeId - Target node ID
     * @param {number} pathIndex - Index of the pathway (default: 0 for best)
     * @returns {Array} Reaction profile data points
     */
    getReactionProfile(nodeId, pathIndex = 0) {
        if (!this.pathwayProfiles[nodeId] || !this.pathwayProfiles[nodeId][pathIndex]) {
            console.warn(`No reaction profile found for node ${nodeId}, path ${pathIndex}`);
            return [];
        }

        return this.pathwayProfiles[nodeId][pathIndex];
    }

    /**
     * Get pathway barrier height (Phase 2.2 enhancement)
     * @param {number} nodeId - Target node ID
     * @param {number} pathIndex - Index of the pathway (default: 0 for best)
     * @returns {number} Barrier height in kcal/mol
     */
    getPathwayBarrier(nodeId, pathIndex = 0) {
        if (!this.pathwayBarriers[nodeId] || this.pathwayBarriers[nodeId][pathIndex] === undefined) {
            console.warn(`No barrier data found for node ${nodeId}, path ${pathIndex}`);
            return Infinity;
        }

        return this.pathwayBarriers[nodeId][pathIndex];
    }

    /**
     * Get pathway length (Phase 2.2 enhancement)
     * @param {number} nodeId - Target node ID
     * @param {number} pathIndex - Index of the pathway (default: 0 for best)
     * @returns {number} Number of reaction steps
     */
    getPathwayLength(nodeId, pathIndex = 0) {
        const profile = this.getReactionProfile(nodeId, pathIndex);
        return profile.length > 0 ? profile.length - 1 : 0;
    }

    /**
     * Get all pathways to a node
     * @param {number} nodeId - Target node ID
     * @returns {Array} Array of pathway information
     */
    getAllPathways(nodeId) {
        return this.pathwayLookup[nodeId] || [];
    }

    /**
     * Check if pre-computed pathways are available
     * @returns {boolean} True if Python pathways are loaded
     */
    hasPrecomputedPathways() {
        return this.precomputedPathways !== null && Object.keys(this.precomputedPathways).length > 0;
    }

    /**
     * Main pathway analysis function - uses Python pathways when possible
     * @param {Object} networkData - The network data object
     * @param {Array} filteredEdges - Currently visible edges after filtering
     * @returns {Object} Best path analysis results
     */
    updateNetworkWithBestPaths(networkData, filteredEdges) {
        console.log("📋 updateNetworkWithBestPaths called with:");
        console.log("  - networkData.bestPathEnergyMap keys:", Object.keys(networkData.bestPathEnergyMap || {}).length);
        console.log("  - networkData.bestPathBarrierMap keys:", Object.keys(networkData.bestPathBarrierMap || {}).length);
        console.log("  - this.precomputedPathways keys:", Object.keys(this.precomputedPathways || {}).length);

        // Always use pre-computed pathways from Python when available
        // Only fall back to JavaScript BFS if explicitly filtering pathways or no Python data
        if (this.hasPrecomputedPathways()) {
            console.log("📊 Using pre-computed Python pathways (optimized path)");
            const results = this.buildResultsFromPythonPathways(networkData);

            // Debug logging
            const totalPaths = Object.keys(results.bestPathEnergyMap).length;
            const infinityPaths = Object.values(results.bestPathEnergyMap).filter(e => e === Infinity).length;
            const realPaths = Object.values(results.bestPathEnergyMap).filter(e => e !== Infinity).length;
            console.log(`✅ ${totalPaths} paths processed, ${infinityPaths} with Infinity values, ${realPaths} with real values`);

            // Sample a few values for debugging
            const sampleKeys = Object.keys(results.bestPathEnergyMap).slice(0, 5);
            for (const key of sampleKeys) {
                console.log(`  Sample node ${key}: energy=${results.bestPathEnergyMap[key]}, barrier=${results.bestPathBarrierMap[key]}`);
            }

            return results;
        } else {
            console.log("🔄 Computing pathways with JavaScript BFS (no Python data available)");
            return this.computePathwaysWithBFS(networkData, filteredEdges);
        }
    }

    /**
     * Check if filtered edges represent the full network (no filtering applied)
     * @param {Object} networkData - Network data
     * @param {Array} filteredEdges - Current edge set
     * @returns {boolean} True if using all edges
     */
    isUsingAllEdges(networkData, filteredEdges) {
        const originalEdgeCount = (networkData.originalEdgeData || []).length;
        const filteredEdgeCount = filteredEdges.length;

        // If counts match, likely using all edges (exact comparison could be added)
        return originalEdgeCount === filteredEdgeCount;
    }

    /**
     * Build results using pre-computed Python pathways
     * @param {Object} networkData - Network data with Python pathways
     * @returns {Object} Path analysis results
     */
    buildResultsFromPythonPathways(networkData) {
        console.log("🔧 Building results from Python pathways...");

        // Create empty maps
        const bestPathEnergyMap = {};
        const bestPathBarrierMap = {};
        const bestPathMap = {};

        // Use best path data computed by Python if available
        // Check for data in networkData (which is NetworkData instance)
        if (networkData.bestPathEnergyMap && Object.keys(networkData.bestPathEnergyMap).length > 0) {
            console.log("✅ Using pre-computed best path maps from Python");
            console.log(`📊 Found ${Object.keys(networkData.bestPathEnergyMap).length} energy entries`);
            console.log(`📊 Found ${Object.keys(networkData.bestPathBarrierMap || {}).length} barrier entries`);

            // Sample some values for verification
            const sampleKeys = Object.keys(networkData.bestPathEnergyMap).slice(0, 3);
            for (const key of sampleKeys) {
                console.log(`  Sample from networkData - Node ${key}: energy=${networkData.bestPathEnergyMap[key]}, barrier=${networkData.bestPathBarrierMap[key]}`);
            }

            Object.assign(bestPathEnergyMap, networkData.bestPathEnergyMap);
            Object.assign(bestPathBarrierMap, networkData.bestPathBarrierMap);
            Object.assign(bestPathMap, networkData.bestPathMap || {});
        } else {
            console.log("⚠️ NO pre-computed energy maps - using barrier_height only (data-only mode)");
            // Use ONLY precomputed barrier_height - NO energy calculations from profiles
            const pathwayCount = Object.keys(this.precomputedPathways).length;
            console.log(`📊 Processing ${pathwayCount} precomputed pathways for barriers only...`);

            for (const [nodeIdStr, pathways] of Object.entries(this.precomputedPathways)) {
                const nodeId = parseInt(nodeIdStr);

                if (pathways && pathways.length > 0) {
                    // Use the first pathway (best one, as sorted by Python)
                    const bestPathway = pathways[0];

                    // CRITICAL: DO NOT calculate energies from profiles - causes -200+ values
                    // JavaScript should never calculate energies - only use precomputed data
                    bestPathEnergyMap[nodeId] = null; // Indicates no precomputed energy available

                    // Use only precomputed barrier_height (this is safe)
                    bestPathBarrierMap[nodeId] = bestPathway.barrier_height || 0;
                    bestPathMap[nodeId] = bestPathway.path || [];
                } else {
                    bestPathEnergyMap[nodeId] = null;
                    bestPathBarrierMap[nodeId] = Infinity;
                    bestPathMap[nodeId] = [];
                }
            }
        }

        // Update networkData without overwriting original values
        if (!networkData.bestPathEnergyMap) networkData.bestPathEnergyMap = {};
        if (!networkData.bestPathBarrierMap) networkData.bestPathBarrierMap = {};
        if (!networkData.bestPathMap) networkData.bestPathMap = {};

        Object.assign(networkData.bestPathEnergyMap, bestPathEnergyMap);
        Object.assign(networkData.bestPathBarrierMap, bestPathBarrierMap);
        Object.assign(networkData.bestPathMap, bestPathMap);

        return {
            bestPathEnergyMap,
            bestPathBarrierMap,
            bestPathMap
        };
    }

    /**
     * Compute pathways using JavaScript BFS (fallback for filtering)
     * @param {Object} networkData - Network data
     * @param {Array} filteredEdges - Filtered edges
     * @returns {Object} Path analysis results
     */
    computePathwaysWithBFS(networkData, filteredEdges) {
        // Use the original BFS implementation
        const allShortestPaths = this.findAllShortestPaths(networkData, filteredEdges);
        const bestPathResults = this.findBestEnergyPaths(networkData, filteredEdges, allShortestPaths);

        // Convert Map results to plain objects for consistency
        const bestPathEnergyMap = {};
        const bestPathBarrierMap = {};
        const bestPathMap = {};

        for (const [nodeId, pathInfo] of bestPathResults.entries()) {
            bestPathEnergyMap[nodeId] = pathInfo.bestEnergy;
            bestPathBarrierMap[nodeId] = pathInfo.bestBarrier;
            bestPathMap[nodeId] = pathInfo.bestPath;
        }

        return {
            bestPathEnergyMap,
            bestPathBarrierMap,
            bestPathMap
        };
    }

    /**
     * Find all shortest paths from substrate to all nodes (by number of steps)
     * [LEGACY BFS IMPLEMENTATION - used when filtering is applied]
     * @param {Object} networkData - The network data object
     * @param {Array} filteredEdges - Currently visible edges after filtering
     * @returns {Object} Object of node_id to path information
     */
    findAllShortestPaths(networkData, filteredEdges) {
        const substrateId = networkData.substrateId;

        // Build adjacency graph from filtered edges
        const graph = new Map();

        filteredEdges.forEach(edge => {
            // Handle both edge formats: network edges (begin/end) and vis.js edges (from/to)
            const beginId = edge.begin || edge.from;
            const endId = edge.end || edge.to;

            if (!graph.has(beginId)) {
                graph.set(beginId, []);
            }
            // Store full edge data for later energy/barrier calculations (matching Python)
            graph.get(beginId).push({
                node: endId,
                edge_data: edge
            });
        });

        // Initialize results for all nodes
        const results = {};
        for (const node_id in networkData.relEnergyMap || {}) {
            const nodeId = parseInt(node_id);
            results[nodeId] = {
                distance: Infinity,
                all_paths: [],
                reachable: false
            };
        }

        // Substrate initialization (matching Python exactly)
        results[substrateId] = {
            distance: 0,
            all_paths: [[substrateId]],  // Python: [[substrate_id]]
            reachable: true
        };

        // BFS to find shortest paths (matching Python deque structure)
        const queue = [{ distance: 0, node: substrateId, path: [substrateId] }];
        const shortestDistances = new Map([[substrateId, 0]]);

        while (queue.length > 0) {
            const { distance: currentDist, node: currentNode, path: currentPath } = queue.shift();

            // Skip if we've moved beyond the shortest distance for this node
            if (shortestDistances.has(currentNode) && currentDist > shortestDistances.get(currentNode)) {
                continue;
            }

            // Explore neighbors
            const neighbors = graph.get(currentNode) || [];
            for (const neighborInfo of neighbors) {
                const neighborNode = neighborInfo.node;

                // Avoid cycles in individual paths (matching Python)
                if (currentPath.includes(neighborNode)) {
                    continue;
                }

                const newDist = currentDist + 1;
                const newPath = [...currentPath, neighborNode];  // Python: current_path + [neighbor_node]

                // Check if this is a shortest path
                if (!shortestDistances.has(neighborNode)) {
                    // First time reaching this node
                    shortestDistances.set(neighborNode, newDist);
                    results[neighborNode] = {
                        distance: newDist,
                        all_paths: [newPath],
                        reachable: true
                    };
                    queue.push({ distance: newDist, node: neighborNode, path: newPath });

                } else if (newDist === shortestDistances.get(neighborNode)) {
                    // Found another path of the same shortest length
                    const nodeResult = results[neighborNode];
                    if (!nodeResult.all_paths.some(path => this.arraysEqual(path, newPath))) {
                        nodeResult.all_paths.push(newPath);
                        queue.push({ distance: newDist, node: neighborNode, path: newPath });
                    }
                }
                // If new_dist > shortest_distances[neighbor_node], ignore (longer path)
            }
        }

        return results;
    }

    /**
     * Find the best energy path to each node based on lowest maximum energy barrier
     * @param {Object} networkData - The network data object
     * @param {Array} filteredEdges - Currently visible edges after filtering
     * @param {Object} allShortestPaths - Output from findAllShortestPaths()
     * @returns {Map} Map of node_id to best path information
     */
    findBestEnergyPaths(networkData, filteredEdges, allShortestPaths) {
        const substrateId = networkData.substrateId;

        // Build edge energy lookup using pre-calculated energies from Python
        const edgeEnergies = new Map();

        filteredEdges.forEach(edge => {
            // Handle both edge formats: network edges (begin/end) and vis.js edges (from/to)
            const beginId = edge.begin || edge.from;
            const endId = edge.end || edge.to;
            const key = `${beginId}-${endId}`;

            // Try to find reaction energy in various possible property names
            const reactionEnergy = this.extractReactionEnergy(edge);
            edgeEnergies.set(key, reactionEnergy);
        });

        const results = new Map();

        for (const [nodeIdStr, pathData] of Object.entries(allShortestPaths)) {
            const nodeId = parseInt(nodeIdStr);

            if (!pathData.reachable) {
                results.set(nodeId, {
                    best_path: [],
                    max_reaction_energy: Infinity,
                    cumulative_energy: Infinity,
                    edge_energies: [],
                    reachable: false,
                    final_node_energy: Infinity
                });
                continue;
            }

            let bestMaxEnergy = Infinity;
            let bestPath = null;
            let bestCumulative = 0.0;
            let bestEdgeEnergies = [];

            // Evaluate each shortest path
            for (const path of pathData.all_paths) {
                if (path.length === 1) { // Substrate itself (matching Python)
                    results.set(nodeId, {
                        best_path: path,
                        max_reaction_energy: 0.0,
                        cumulative_energy: 0.0,
                        edge_energies: [],
                        reachable: true,
                        final_node_energy: 0.0 // substrate is at energy 0
                    });
                    break;
                }

                // Calculate path energies using edge energies (matching Python)
                const pathEdgeEnergies = [];
                let pathValid = true;

                for (let i = 0; i < path.length - 1; i++) {
                    const fromNode = path[i];
                    const toNode = path[i + 1];
                    const edgeKey = `${fromNode}-${toNode}`;

                    // Get reaction energy for this step
                    const edgeEnergy = edgeEnergies.get(edgeKey);
                    if (edgeEnergy === undefined) {
                        pathValid = false;
                        break;
                    }
                    pathEdgeEnergies.push(edgeEnergy);
                }

                if (!pathValid) {
                    continue;
                }

                // Calculate cumulative energies along the path (matching Python)
                const cumulativeEnergies = [0.0]; // Start at substrate (energy 0)
                for (const edgeEnergy of pathEdgeEnergies) {
                    cumulativeEnergies.push(cumulativeEnergies[cumulativeEnergies.length - 1] + edgeEnergy);
                }

                // Calculate max uphill energy barrier with propagation (matching Python exactly)
                // This finds the maximum energy barrier that must be overcome along any path
                let maxUphillBarrier = 0.0;
                let runningMin = 0.0; // substrate energy (0)
                let runningMaxBarrier = 0.0; // Track the max barrier encountered so far

                for (let i = 1; i < cumulativeEnergies.length; i++) { // Skip substrate
                    const cumulativeEnergy = cumulativeEnergies[i];

                    // Update running minimum (lowest point seen so far)
                    runningMin = Math.min(runningMin, cumulativeEnergy);

                    // Calculate uphill barrier from lowest point to current point
                    const uphillBarrier = cumulativeEnergy - runningMin;

                    // The barrier for this node is the max of:
                    // 1. The highest barrier encountered previously
                    // 2. The current uphill barrier from the global minimum
                    runningMaxBarrier = Math.max(runningMaxBarrier, uphillBarrier);

                    // Update the overall max barrier for this path
                    maxUphillBarrier = Math.max(maxUphillBarrier, runningMaxBarrier);
                }

                // Final cumulative energy (relative to substrate)
                const finalCumulative = cumulativeEnergies[cumulativeEnergies.length - 1];

                // Here the relative energy calculation for neutral species should be inserted, it also needs to keep track of the eliminated products which is not gonna happen

                // Check if this path is better (lower max uphill barrier)
                if (maxUphillBarrier < bestMaxEnergy) {
                    bestMaxEnergy = maxUphillBarrier;
                    bestPath = path;
                    bestCumulative = finalCumulative;
                    bestEdgeEnergies = pathEdgeEnergies;
                }
            }

            // Handle case where no valid path was found
            if (bestPath === null) {
                results.set(nodeId, {
                    best_path: [],
                    max_reaction_energy: Infinity,
                    cumulative_energy: Infinity,
                    edge_energies: [],
                    reachable: false,
                    final_node_energy: Infinity
                });
                continue;
            }

            results.set(nodeId, {
                best_path: bestPath,
                max_reaction_energy: bestMaxEnergy, // Max uphill barrier with propagation
                cumulative_energy: bestCumulative,
                edge_energies: bestEdgeEnergies,
                reachable: true,
                final_node_energy: bestCumulative // Final energy relative to substrate
            });
        }

        return results;
    }

    /**
     * Extract reaction energy from edge object
     * @param {Object} edge - Edge object that may contain energy data in various properties
     * @returns {number} Reaction energy or Infinity if not found
     */
    extractReactionEnergy(edge) {
        // Try different property names where Python might have stored the calculated energy
        if (edge.rxn_energy !== undefined) {
            return edge.rxn_energy;
        } else if (edge.reaction_energy !== undefined) {
            return edge.reaction_energy;
        } else if (edge.energy !== undefined) {
            return edge.energy;
        } else if (edge.barrier !== undefined) {
            return edge.barrier;
        } else if (edge.originalData && edge.originalData.rxn_energy !== undefined) {
            return edge.originalData.rxn_energy;
        } else if (edge.originalData && edge.originalData.reaction_energy !== undefined) {
            return edge.originalData.reaction_energy;
        } else if (edge.originalData && edge.originalData.energy !== undefined) {
            return edge.originalData.energy;
        } else if (edge.originalData && edge.originalData.barrier !== undefined) {
            return edge.originalData.barrier;
        }

        // Default to Infinity if no energy found (matching Python behavior)
        return Infinity;
    }

    /**
     * Helper function to compare arrays
     */
    arraysEqual(a, b) {
        if (a.length !== b.length) return false;
        for (let i = 0; i < a.length; i++) {
            if (a[i] !== b[i]) return false;
        }
        return true;
    }

    /**
     * Get analysis results for debugging
     */
    getResults() {
        return this.results;
    }

    // ===== PHASE 2.3: TRANSITION STATE SUPPORT =====

    /**
     * Get transition state information for a specific reaction step
     * @param {number} nodeId - Target node ID
     * @param {number} pathIndex - Pathway index
     * @param {number} stepIndex - Step within the pathway
     * @returns {Object|null} Transition state data or null if not available
     */
    getTransitionStateInfo(nodeId, pathIndex = 0, stepIndex = 0) {
        const profile = this.getReactionProfile(nodeId, pathIndex);

        if (!profile || stepIndex >= profile.length - 1) {
            return null;
        }

        // Look for transition state between stepIndex and stepIndex + 1
        const currentStep = profile[stepIndex];
        const nextStep = profile[stepIndex + 1];

        // Check if there's transition state data in the network
        const networkData = window.networkMain?.networkData;
        if (!networkData?.originalEdgeData) {
            return null;
        }

        // Find the edge corresponding to this step
        const edge = networkData.originalEdgeData.find(e =>
            e.begin === currentStep.nodeId && e.end === nextStep.nodeId
        );

        if (edge?.ts) {
            return {
                molBlock: edge.ts,
                energy: edge.ts_energy,
                barrier: edge.barrier || 0,
                stepIndex: stepIndex,
                beginNodeId: currentStep.nodeId,
                endNodeId: nextStep.nodeId
            };
        }

        return null;
    }

    /**
     * Check if a pathway has transition states
     * @param {number} nodeId - Target node ID
     * @param {number} pathIndex - Pathway index
     * @returns {boolean} True if pathway contains transition states
     */
    hasTransitionStates(nodeId, pathIndex = 0) {
        const profile = this.getReactionProfile(nodeId, pathIndex);

        if (!profile || profile.length < 2) {
            return false;
        }

        // Check each step for transition state data
        for (let i = 0; i < profile.length - 1; i++) {
            if (this.getTransitionStateInfo(nodeId, pathIndex, i)) {
                return true;
            }
        }

        return false;
    }

    /**
     * Get all transition states in a pathway
     * @param {number} nodeId - Target node ID
     * @param {number} pathIndex - Pathway index
     * @returns {Array} Array of transition state information
     */
    getAllTransitionStates(nodeId, pathIndex = 0) {
        const profile = this.getReactionProfile(nodeId, pathIndex);
        const transitionStates = [];

        if (!profile || profile.length < 2) {
            return transitionStates;
        }

        for (let i = 0; i < profile.length - 1; i++) {
            const tsInfo = this.getTransitionStateInfo(nodeId, pathIndex, i);
            if (tsInfo) {
                transitionStates.push(tsInfo);
            }
        }

        return transitionStates;
    }

    /**
     * Enhanced reaction profile with transition state markers
     * @param {number} nodeId - Target node ID
     * @param {number} pathIndex - Pathway index
     * @returns {Array} Enhanced profile with TS markers
     */
    getEnhancedReactionProfile(nodeId, pathIndex = 0) {
        const profile = this.getReactionProfile(nodeId, pathIndex);

        if (!profile) {
            return [];
        }

        const enhancedProfile = [];

        for (let i = 0; i < profile.length; i++) {
            const point = { ...profile[i] };

            // Add minimum marker
            point.type = 'minimum';
            enhancedProfile.push(point);

            // Add transition state if it exists before next step
            if (i < profile.length - 1) {
                const tsInfo = this.getTransitionStateInfo(nodeId, pathIndex, i);
                if (tsInfo) {
                    enhancedProfile.push({
                        nodeId: `TS${i + 1}`,
                        energy: tsInfo.energy || (point.energy + tsInfo.barrier),
                        step: i + 0.5,
                        reactionCoordinate: point.reactionCoordinate + 0.5,
                        type: 'transition_state',
                        transitionStateInfo: tsInfo
                    });
                }
            }
        }

        return enhancedProfile;
    }
}

// Export for use in other modules
window.PathAnalyzer = PathAnalyzer;

// Initialize global path analyzer instance
window.pathAnalyzer = new PathAnalyzer();
