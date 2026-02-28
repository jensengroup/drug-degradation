
/*
 * Network Pathway Optimization Module
 * Implements optimal pathway finding algorithms for chemical reaction networks
 */

class PathwayOptimizationMixin {
    /**
     * Find the optimal (lowest energy barrier) pathway to a target node
     * @param {number} targetNodeId - The target node to reach
     * @param {Object} options - Optimization options
     * @returns {Object} Pathway with energy profile and metadata
     */
    findOptimalPathway(targetNodeId, options = {}) {
        const {
            maxPathLength = 8,
            energyThreshold = 100,  // kcal/mol
            algorithm = 'dijkstra',
            considerKinetics = true,
            allowCycles = false
        } = options;

        if (!this.hasNode(targetNodeId)) {
            throw new Error(`Target node ${targetNodeId} not found in network`);
        }

        console.log(`🎯 Finding optimal pathway to node ${targetNodeId}`);
        console.log(`   Algorithm: ${algorithm}, Max length: ${maxPathLength}`);

        switch (algorithm) {
            case 'dijkstra':
                return this._dijkstraOptimalPath(targetNodeId, options);
            case 'a_star':
                return this._aStarOptimalPath(targetNodeId, options);
            case 'dynamic_programming':
                return this._dynamicProgrammingPath(targetNodeId, options);
            default:
                throw new Error(`Unknown optimization algorithm: ${algorithm}`);
        }
    }

    /**
     * Dijkstra's algorithm for finding shortest energy path
     */
    _dijkstraOptimalPath(targetNodeId, options) {
        const { maxPathLength, energyThreshold, considerKinetics } = options;

        // Priority queue: [distance, nodeId, path, profile]
        const pq = new MinHeap((a, b) => a[0] - b[0]);
        const distances = new Map();
        const visited = new Set();

        // Initialize with substrate
        pq.push([0, this.substrateId, [this.substrateId], [{
            node: this.substrateId,
            energy: 0,
            type: 'minimum',
            reaction_coordinate: 0.0
        }]]);
        distances.set(this.substrateId, 0);

        let iterations = 0;
        const maxIterations = 10000;

        while (!pq.isEmpty() && iterations < maxIterations) {
            iterations++;
            const [currentDistance, currentNode, currentPath, currentProfile] = pq.pop();

            // Skip if already visited with better distance
            if (visited.has(currentNode)) continue;
            visited.add(currentNode);

            // Found target
            if (currentNode === targetNodeId) {
                const barrier = this.findMaxBarrier(currentProfile);
                console.log(`✅ Found optimal pathway in ${iterations} iterations`);
                console.log(`   Path: ${currentPath.join(' → ')}`);
                console.log(`   Barrier: ${barrier.toFixed(2)} kcal/mol`);

                return {
                    path: currentPath,
                    profile: currentProfile,
                    barrier_height: barrier,
                    algorithm: 'dijkstra',
                    iterations: iterations,
                    metadata: {
                        path_length: currentPath.length,
                        energy_threshold: energyThreshold,
                        consider_kinetics: considerKinetics
                    }
                };
            }

            // Skip if path too long
            if (currentPath.length >= maxPathLength) continue;

            // Skip if current barrier already too high
            if (currentDistance > energyThreshold) continue;

            // Explore neighbors
            const neighbors = this.getNeighbors(currentNode);
            for (const neighborId of neighbors) {
                if (visited.has(neighborId)) continue;

                // Calculate edge energy (transition state barrier)
                const edgeEnergy = this._calculateEdgeEnergy(currentNode, neighborId, considerKinetics);
                if (edgeEnergy === null) continue;

                const newDistance = Math.max(currentDistance, edgeEnergy);

                // Skip if this is a worse path to this neighbor
                if (distances.has(neighborId) && newDistance >= distances.get(neighborId)) {
                    continue;
                }

                distances.set(neighborId, newDistance);

                // Create new path and profile
                const newPath = [...currentPath, neighborId];
                const newProfile = this._extendProfile(currentProfile, currentNode, neighborId, edgeEnergy);

                pq.push([newDistance, neighborId, newPath, newProfile]);
            }
        }

        console.log(`❌ No pathway found to node ${targetNodeId} within constraints`);
        console.log(`   Iterations: ${iterations}, Visited: ${visited.size} nodes`);

        return null;
    }

    /**
     * A* algorithm with energy-based heuristic
     */
    _aStarOptimalPath(targetNodeId, options) {
        const { maxPathLength, energyThreshold } = options;

        // Priority queue: [f_score, g_score, nodeId, path, profile]
        const pq = new MinHeap((a, b) => a[0] - b[0]);
        const gScore = new Map();
        const visited = new Set();

        // Heuristic: estimate minimum energy to target
        const heuristic = (nodeId) => {
            // Simple heuristic: use chemical similarity or direct energy if available
            return 0; // Conservative estimate
        };

        // Initialize
        const startG = 0;
        const startF = startG + heuristic(this.substrateId);
        pq.push([startF, startG, this.substrateId, [this.substrateId], [{
            node: this.substrateId,
            energy: 0,
            type: 'minimum',
            reaction_coordinate: 0.0
        }]]);
        gScore.set(this.substrateId, startG);

        let iterations = 0;
        const maxIterations = 10000;

        while (!pq.isEmpty() && iterations < maxIterations) {
            iterations++;
            const [fScore, gScore_current, currentNode, currentPath, currentProfile] = pq.pop();

            if (visited.has(currentNode)) continue;
            visited.add(currentNode);

            // Found target
            if (currentNode === targetNodeId) {
                const barrier = this.findMaxBarrier(currentProfile);
                console.log(`✅ A* found optimal pathway in ${iterations} iterations`);

                return {
                    path: currentPath,
                    profile: currentProfile,
                    barrier_height: barrier,
                    algorithm: 'a_star',
                    iterations: iterations,
                    metadata: {
                        path_length: currentPath.length,
                        f_score: fScore,
                        g_score: gScore_current
                    }
                };
            }

            // Explore neighbors
            const neighbors = this.getNeighbors(currentNode);
            for (const neighborId of neighbors) {
                if (visited.has(neighborId)) continue;
                if (currentPath.length >= maxPathLength) continue;

                const edgeEnergy = this._calculateEdgeEnergy(currentNode, neighborId, true);
                if (edgeEnergy === null) continue;

                const tentativeG = Math.max(gScore_current, edgeEnergy);

                if (gScore.has(neighborId) && tentativeG >= gScore.get(neighborId)) {
                    continue;
                }

                gScore.set(neighborId, tentativeG);
                const f = tentativeG + heuristic(neighborId);

                const newPath = [...currentPath, neighborId];
                const newProfile = this._extendProfile(currentProfile, currentNode, neighborId, edgeEnergy);

                pq.push([f, tentativeG, neighborId, newPath, newProfile]);
            }
        }

        return null;
    }

    /**
     * Calculate the energy barrier for traversing an edge
     */
    _calculateEdgeEnergy(fromNode, toNode, considerKinetics = true) {
        const edge = this.findEdge(fromNode, toNode);
        if (!edge) return null;

        try {
            // Use precomputed barrier/energy data if available for the target node
            if (considerKinetics) {
                // Return barrier from precomputed data
                const barrier = this.getBestPathwayBarrier(toNode);
                return barrier !== Infinity ? barrier : null;
            } else {
                // Return energy difference from precomputed data
                const energy = this.getBestPathwayEnergy(toNode);
                return energy !== Infinity ? Math.abs(energy) : null;
            }
        } catch (error) {
            console.warn(`Error getting edge energy ${fromNode}→${toNode}:`, error);
            return null;
        }
    }

    /**
     * Extend a reaction profile with a new edge (using precomputed data)
     */
    _extendProfile(currentProfile, fromNode, toNode, edgeEnergy) {
        try {
            // Get precomputed profile for the target node
            const targetProfile = this.getPrecomputedProfile(toNode, 0);

            if (!targetProfile || targetProfile.length === 0) {
                console.warn(`Could not get precomputed profile for node ${toNode}`);
                return currentProfile;
            }

            // Merge profiles (skip the first point of target profile to avoid duplication)
            const lastCoord = currentProfile[currentProfile.length - 1].reactionCoordinate || 0;
            const mergedProfile = [...currentProfile];

            for (let i = 1; i < targetProfile.length; i++) {
                const point = { ...targetProfile[i] };
                point.reactionCoordinate = lastCoord + i;
                mergedProfile.push(point);
            }

            return mergedProfile;
        } catch (error) {
            console.warn(`Error extending profile with ${fromNode}→${toNode}:`, error);
            return currentProfile;
        }
    }

    /**
     * Get neighboring nodes from current node
     */
    getNeighbors(nodeId) {
        const neighbors = [];
        for (const edge of this.edges) {
            if (edge.begin === nodeId) {
                neighbors.push(edge.end);
            } else if (edge.end === nodeId) {
                neighbors.push(edge.begin);
            }
        }
        return neighbors;
    }

    /**
     * Check if a node exists in the network
     */
    hasNode(nodeId) {
        return this.nodes.some(node => node.id === nodeId);
    }

    /**
     * Find edge between two nodes
     */
    findEdge(fromNode, toNode) {
        return this.edges.find(edge =>
            (edge.begin === fromNode && edge.end === toNode) ||
            (edge.begin === toNode && edge.end === fromNode)
        );
    }

    /**
     * Find multiple optimal pathways to target
     */
    findMultipleOptimalPathways(targetNodeId, options = {}) {
        const {
            maxPathways = 5,
            energyTolerance = 5.0,  // kcal/mol tolerance from best pathway
            diversityThreshold = 0.5  // Jaccard similarity threshold
        } = options;

        console.log(`🎯 Finding up to ${maxPathways} optimal pathways to node ${targetNodeId}`);

        const pathways = [];
        const usedNodes = new Set();

        // Find first optimal pathway
        const firstPath = this.findOptimalPathway(targetNodeId, options);
        if (!firstPath) return [];

        pathways.push(firstPath);
        const bestBarrier = firstPath.barrier_height;

        // Track nodes used in first path for diversity
        firstPath.path.forEach(node => usedNodes.add(node));

        // Find additional diverse pathways
        for (let i = 1; i < maxPathways; i++) {
            // Modify options to encourage diversity
            const diverseOptions = {
                ...options,
                energyThreshold: bestBarrier + energyTolerance,
                penalizeUsedNodes: true,
                usedNodes: new Set(usedNodes)
            };

            const pathway = this._findDiversePathway(targetNodeId, diverseOptions, pathways);
            if (pathway && pathway.barrier_height <= bestBarrier + energyTolerance) {
                // Check diversity
                const diversity = this._calculatePathwayDiversity(pathway, pathways);
                if (diversity >= diversityThreshold) {
                    pathways.push(pathway);
                    pathway.path.forEach(node => usedNodes.add(node));
                }
            }
        }

        // Sort by barrier height
        pathways.sort((a, b) => a.barrier_height - b.barrier_height);

        console.log(`✅ Found ${pathways.length} diverse optimal pathways`);
        pathways.forEach((pathway, i) => {
            console.log(`   ${i + 1}. ${pathway.path.join('→')} (${pathway.barrier_height.toFixed(2)} kcal/mol)`);
        });

        return pathways;
    }

    /**
     * Calculate diversity between pathways using Jaccard similarity
     */
    _calculatePathwayDiversity(newPathway, existingPathways) {
        const newNodes = new Set(newPathway.path);

        let maxSimilarity = 0;
        for (const existing of existingPathways) {
            const existingNodes = new Set(existing.path);
            const intersection = new Set([...newNodes].filter(x => existingNodes.has(x)));
            const union = new Set([...newNodes, ...existingNodes]);
            const similarity = intersection.size / union.size;
            maxSimilarity = Math.max(maxSimilarity, similarity);
        }

        return 1 - maxSimilarity; // Convert similarity to diversity
    }

    /**
     * Validate pathway optimization capabilities
     */
    validatePathwayOptimization() {
        const validation = {
            methods: {
                findOptimalPathway: typeof this.findOptimalPathway === 'function',
                findMultipleOptimalPathways: typeof this.findMultipleOptimalPathways === 'function',
                dijkstra: typeof this._dijkstraOptimalPath === 'function',
                astar: typeof this._aStarOptimalPath === 'function'
            },
            utilities: {
                calculateEdgeEnergy: typeof this._calculateEdgeEnergy === 'function',
                extendProfile: typeof this._extendProfile === 'function',
                getNeighbors: typeof this.getNeighbors === 'function',
                hasNode: typeof this.hasNode === 'function',
                findEdge: typeof this.findEdge === 'function'
            },
            data_structures: {
                nodes: Array.isArray(this.nodes),
                edges: Array.isArray(this.edges),
                edgeMap: this.edgeMap instanceof Map
            }
        };

        const allValid = Object.values(validation.methods).every(v => v) &&
            Object.values(validation.utilities).every(v => v) &&
            Object.values(validation.data_structures).every(v => v);

        validation.overall = allValid;
        return validation;
    }
}

/**
 * Min-heap implementation for priority queue
 */
class MinHeap {
    constructor(compareFunction = (a, b) => a - b) {
        this.heap = [];
        this.compare = compareFunction;
    }

    push(element) {
        this.heap.push(element);
        this._heapifyUp(this.heap.length - 1);
    }

    pop() {
        if (this.heap.length === 0) return undefined;
        if (this.heap.length === 1) return this.heap.pop();

        const root = this.heap[0];
        this.heap[0] = this.heap.pop();
        this._heapifyDown(0);
        return root;
    }

    isEmpty() {
        return this.heap.length === 0;
    }

    _heapifyUp(index) {
        if (index === 0) return;

        const parentIndex = Math.floor((index - 1) / 2);
        if (this.compare(this.heap[index], this.heap[parentIndex]) < 0) {
            [this.heap[index], this.heap[parentIndex]] = [this.heap[parentIndex], this.heap[index]];
            this._heapifyUp(parentIndex);
        }
    }

    _heapifyDown(index) {
        const leftChild = 2 * index + 1;
        const rightChild = 2 * index + 2;
        let smallest = index;

        if (leftChild < this.heap.length &&
            this.compare(this.heap[leftChild], this.heap[smallest]) < 0) {
            smallest = leftChild;
        }

        if (rightChild < this.heap.length &&
            this.compare(this.heap[rightChild], this.heap[smallest]) < 0) {
            smallest = rightChild;
        }

        if (smallest !== index) {
            [this.heap[index], this.heap[smallest]] = [this.heap[smallest], this.heap[index]];
            this._heapifyDown(smallest);
        }
    }
}

// Apply mixin to ReactionNetworkJS
Object.assign(ReactionNetworkJS.prototype, PathwayOptimizationMixin.prototype);
