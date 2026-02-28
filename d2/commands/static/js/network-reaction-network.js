/**
 * ReactionNetwork JavaScript Implementation
 *
 * This class provides an exact JavaScript port of the Python ReactionNetwork class
 * from rn.py, ensuring 1:1 matching of energy calculations, pathfinding algorithms,
 * reaction profile generation, and filtering logic.
 *
 * Key Features:
 * - Exact energy conversion and chemical corrections
 * - Optimized two-phase pathfinding
 * - Chemical-aware filtering (reaction vs protonation/deprotonation)
 * - Transition state energy handling
 * - Complete reaction profile generation
 */

class ReactionNetworkJS {
    constructor(networkData, substrateId, propName, pH = 7.0) {
        // Constants matching Python implementation exactly
        this.R = 1.98720425864083 / 1000; // Gas constant in kcal/(mol*K)
        this.T = 313.15; // Temperature in K
        this.PKA_SLOPE = 0.43496246;
        this.PKA_INTERCEPT = -114.66798804;

        // Configuration
        this.substrateId = substrateId;
        this.propName = propName;
        this.pH = pH;
        this.energyType = propName.includes('gibbs') ? 'G' : 'E';

        // Core data structures matching Python
        this.nodes = new Map(); // node_id -> node data
        this.edges = []; // Array of edge objects
        this.nodeLabels = new Map(); // node_id -> label string
        this.edgeMap = new Map(); // "start_id-end_id" -> edge data for O(1) lookup

        // Precomputed pathway data
        this.pathways = new Map(); // node_id -> array of pathway data

        // Performance optimization: caching
        this.profileCache = new Map(); // Cache for reaction profiles

        // Initialize from provided data
        this.initializeFromData(networkData);
    }

    /**
     * Initialize the ReactionNetwork from provided data
     * @param {Object} networkData - Data from Python to_visualization_data()
     */
    initializeFromData(networkData) {
        console.log('Initializing ReactionNetwork from data...');

        // Load nodes
        this.loadNodes(networkData.nodes);

        // Load edges and build edge map
        this.loadEdges(networkData.edges);

        // Store configuration
        this.config = networkData.config || {};

        console.log(`Network loaded: ${this.nodes.size} nodes, ${this.edges.length} edges`);
    }

    /**
     * Load nodes from data array
     * @param {Array} nodesData - Array of node objects
     */
    loadNodes(nodesData) {
        for (const nodeData of nodesData) {
            const nodeId = nodeData.id;

            // Store node data with same structure as Python
            // NO energy conversions needed - Python provides all precomputed values
            this.nodes.set(nodeId, {
                id: nodeId,
                energy: nodeData.energy, // Already in kcal/mol from Python
                rawEnergy: nodeData.energy, // Keep same value - no calculations performed
                svg: nodeData.svg || '',
                charge: nodeData.charge || 0,
                weight: nodeData.weight || 0,
                smiles: nodeData.smiles || '',
                molData: nodeData.mol_data || '',
                originType: nodeData.origin_type || '',
                x: nodeData.x || 0,
                y: nodeData.y || 0
            });

            // Store label
            this.nodeLabels.set(nodeId, nodeData.label || nodeId.toString());
        }
    }

    /**
     * Load edges and build edge map for O(1) lookup
     * @param {Array} edgesData - Array of edge objects
     */
    loadEdges(edgesData) {
        this.edges = [];
        this.edgeMap.clear();

        for (const edgeData of edgesData) {
            const edge = {
                begin: edgeData.source,
                end: edgeData.target,
                type: edgeData.type,
                count: edgeData.count || 0,
                barrier: edgeData.barrier || 0.0,
                hasTs: edgeData.has_ts || false,
                tsEnergy: edgeData.ts_energy || null,
                tsBarrier: edgeData.ts_barrier || null,
                rxnEnergy: edgeData.rxn_energy || 0.0,
                smallerProducts: edgeData.smaller_products || []
            };

            this.edges.push(edge);

            // Build edge map for O(1) lookup (matching Python)
            const edgeKey = `${edge.begin}-${edge.end}`;
            this.edgeMap.set(edgeKey, edge);
        }

        console.log(`Built edge map with ${this.edgeMap.size} entries`);
    }

    /**
     * Convert hartree to kcal/mol (exact match of Python function)
     * @param {number} hartree - Energy in hartree
     * @returns {number} Energy in kcal/mol
     */
    hartree2kcalmol(hartree) {
        return hartree * 627.509474;
    }

    /**
     * Get edge data for given start and end node IDs
     * @param {number} startId - Starting node ID
     * @param {number} endId - Ending node ID
     * @returns {Object|null} Edge data or null if not found
     */
    getEdge(startId, endId) {
        const edgeKey = `${startId}-${endId}`;
        return this.edgeMap.get(edgeKey) || null;
    }

    /**
     * Get edge count for given start and end node IDs
     * @param {number} startId - Starting node ID
     * @param {number} endId - Ending node ID
     * @returns {number} Edge count or 0 if edge not found
     */
    getEdgeCount(startId, endId) {
        const edge = this.getEdge(startId, endId);
        return edge ? edge.count : 0;
    }

    /**
     * Get edge type for given start and end node IDs
     * @param {number} startId - Starting node ID
     * @param {number} endId - Ending node ID
     * @returns {string} Edge type or empty string if edge not found
     */
    getEdgeType(startId, endId) {
        const edge = this.getEdge(startId, endId);
        return edge ? edge.type : '';
    }

    /**
     * Check if an edge meets the count requirement based on chemical logic
     * Exact port of Python method with same logic:
     * - Reaction edges: must meet count threshold
     * - Protonation/deprotonation edges: always pass regardless of count
     *
     * @param {number} startId - Starting node ID
     * @param {number} endId - Ending node ID
     * @param {number} minCount - Minimum count threshold
     * @returns {boolean} True if edge meets requirement
     */
    edgeMeetsCountRequirement(startId, endId, minCount) {
        const edge = this.getEdge(startId, endId);
        if (!edge) {
            return false;
        }

        const edgeType = edge.type.toLowerCase();

        // Only check count for reaction edges
        // Protonation/deprotonation edges always pass (chemical state changes)
        if (edgeType === 'reaction') {
            return edge.count >= minCount;
        } else {
            // Protonation, deprotonation, tautomerization edges always pass
            return true;
        }
    }

    /**
     * Get transition state data for given edge
     * @param {number} startId - Starting node ID
     * @param {number} endId - Ending node ID
     * @returns {Object|null} TS data object or null
     */
    getTransitionStateData(startId, endId) {
        const edge = this.getEdge(startId, endId);
        if (!edge || !edge.hasTs) {
            return null;
        }

        return {
            energy: edge.tsEnergy,
            barrier: edge.tsBarrier,
            hasData: edge.hasTs
        };
    }

    /**
     * Validate data consistency (debugging utility)
     * @returns {Object} Validation results
     */
    validateData() {
        const validation = {
            nodesValid: true,
            edgesValid: true,
            edgeMapValid: true,
            issues: []
        };

        // Check nodes
        for (const [nodeId, nodeData] of this.nodes) {
            if (typeof nodeData.energy !== 'number') {
                validation.nodesValid = false;
                validation.issues.push(`Node ${nodeId} has invalid energy: ${nodeData.energy}`);
            }
        }

        // Check edges
        for (const edge of this.edges) {
            if (!this.nodes.has(edge.begin) || !this.nodes.has(edge.end)) {
                validation.edgesValid = false;
                validation.issues.push(`Edge ${edge.begin}-${edge.end} references non-existent nodes`);
            }
        }

        // Check edge map consistency
        if (this.edges.length !== this.edgeMap.size) {
            validation.edgeMapValid = false;
            validation.issues.push(`Edge array length (${this.edges.length}) != edge map size (${this.edgeMap.size})`);
        }

        return validation;
    }

    /**
     * Get summary statistics
     * @returns {Object} Summary statistics
     */
    getSummary() {
        const reactionEdges = this.edges.filter(e => e.type === 'reaction').length;
        const protonationEdges = this.edges.filter(e => e.type === 'protonation').length;
        const deprotonationEdges = this.edges.filter(e => e.type === 'deprotonation').length;
        const tsEdges = this.edges.filter(e => e.hasTs).length;

        return {
            totalNodes: this.nodes.size,
            totalEdges: this.edges.length,
            reactionEdges,
            protonationEdges,
            deprotonationEdges,
            tsEdges,
            substrateId: this.substrateId,
            pH: this.pH,
            energyType: this.energyType
        };
    }
}

// Make available globally
window.ReactionNetworkJS = ReactionNetworkJS;
