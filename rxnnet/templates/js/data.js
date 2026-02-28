/**
 * Data management for the network visualization
 */
const Data = {
    raw: null,
    nodes: [],
    edges: [],
    junctionNodes: [],
    filteredNodes: [],
    filteredEdges: [],
    rxnTypes: new Set(),

    init(data) {
        this.raw = data;
        this.parseData();
    },

    parseData() {
        const { svgMap, molEnergyMap, chargeMap, weightMap, relEnergyMap,
                originalEdgeData, substrateId, bestPathBarrierMap, bestPathEnergyMap } = this.raw;

        // Parse nodes
        this.nodes = Object.keys(svgMap).map(id => {
            const nodeId = parseInt(id);
            return {
                id: nodeId,
                svg: svgMap[id],
                energy: molEnergyMap[id],
                charge: chargeMap[id] || 0,
                weight: weightMap[id] || 0,
                relEnergy: relEnergyMap[id] || 0,
                barrier: bestPathBarrierMap[nodeId] !== undefined ? bestPathBarrierMap[nodeId] : null,
                pathEnergy: bestPathEnergyMap[nodeId] !== undefined ? bestPathEnergyMap[nodeId] : null,
                isSubstrate: nodeId === substrateId
            };
        });

        // Parse edges and create hyperedge junctions
        let edgeId = 0;
        let junctionId = -1;
        this.edges = [];
        this.junctionNodes = [];

        originalEdgeData.forEach(e => {
            this.rxnTypes.add(e.type);
            const eliminatedProducts = e.smaller_products || [];

            if (eliminatedProducts.length > 0) {
                // Hyperedge: create junction node for branching
                const jId = junctionId--;
                this.junctionNodes.push({
                    id: jId,
                    isJunction: true,
                    reactantId: e.begin,
                    productIds: [e.end, ...eliminatedProducts]
                });

                // Edge from reactant to junction
                this.edges.push({
                    id: edgeId++,
                    from: e.begin,
                    to: jId,
                    type: e.type,
                    count: e.count,
                    isToJunction: true
                });

                // Edge from junction to primary product
                this.edges.push({
                    id: edgeId++,
                    from: jId,
                    to: e.end,
                    type: e.type,
                    count: e.count,
                    isFromJunction: true,
                    isPrimary: true
                });

                // Edges from junction to eliminated products
                eliminatedProducts.forEach(elimId => {
                    this.edges.push({
                        id: edgeId++,
                        from: jId,
                        to: elimId,
                        type: e.type,
                        count: e.count,
                        isFromJunction: true,
                        isEliminated: true
                    });
                });
            } else {
                // Simple edge: direct connection
                this.edges.push({
                    id: edgeId++,
                    from: e.begin,
                    to: e.end,
                    type: e.type,
                    count: e.count
                });
            }
        });

        this.filteredNodes = [...this.nodes];
        this.filteredEdges = [...this.edges];
    },

    getNode(id) {
        return this.nodes.find(n => n.id === id);
    },

    getPathways(nodeId) {
        return this.raw.pathways[nodeId] || [];
    },

    getBestPath(nodeId) {
        return this.raw.bestPathMap[nodeId] || null;
    },

    filter({ maxBarrier = Infinity, maxEnergy = Infinity, minCount = 1 } = {}) {
        // Filter nodes by barrier and energy
        const visibleNodeIds = new Set();
        this.nodes.forEach(n => {
            const passBarrier = n.isSubstrate || n.barrier === null || n.barrier <= maxBarrier;
            const passEnergy = n.isSubstrate || n.pathEnergy === null || n.pathEnergy <= maxEnergy;
            if (passBarrier && passEnergy) {
                visibleNodeIds.add(n.id);
            }
        });

        // Filter edges by count and connected nodes
        this.filteredEdges = this.edges.filter(e => {
            const countPass = (e.count || 1) >= minCount;
            // For junction edges, check if junction is valid
            if (e.isToJunction || e.isFromJunction) {
                return countPass; // Junction edges handled separately
            }
            const nodesExist = visibleNodeIds.has(e.from) && visibleNodeIds.has(e.to);
            return countPass && nodesExist;
        });

        // Get connected nodes from filtered edges
        const connectedIds = new Set();
        this.filteredEdges.forEach(e => {
            if (e.from > 0) connectedIds.add(e.from);
            if (e.to > 0) connectedIds.add(e.to);
        });

        // Include substrate
        if (this.raw.substrateId) {
            connectedIds.add(this.raw.substrateId);
        }

        this.filteredNodes = this.nodes.filter(n => connectedIds.has(n.id));

        // Get visible junctions
        const visibleJunctions = this.junctionNodes.filter(j =>
            this.filteredEdges.some(e => e.from === j.id || e.to === j.id)
        );

        return {
            nodes: this.filteredNodes,
            edges: this.filteredEdges,
            junctions: visibleJunctions
        };
    },

    reset() {
        this.filteredNodes = [...this.nodes];
        this.filteredEdges = [...this.edges];
        return {
            nodes: this.filteredNodes,
            edges: this.filteredEdges,
            junctions: this.junctionNodes
        };
    }
};
