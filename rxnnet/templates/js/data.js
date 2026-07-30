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

        // Eliminated (byproduct) nodes are produced in the same reaction step
        // as their hyperedge's primary product, so they never get their own
        // entry in bestPathBarrierMap/bestPathEnergyMap (the backend only
        // traces pathways through primary products). Propagate the primary
        // product's values so byproducts get pruned alongside it instead of
        // always passing the filter.
        const nodesById = new Map(this.nodes.map(n => [n.id, n]));
        originalEdgeData.forEach(e => {
            const eliminatedProducts = e.smaller_products || [];
            if (eliminatedProducts.length === 0) return;
            const primaryNode = nodesById.get(e.end);
            if (!primaryNode) return;

            eliminatedProducts.forEach(elimId => {
                const elimNode = nodesById.get(elimId);
                if (!elimNode) return;
                if (elimNode.barrier === null ||
                    (primaryNode.barrier !== null && primaryNode.barrier < elimNode.barrier)) {
                    elimNode.barrier = primaryNode.barrier;
                }
                if (elimNode.pathEnergy === null ||
                    (primaryNode.pathEnergy !== null && primaryNode.pathEnergy < elimNode.pathEnergy)) {
                    elimNode.pathEnergy = primaryNode.pathEnergy;
                }
            });
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

    filter({ maxPathEnergy = Infinity, minCount = 1 } = {}) {
        // Filter nodes by best-pathway barrier (max path energy)
        const visibleNodeIds = new Set();
        this.nodes.forEach(n => {
            const passPathEnergy = n.isSubstrate || n.barrier === null || n.barrier <= maxPathEnergy;
            if (passPathEnergy) {
                visibleNodeIds.add(n.id);
            }
        });

        // Filter edges by count and connected nodes
        this.filteredEdges = this.edges.filter(e => {
            const countPass = (e.count || 1) >= minCount;
            // Junction (hyperedge) legs only have one real node endpoint;
            // the junction id itself is virtual and not in visibleNodeIds.
            if (e.isToJunction) {
                return countPass && visibleNodeIds.has(e.from);
            }
            if (e.isFromJunction) {
                return countPass && visibleNodeIds.has(e.to);
            }
            const nodesExist = visibleNodeIds.has(e.from) && visibleNodeIds.has(e.to);
            return countPass && nodesExist;
        });

        // A junction (the little grey hyperedge-split point) is only meaningful
        // if it still has both a reactant coming in and at least one product
        // going out. If every product on one side got pruned, drop the junction
        // itself along with its now-dangling stub edge(s).
        const junctionIncoming = new Set();
        const junctionOutgoing = new Set();
        this.filteredEdges.forEach(e => {
            if (e.isToJunction) junctionIncoming.add(e.to);
            if (e.isFromJunction) junctionOutgoing.add(e.from);
        });
        const validJunctions = new Set(
            [...junctionIncoming].filter(id => junctionOutgoing.has(id))
        );

        this.filteredEdges = this.filteredEdges.filter(e => {
            if (e.isToJunction) return validJunctions.has(e.to);
            if (e.isFromJunction) return validJunctions.has(e.from);
            return true;
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
        const visibleJunctions = this.junctionNodes.filter(j => validJunctions.has(j.id));

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
