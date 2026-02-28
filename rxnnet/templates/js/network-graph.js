/**
 * Network graph visualization using vis.js
 */
const NetworkGraph = {
    network: null,
    container: null,
    nodesDataset: null,
    edgesDataset: null,
    selectedNodes: new Set(),
    selectedNodesInOrder: [],
    originalNodeColors: new Map(),
    originalEdgeColors: new Map(),

    init(containerId) {
        this.container = document.getElementById(containerId);
        this.nodesDataset = new vis.DataSet();
        this.edgesDataset = new vis.DataSet();

        const options = {
            nodes: {
                shape: 'circle',
                size: 16,
                font: { size: 12, color: '#000000' },
                borderWidth: 2
            },
            edges: {
                width: 2,
                arrows: { to: { enabled: true, scaleFactor: 0.8 } },
                smooth: false // Straight edges for cleaner hyperedge appearance
            },
            physics: {
                enabled: true,
                solver: 'forceAtlas2Based',
                forceAtlas2Based: {
                    gravitationalConstant: -50,
                    centralGravity: 0.01,
                    springLength: 100,
                    springConstant: 0.08
                },
                stabilization: { iterations: 150 }
            },
            interaction: {
                multiselect: true,
                selectConnectedEdges: false,
                hover: true,
                tooltipDelay: 200
            }
        };

        this.network = new vis.Network(
            this.container,
            { nodes: this.nodesDataset, edges: this.edgesDataset },
            options
        );

        this.bindEvents();
    },

    bindEvents() {
        this.network.on('click', (params) => {
            if (params.nodes.length > 0) {
                const nodeId = params.nodes[0];
                // Skip junction nodes
                if (nodeId < 0) return;

                const ctrlKey = params.event.srcEvent && (params.event.srcEvent.ctrlKey || params.event.srcEvent.metaKey);

                if (ctrlKey) {
                    this.toggleNodeSelection(nodeId);
                } else {
                    this.selectNode(nodeId);
                }
            } else {
                this.clearSelection();
            }
        });

        this.network.on('doubleClick', (params) => {
            if (params.nodes.length > 0) {
                const nodeId = params.nodes[0];
                if (nodeId < 0) return; // Skip junctions
                this.selectPathway(nodeId);
            }
        });
    },

    update(nodes, edges, junctions = []) {
        // Convert molecule nodes to vis.js format
        const visNodes = nodes.map(n => {
            const node = {
                id: n.id,
                label: n.isSubstrate ? 'Substrate' : n.id.toString(),
                color: this.getNodeColor(n),
                shape: n.isSubstrate ? 'ellipse' : 'circle',
                font: n.isSubstrate ? { color: 'white' } : {},
                title: this.getNodeTooltip(n)
            };
            this.originalNodeColors.set(n.id, node.color);
            return node;
        });

        // Add junction nodes (small dots that stay between reactant and products)
        junctions.forEach(j => {
            visNodes.push({
                id: j.id,
                label: '',
                shape: 'dot',
                size: 4,
                color: { background: '#666', border: '#666' },
                title: `Hyperedge: ${j.reactantId} → [${j.productIds.join(', ')}]`,
                fixed: false,
                physics: true
            });
        });

        // Convert edges
        const visEdges = edges.map(e => {
            const edge = {
                id: e.id,
                from: e.from,
                to: e.to,
                color: this.getEdgeColor(e),
                width: this.getEdgeWidth(e),
                title: this.getEdgeTooltip(e)
            };

            if (e.isToJunction) {
                edge.arrows = { to: { enabled: false } };
                edge.length = 50; // Half length for reactant->junction
            }

            if (e.isFromJunction) {
                edge.length = 50; // Half length for junction->products
            }

            if (e.isEliminated) {
                edge.dashes = [4, 4];
            }

            this.originalEdgeColors.set(e.id, edge.color);
            return edge;
        });

        this.nodesDataset.clear();
        this.edgesDataset.clear();
        this.nodesDataset.add(visNodes);
        this.edgesDataset.add(visEdges);
    },

    getNodeColor(node) {
        if (node.isSubstrate) return CONFIG.NODE_COLORS.substrate;
        if (this.selectedNodes.has(node.id)) return CONFIG.NODE_COLORS.selected;
        return CONFIG.NODE_COLORS.default;
    },

    getEdgeColor(edge) {
        const type = (edge.type || '').toLowerCase();
        return CONFIG.EDGE_COLORS[type] || CONFIG.EDGE_COLORS.default;
    },

    getEdgeWidth(edge) {
        if (edge.type === 'reaction' || edge.type === 'mtd-reaction') {
            return Math.min(Math.max(1, (edge.count || 1)), 10) / 3;
        }
        return 1.5;
    },

    getNodeTooltip(node) {
        let tip = `Node ${node.id}`;
        tip += `\n∆E: ${node.relEnergy.toFixed(2)} kcal/mol`;
        tip += `\nCharge: ${node.charge}`;
        if (node.barrier !== null) {
            tip += `\nBarrier: ${node.barrier.toFixed(2)} kcal/mol`;
        }
        return tip;
    },

    getEdgeTooltip(edge) {
        let tip = edge.type.charAt(0).toUpperCase() + edge.type.slice(1);
        if (edge.count > 1) {
            tip += ` (${edge.count}x)`;
        }
        return tip;
    },

    selectNode(nodeId) {
        this.clearSelection();
        this.selectedNodes.add(nodeId);
        this.selectedNodesInOrder = [nodeId];
        this.highlightSelection();
        UI.updateMoleculeDisplay();
    },

    toggleNodeSelection(nodeId) {
        if (this.selectedNodes.has(nodeId)) {
            this.selectedNodes.delete(nodeId);
            this.selectedNodesInOrder = this.selectedNodesInOrder.filter(id => id !== nodeId);
        } else {
            this.selectedNodes.add(nodeId);
            this.selectedNodesInOrder.push(nodeId);
        }
        this.highlightSelection();
        UI.updateMoleculeDisplay();
    },

    selectPathway(targetNodeId) {
        const path = Data.getBestPath(targetNodeId);
        if (!path || path.length === 0) {
            // No pathway, just select the node
            this.selectNode(targetNodeId);
            return;
        }

        this.clearSelection();
        path.forEach(nodeId => {
            this.selectedNodes.add(nodeId);
            this.selectedNodesInOrder.push(nodeId);
        });
        this.highlightSelection();
        this.highlightPathEdges(path);
        UI.updateMoleculeDisplay();
        UI.showReactionProfile(path);
    },

    clearSelection() {
        this.selectedNodes.clear();
        this.selectedNodesInOrder = [];
        this.resetHighlights();
        UI.updateMoleculeDisplay();
        UI.hideReactionProfile();
    },

    highlightSelection() {
        const updates = [];
        this.nodesDataset.forEach(node => {
            if (node.id > 0) { // Skip junctions
                const isSelected = this.selectedNodes.has(node.id);
                const originalColor = this.originalNodeColors.get(node.id);
                updates.push({
                    id: node.id,
                    color: isSelected ? CONFIG.SELECTED_HIGHLIGHT_COLOR : originalColor,
                    borderWidth: isSelected ? 3 : 2
                });
            }
        });
        this.nodesDataset.update(updates);
    },

    highlightPathEdges(path) {
        const edgeUpdates = [];
        for (let i = 0; i < path.length - 1; i++) {
            const from = path[i];
            const to = path[i + 1];
            // Find edge between these nodes
            this.edgesDataset.forEach(edge => {
                if ((edge.from === from && edge.to === to) ||
                    (edge.from === to && edge.to === from)) {
                    edgeUpdates.push({
                        id: edge.id,
                        color: CONFIG.SELECTION_COLOR,
                        width: 3
                    });
                }
            });
        }
        this.edgesDataset.update(edgeUpdates);
    },

    resetHighlights() {
        // Reset nodes
        const nodeUpdates = [];
        this.nodesDataset.forEach(node => {
            if (node.id > 0) {
                nodeUpdates.push({
                    id: node.id,
                    color: this.originalNodeColors.get(node.id),
                    borderWidth: 2
                });
            }
        });
        this.nodesDataset.update(nodeUpdates);

        // Reset edges
        const edgeUpdates = [];
        this.edgesDataset.forEach(edge => {
            edgeUpdates.push({
                id: edge.id,
                color: this.originalEdgeColors.get(edge.id),
                width: this.getEdgeWidth(Data.edges.find(e => e.id === edge.id) || {})
            });
        });
        this.edgesDataset.update(edgeUpdates);
    },

    focusOnNode(nodeId) {
        this.network.focus(nodeId, {
            scale: 1.5,
            animation: {
                duration: CONFIG.ANIMATION.FOCUS_DURATION,
                easingFunction: 'easeInOutQuad'
            }
        });
    }
};
