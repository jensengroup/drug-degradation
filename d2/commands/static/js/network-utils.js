/**
 * Network statistics and utility functions
 */
class NetworkStats {
    constructor(network) {
        this.network = network;
    }

    /**
     * Update graph statistics display
     */
    updateGraphStats() {
        const nodes = this.network.body.data.nodes;
        const edges = this.network.body.data.edges;
        const currentNodes = nodes.get();
        const currentEdges = edges.get();

        // Update node count
        document.getElementById('nodeCount').textContent = currentNodes.length;

        // Count edges by type
        const edgeTypeCounts = {};
        currentEdges.forEach(edge => {
            // Get edge type from original edge data
            const edgeDataInfo = window.edgeData ? window.edgeData.find(e => e.id === edge.id) : null;
            const edgeType = edgeDataInfo ? edgeDataInfo.type : 'unknown';
            edgeTypeCounts[edgeType] = (edgeTypeCounts[edgeType] || 0) + 1;
        });

        // Update total edge count
        document.getElementById('edgeCount').textContent = currentEdges.length;

        // Update edge breakdown
        const edgeBreakdown = document.getElementById('edgeBreakdown');
        const breakdownLines = [];
        for (const [type, count] of Object.entries(edgeTypeCounts).sort()) {
            const capitalizedType = type.charAt(0).toUpperCase() + type.slice(1);
            const color = this.getEdgeColor(type);
            breakdownLines.push(`<div class="edge-type-item clickable" data-edge-type="${type}" title="Click to view ${capitalizedType} edge indices">
                <span style="display: inline-block; width: 8px; height: 8px; border-radius: 50%; background-color: ${color}; margin-right: 4px;"></span>
                ${capitalizedType}: ${count}
            </div>`);
        }
        edgeBreakdown.innerHTML = breakdownLines.join('') || 'No edges';

        // Add click event listeners to edge type items
        const edgeTypeItems = document.querySelectorAll('.edge-type-item');
        edgeTypeItems.forEach(item => {
            item.addEventListener('click', () => {
                const edgeType = item.getAttribute('data-edge-type');
                window.networkApp.ui.showIndicesModal('edges', edgeType);
            });
        });
    }

    /**
     * Get edge color based on reaction type
     */
    getEdgeColor(rxnType) {
        const type = rxnType.toLowerCase();
        return CONFIG.EDGE_COLORS[type] || CONFIG.EDGE_COLORS.default;
    }
}

/**
 * Utility functions
 */
const NetworkUtils = {
    /**
     * Update node tooltips only
     */
    updateNodeTooltipsOnly(network, networkData) {
        const nodes = network.body.data.nodes;
        const currentNodes = nodes.get();
        const updatedNodes = [];

        currentNodes.forEach(node => {
            const relEnergy = networkData.relEnergyMap[node.id] || 0;

            let title;
            if (node.id == networkData.substrateId) {
                title = `Substrate (Node ${node.id})\n∆${networkData.energyType}: 0.00 kcal/mol`;
            } else {
                title = `Node ${node.id}\n∆${networkData.energyType}: ${relEnergy.toFixed(2)} kcal/mol`;
            }

            updatedNodes.push({
                id: node.id,
                title: title
            });
        });

        if (updatedNodes.length > 0) {
            nodes.update(updatedNodes);
        }
    },

    /**
     * Update edge tooltips with static information
     */
    updateEdgeTooltips(network, networkData) {
        const edges = network.body.data.edges;
        const currentEdges = edges.get();
        const currentEdgeIds = new Set(currentEdges.map(e => e.id));
        const updatedEdges = [];

        // Use original edge data
        let edgeDataToUse = networkData.originalEdgeData;

        // Only update tooltips for edges that actually exist in the network
        edgeDataToUse.forEach(edgeInfo => {
            if (!currentEdgeIds.has(edgeInfo.id)) {
                return;
            }

            const productFragments = edgeInfo.smaller_products || [];

            // Create tooltip text with energy and barrier information
            let tooltipText = edgeInfo.type.charAt(0).toUpperCase() + edgeInfo.type.slice(1);
            if (edgeInfo.type.toLowerCase() === "reaction" && edgeInfo.count) {
                tooltipText += ` (${edgeInfo.count}x)`;
            }
            if (productFragments.length > 0) {
                tooltipText += `\nEliminated Products: ${productFragments.join(', ')}`;
            }
            // Include energy information if available
            if (edgeInfo.rxn_energy !== undefined && edgeInfo.rxn_energy !== null) {
                tooltipText += `\n∆${networkData.energyType} ${edgeInfo.rxn_energy.toFixed(2)} kcal/mol`;
            }
            // Include barrier information if available
            if (edgeInfo.barrier !== undefined && edgeInfo.barrier !== null) {
                tooltipText += `\nBarrier: ${edgeInfo.barrier.toFixed(2)} kcal/mol`;
            }

            updatedEdges.push({
                id: edgeInfo.id,
                title: tooltipText
            });
        });

        if (updatedEdges.length > 0) {
            edges.update(updatedEdges);
        }
    },

    /**
     * Create updated node data
     */
    createUpdatedNodeData(nodeId, originalNodeData, networkData, selectedNodes) {
        const nodeInfo = originalNodeData.find(n => n.id === nodeId);
        if (!nodeInfo) {
            console.warn(`Node info not found for nodeId: ${nodeId}`);
            return null;
        }

        // Use simple node label
        let nodeLabel = nodeId.toString();

        // Create updated tooltip
        const relEnergy = networkData.relEnergyMap[nodeId] || 0;

        let updatedTitle;
        if (nodeId == networkData.substrateId) {
            updatedTitle = `Substrate (Node ${nodeId})\n∆${networkData.energyType}: 0.00 kcal/mol`;
        } else {
            updatedTitle = `Node ${nodeId}\n∆${networkData.energyType}: ${relEnergy.toFixed(2)} kcal/mol`;
        }

        const updatedNodeInfo = { ...nodeInfo, title: updatedTitle, label: nodeLabel };

        // Restore selection state for re-added nodes
        if (selectedNodes.has(nodeId)) {
            return { ...updatedNodeInfo, color: CONFIG.SELECTED_HIGHLIGHT_COLOR };
        } else {
            return updatedNodeInfo;
        }
    },

    /**
     * Validate network state
     */
    validateNetworkState(network) {
        const nodes = network.body.data.nodes;
        const edges = network.body.data.edges;
        const postUpdateNodes = nodes.get();
        const postUpdateEdges = edges.get();
        const postUpdateNodeIds = new Set(postUpdateNodes.map(n => n.id));

        const orphanedEdges = [];
        const malformedEdges = [];

        postUpdateEdges.forEach(edge => {
            const edgeFrom = edge.from || edge.source || edge.fromId;
            const edgeTo = edge.to || edge.target || edge.toId;

            if (edgeFrom === undefined || edgeTo === undefined) {
                malformedEdges.push(edge.id);
                console.warn(`🔧 MALFORMED EDGE: ${edge.id} has undefined endpoints`);
                return;
            }

            if (!postUpdateNodeIds.has(edgeFrom) || !postUpdateNodeIds.has(edgeTo)) {
                orphanedEdges.push(edge.id);
                console.error(`🚨 ORPHANED EDGE: ${edge.id} (${edgeFrom} → ${edgeTo})`);
            }
        });

        // Clean up problematic edges
        const edgesToCleanup = [...malformedEdges, ...orphanedEdges];
        if (edgesToCleanup.length > 0) {
            edges.remove(edgesToCleanup);
        }
    }
};
