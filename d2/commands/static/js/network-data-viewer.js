/**
 * Network Data Viewer - handles viewing and exporting current network nodes and edges
 */
class NetworkDataViewer {
    constructor(networkApp) {
        this.app = networkApp;
        this.setupEventListeners();
    }

    setupEventListeners() {
        // Data viewing functionality removed

        // Export buttons (if they exist)
        const exportNodesBtn = document.getElementById('exportNodesButton');
        const exportEdgesBtn = document.getElementById('exportEdgesButton');
        const copyNetworkBtn = document.getElementById('copyNetworkDataButton');

        if (exportNodesBtn) {
            exportNodesBtn.addEventListener('click', () => this.exportNodes());
        }
        if (exportEdgesBtn) {
            exportEdgesBtn.addEventListener('click', () => this.exportEdges());
        }
        if (copyNetworkBtn) {
            copyNetworkBtn.addEventListener('click', () => this.copyNetworkData());
        }

        // Search functionality
        document.getElementById('nodeListSearch').addEventListener('input', (e) => this.filterNodeList(e.target.value));
        document.getElementById('edgeListSearch').addEventListener('input', (e) => this.filterEdgeList(e.target.value));
    }

    showNodeList() {
        const currentNodes = this.getCurrentNodes();
        this.populateNodeList(currentNodes);
        document.getElementById('nodeListModal').style.display = 'block';
    }

    showEdgeList() {
        const currentEdges = this.getCurrentEdges();
        this.populateEdgeList(currentEdges);
        document.getElementById('edgeListModal').style.display = 'block';
    }

    getCurrentNodes() {
        const nodes = this.app.network.body.data.nodes.get();
        return nodes.map(node => ({
            id: node.id,
            label: node.label,
            energy: this.app.networkData.relEnergyMap[node.id] || 0,
            charge: this.app.networkData.molChargeMap[node.id] || 0,
            weight: this.app.networkData.molWeightMap[node.id] || 0,
            isStereoisomer: node.label && node.label.includes('*'),
            groupInfo: this.app.stereoisomerManager?.getStereoisomerGroupInfo(node.id)
        })).sort((a, b) => a.id - b.id);
    }

    getCurrentEdges() {
        const edges = this.app.network.body.data.edges.get();
        // Get edge data from the current consolidated or original edges
        let edgeData = [];
        if (this.app.stereoisomerManager && this.app.stereoisomerManager.getCurrentEdges()) {
            edgeData = this.app.stereoisomerManager.getCurrentEdges().edgeData;
        } else if (this.app.edgeData) {
            edgeData = this.app.edgeData;
        }

        return edges.map(edge => {
            const edgeInfo = edgeData.find(e => e.id === edge.id);
            return {
                id: edge.id,
                from: edge.from,
                to: edge.to,
                type: edgeInfo?.type || 'unknown',
                count: edgeInfo?.count || 0,
                deltaE: edgeInfo?.deltaE || 0,
                label: `${edge.from} → ${edge.to}`,
                isConsolidated: edge.id.includes('consolidated_')
            };
        }).sort((a, b) => {
            if (a.from !== b.from) return a.from - b.from;
            return a.to - b.to;
        });
    }

    populateNodeList(nodes) {
        const container = document.getElementById('nodeListContent');
        const countElement = document.getElementById('nodeListCount');

        countElement.textContent = `${nodes.length} nodes`;

        container.innerHTML = nodes.map(node => `
            <div class="data-item" data-node-id="${node.id}">
                <div class="data-item-main">
                    <span class="data-item-id">Node ${node.id}${node.isStereoisomer ? '*' : ''}</span>
                    <div class="data-item-details">
                        ${node.groupInfo && node.groupInfo.groupSize > 1 ?
                `Group of ${node.groupInfo.groupSize} stereoisomers` :
                `Charge: ${node.charge > 0 ? '+' : ''}${node.charge}, MW: ${node.weight} Da`
            }
                    </div>
                    <span class="data-item-energy">ΔE: ${node.energy.toFixed(2)} kcal/mol</span>
                </div>
                <div class="data-item-actions">
                    <button class="data-item-button" onclick="networkDataViewer.focusNode(${node.id})">Focus</button>
                    <button class="data-item-button" onclick="networkDataViewer.selectNode(${node.id})">Select</button>
                </div>
            </div>
        `).join('');
    }

    populateEdgeList(edges) {
        const container = document.getElementById('edgeListContent');
        const countElement = document.getElementById('edgeListCount');

        countElement.textContent = `${edges.length} edges`;

        container.innerHTML = edges.map(edge => `
            <div class="data-item" data-edge-id="${edge.id}">
                <div class="data-item-main">
                    <span class="data-item-id">${edge.label}</span>
                    <div class="data-item-details">
                        ${edge.type}${edge.isConsolidated ? ' (consolidated)' : ''}
                        ${edge.count > 0 ? `, Count: ${edge.count}` : ''}
                    </div>
                    <span class="data-item-energy">ΔE: ${edge.deltaE.toFixed(2)} kcal/mol</span>
                </div>
                <div class="data-item-actions">
                    <button class="data-item-button" onclick="networkDataViewer.focusEdge('${edge.id}')">Focus</button>
                    <button class="data-item-button" onclick="networkDataViewer.selectEdge('${edge.id}')">Select</button>
                </div>
            </div>
        `).join('');
    }

    filterNodeList(searchTerm) {
        const items = document.querySelectorAll('#nodeListContent .data-item');
        const term = searchTerm.toLowerCase();

        let visibleCount = 0;
        items.forEach(item => {
            const text = item.textContent.toLowerCase();
            const isVisible = text.includes(term);
            item.style.display = isVisible ? 'flex' : 'none';
            if (isVisible) visibleCount++;
        });

        document.getElementById('nodeListCount').textContent =
            `${visibleCount} of ${items.length} nodes${searchTerm ? ' (filtered)' : ''}`;
    }

    filterEdgeList(searchTerm) {
        const items = document.querySelectorAll('#edgeListContent .data-item');
        const term = searchTerm.toLowerCase();

        let visibleCount = 0;
        items.forEach(item => {
            const text = item.textContent.toLowerCase();
            const isVisible = text.includes(term);
            item.style.display = isVisible ? 'flex' : 'none';
            if (isVisible) visibleCount++;
        });

        document.getElementById('edgeListCount').textContent =
            `${visibleCount} of ${items.length} edges${searchTerm ? ' (filtered)' : ''}`;
    }

    focusNode(nodeId) {
        // Focus on the node in the network
        this.app.network.focus(nodeId, {
            scale: 1.5,
            animation: {
                duration: 1000,
                easingFunction: 'easeInOutQuad'
            }
        });

        // Highlight the node temporarily
        this.app.network.selectNodes([nodeId]);
        setTimeout(() => {
            this.app.network.unselectAll();
        }, 2000);
    }

    selectNode(nodeId) {
        // Select the node and keep it selected
        this.app.network.selectNodes([nodeId]);
        this.app.network.focus(nodeId, {
            scale: 1.5,
            animation: {
                duration: 1000,
                easingFunction: 'easeInOutQuad'
            }
        });
    }

    focusEdge(edgeId) {
        const edges = this.app.network.body.data.edges;
        const edge = edges.get(edgeId);

        if (edge) {
            // Focus on the edge by centering between its nodes
            this.app.network.fit({
                nodes: [edge.from, edge.to],
                animation: {
                    duration: 1000,
                    easingFunction: 'easeInOutQuad'
                }
            });

            // Select the edge temporarily
            this.app.network.selectEdges([edgeId]);
            setTimeout(() => {
                this.app.network.unselectAll();
            }, 2000);
        }
    }

    selectEdge(edgeId) {
        const edges = this.app.network.body.data.edges;
        const edge = edges.get(edgeId);

        if (edge) {
            // Select the edge and keep it selected
            this.app.network.selectEdges([edgeId]);
            this.app.network.fit({
                nodes: [edge.from, edge.to],
                animation: {
                    duration: 1000,
                    easingFunction: 'easeInOutQuad'
                }
            });
        }
    }

    exportNodes() {
        const nodes = this.getCurrentNodes();
        const csvData = this.nodesToCSV(nodes);
        this.downloadCSV(csvData, 'network_nodes.csv');
    }

    exportEdges() {
        const edges = this.getCurrentEdges();
        const csvData = this.edgesToCSV(edges);
        this.downloadCSV(csvData, 'network_edges.csv');
    }

    nodesToCSV(nodes) {
        const headers = ['Node_ID', 'Label', 'Energy_kcal_mol', 'Charge', 'Molecular_Weight_Da', 'Is_Stereoisomer_Group', 'Group_Size'];
        const rows = nodes.map(node => [
            node.id,
            node.label,
            node.energy.toFixed(6),
            node.charge,
            node.weight,
            node.isStereoisomer ? 'Yes' : 'No',
            node.groupInfo?.groupSize || 1
        ]);

        return [headers, ...rows].map(row =>
            row.map(cell => `"${cell}"`).join(',')
        ).join('\n');
    }

    edgesToCSV(edges) {
        const headers = ['Edge_ID', 'From_Node', 'To_Node', 'Type', 'Count', 'DeltaE_kcal_mol', 'Is_Consolidated'];
        const rows = edges.map(edge => [
            edge.id,
            edge.from,
            edge.to,
            edge.type,
            edge.count,
            edge.deltaE.toFixed(6),
            edge.isConsolidated ? 'Yes' : 'No'
        ]);

        return [headers, ...rows].map(row =>
            row.map(cell => `"${cell}"`).join(',')
        ).join('\n');
    }

    downloadCSV(csvData, filename) {
        const blob = new Blob([csvData], { type: 'text/csv;charset=utf-8;' });
        const link = document.createElement('a');
        if (link.download !== undefined) {
            const url = URL.createObjectURL(blob);
            link.setAttribute('href', url);
            link.setAttribute('download', filename);
            link.style.visibility = 'hidden';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }
    }

    showNetworkData() {
        // Implement the logic to gather and display network-wide data
        const networkData = this.getNetworkData();
        this.populateNetworkData(networkData);
        document.getElementById('networkDataModal').style.display = 'block';
    }

    getNetworkData() {
        // Placeholder for network data retrieval logic
        return {
            totalNodes: this.app.network.body.data.nodes.length,
            totalEdges: this.app.network.body.data.edges.length,
            // Add more aggregated data as needed
        };
    }

    populateNetworkData(data) {
        const container = document.getElementById('networkDataContent');
        container.innerHTML = `
            <div>Total Nodes: ${data.totalNodes}</div>
            <div>Total Edges: ${data.totalEdges}</div>
            <!-- Add more data fields as needed -->
        `;
    }

    copyNetworkData() {
        const data = this.getNetworkData();
        const text = `Total Nodes: ${data.totalNodes}\nTotal Edges: ${data.totalEdges}`;
        navigator.clipboard.writeText(text).then(() => {
            alert('Network data copied to clipboard');
        }, () => {
            alert('Failed to copy network data');
        });
    }
}

// Global functions for modal control
function hideNodeListModal() {
    document.getElementById('nodeListModal').style.display = 'none';
}

function hideEdgeListModal() {
    document.getElementById('edgeListModal').style.display = 'none';
}

function hideNetworkDataModal() {
    document.getElementById('networkDataModal').style.display = 'none';
}

// Global variable for access from HTML
let networkDataViewer = null;
