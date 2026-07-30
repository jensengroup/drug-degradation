/**
 * UI interactions and updates
 */
const UI = {
    currentPath: null,
    profileExpanded: true,

    init() {
        this.bindControls();
        this.updateStats();
    },

    bindControls() {
        // Filter sliders
        ['countThreshold', 'pathEnergyThreshold'].forEach(id => {
            const slider = document.getElementById(id);
            const valueSpan = document.getElementById(id + 'Value');
            if (slider && valueSpan) {
                slider.addEventListener('input', () => {
                    valueSpan.textContent = slider.value;
                });
                slider.addEventListener('change', () => this.applyFilters());
            }
        });

        // Disable filters checkbox
        document.getElementById('disableFilters')?.addEventListener('change', (e) => {
            if (e.target.checked) {
                this.disableFilters();
            } else {
                this.applyFilters();
            }
        });

        // Search
        document.getElementById('searchButton')?.addEventListener('click', () => this.searchNode());
        document.getElementById('nodeSearch')?.addEventListener('keypress', (e) => {
            if (e.key === 'Enter') this.searchNode();
        });

        // Grid view
        document.getElementById('gridViewButton')?.addEventListener('click', () => this.showGridModal());
        document.getElementById('updateGridButton')?.addEventListener('click', () => this.updateGridView());

        // Help
        document.getElementById('helpButton')?.addEventListener('click', () => this.showHelpModal());

        // Modal close events
        document.querySelectorAll('.modal').forEach(modal => {
            modal.addEventListener('click', (e) => {
                if (e.target === modal) this.hideModals();
            });
        });

        document.addEventListener('keydown', (e) => {
            if (e.key === 'Escape') this.hideModals();
        });

        // Chart modal click
        document.getElementById('reactionProfileChart')?.addEventListener('click', () => {
            if (this.currentPath) this.showChartModal();
        });
    },

    applyFilters() {
        const pathEnergyRaw = parseFloat(document.getElementById('pathEnergyThreshold')?.value);
        const maxPathEnergy = isNaN(pathEnergyRaw) ? Infinity : pathEnergyRaw;
        const countRaw = parseInt(document.getElementById('countThreshold')?.value);
        const minCount = isNaN(countRaw) ? 1 : countRaw;

        const { nodes, edges, junctions } = Data.filter({ maxPathEnergy, minCount });
        NetworkGraph.update(nodes, edges, junctions);
        this.updateStats();
    },

    disableFilters() {
        const { nodes, edges, junctions } = Data.reset();
        NetworkGraph.update(nodes, edges, junctions);
        this.updateStats();
    },

    updateStats() {
        document.getElementById('nodeCount').textContent = Data.filteredNodes.length;
        document.getElementById('edgeCount').textContent = Data.filteredEdges.length;
    },

    searchNode() {
        const input = document.getElementById('nodeSearch');
        const result = document.getElementById('searchResult');
        const nodeId = parseInt(input?.value);

        if (!nodeId || isNaN(nodeId)) {
            result.textContent = 'Please enter a valid node ID';
            return;
        }

        const node = Data.getNode(nodeId);
        if (!node) {
            result.textContent = `Node ${nodeId} not found`;
            return;
        }

        const isVisible = Data.filteredNodes.some(n => n.id === nodeId);
        if (!isVisible) {
            result.textContent = `Node ${nodeId} is filtered out`;
            return;
        }

        result.textContent = `Found node ${nodeId}`;
        NetworkGraph.selectNode(nodeId);
        NetworkGraph.focusOnNode(nodeId);
    },

    updateMoleculeDisplay() {
        const container = document.getElementById('moleculeContainer');
        if (!container) return;

        const selectedNodes = NetworkGraph.selectedNodesInOrder;

        if (selectedNodes.length === 0) {
            container.innerHTML = `
                <div class="no-molecules-message">
                    <p>Click a node to select</p>
                    <p class="no-molecules-hint">Ctrl+click to multi-select<br>Double-click for pathway</p>
                </div>
            `;
            return;
        }

        container.innerHTML = selectedNodes.map(nodeId => {
            const node = Data.getNode(nodeId);
            if (!node) return '';

            const title = node.isSubstrate ? `Substrate (${nodeId})` : `Node ${nodeId}`;
            const energy = node.relEnergy.toFixed(2);

            return `
                <div class="molecule-item">
                    <div class="molecule-header">
                        <span>${title}</span>
                        <span class="molecule-energy-info">∆E: ${energy} kcal/mol</span>
                    </div>
                    <div class="molecule-svg" onclick="UI.showMoleculeModal(${nodeId})">
                        ${node.svg || '<p style="color: #999;">No structure</p>'}
                    </div>
                    <div style="padding: 8px; text-align: center;">
                        <button class="btn" onclick="UI.showMoleculeModal(${nodeId})">Expand</button>
                        <button class="btn btn-secondary" onclick="UI.showPathwayModal(${nodeId})">Pathways</button>
                    </div>
                </div>
            `;
        }).join('');
    },

    showReactionProfile(path) {
        const pathways = Data.getPathways(path[path.length - 1]);
        if (!pathways || pathways.length === 0) return;

        // Find the pathway matching our path
        let profile = null;
        for (const p of pathways) {
            if (JSON.stringify(p.path) === JSON.stringify(path)) {
                profile = p.profile;
                break;
            }
        }

        if (!profile) {
            // Use best pathway
            const best = pathways.reduce((a, b) => a.barrier_height < b.barrier_height ? a : b);
            profile = best.profile;
        }

        this.currentPath = path;

        document.getElementById('reactionProfileContainer')?.classList.add('visible');
        ChartModule.createSidebarChart(profile);
    },

    hideReactionProfile() {
        this.currentPath = null;
        document.getElementById('reactionProfileContainer')?.classList.remove('visible');
        ChartModule.destroy('sidebar');
    },

    toggleReactionProfile() {
        const chart = document.getElementById('reactionProfileChart');
        const toggle = document.getElementById('reactionProfileToggle');
        if (chart && toggle) {
            this.profileExpanded = !this.profileExpanded;
            chart.style.display = this.profileExpanded ? 'block' : 'none';
            toggle.textContent = this.profileExpanded ? '▼' : '▶';
        }
    },

    showMoleculeModal(nodeId) {
        const node = Data.getNode(nodeId);
        if (!node) return;

        const modal = document.getElementById('moleculeModal');
        const content = document.getElementById('modalMoleculeContent');

        const title = node.isSubstrate ? `Substrate (Node ${nodeId})` : `Node ${nodeId}`;

        content.innerHTML = `
            <h2 style="margin: 0 0 10px 0; color: #333; font-size: 18px;">${title}</h2>
            <div style="margin-bottom: 10px; color: #666; font-size: 14px;">
                <div>∆E: ${node.relEnergy.toFixed(2)} kcal/mol</div>
                <div>Charge: ${node.charge} | Weight: ${node.weight.toFixed(1)} Da</div>
            </div>
            <div style="max-height: 70vh; overflow: auto;">
                ${node.svg || '<p>No structure available</p>'}
            </div>
        `;

        // Scale SVG
        const svgElement = content.querySelector('svg');
        if (svgElement) {
            svgElement.style.cssText = 'max-width: 90vw; height: auto; max-height: 60vh;';
        }

        modal.classList.add('active');
    },

    showChartModal() {
        if (!this.currentPath) return;

        const targetId = this.currentPath[this.currentPath.length - 1];
        const pathways = Data.getPathways(targetId);
        if (!pathways || pathways.length === 0) return;

        const best = pathways.reduce((a, b) => a.barrier_height < b.barrier_height ? a : b);

        document.getElementById('modalChartTitle').textContent =
            `Reaction Profile (${this.currentPath.length} nodes)`;

        document.getElementById('chartModal').classList.add('active');

        setTimeout(() => {
            ChartModule.createModalChart(best.profile);
        }, 100);
    },

    showPathwayModal(nodeId) {
        const pathways = Data.getPathways(nodeId);
        if (!pathways || pathways.length === 0) {
            alert('No pathways found for this node');
            return;
        }

        // For now, just show the best pathway as a chart
        const best = pathways.reduce((a, b) => a.barrier_height < b.barrier_height ? a : b);

        document.getElementById('modalChartTitle').textContent =
            `Best Pathway to Node ${nodeId} (barrier: ${best.barrier_height.toFixed(2)} kcal/mol)`;

        document.getElementById('chartModal').classList.add('active');

        setTimeout(() => {
            ChartModule.createModalChart(best.profile);
        }, 100);
    },

    showGridModal() {
        document.getElementById('gridViewModal').classList.add('active');
        this.updateGridView();
    },

    hideGridModal() {
        document.getElementById('gridViewModal').classList.remove('active');
    },

    updateGridView() {
        const container = document.getElementById('gridContainer');
        if (!container) return;

        const maxEnergy = parseFloat(document.getElementById('gridEnergyThreshold')?.value) || Infinity;
        const neutralOnly = document.getElementById('neutralOnlyCheckbox')?.checked;
        const minWeight = parseFloat(document.getElementById('minMolWeight')?.value) || 0;
        const maxWeight = parseFloat(document.getElementById('maxMolWeight')?.value) || Infinity;

        const filtered = Data.nodes.filter(n => {
            if (n.isSubstrate) return false;
            if (n.relEnergy > maxEnergy) return false;
            if (neutralOnly && n.charge !== 0) return false;
            if (n.weight < minWeight || n.weight > maxWeight) return false;
            return true;
        });

        // Sort by energy
        filtered.sort((a, b) => a.relEnergy - b.relEnergy);

        if (filtered.length === 0) {
            container.innerHTML = '<p style="text-align: center; color: #999;">No molecules match the filters</p>';
            return;
        }

        container.innerHTML = filtered.map(node => `
            <div class="grid-item" onclick="UI.showMoleculeModal(${node.id})">
                <div class="grid-item-title">Node ${node.id}</div>
                <div class="grid-item-energy">∆E: ${node.relEnergy.toFixed(2)} kcal/mol</div>
                <div class="grid-item-svg">${node.svg || ''}</div>
            </div>
        `).join('');
    },

    showHelpModal() {
        document.getElementById('helpModal').classList.add('active');
    },

    hideHelpModal() {
        document.getElementById('helpModal').classList.remove('active');
    },

    hideModals() {
        document.querySelectorAll('.modal').forEach(m => m.classList.remove('active'));
    }
};
