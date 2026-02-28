/**
 * UI components and modal management
 */
class NetworkUI {
    constructor(networkData, network) {
        this.networkData = networkData;
        this.network = network;
        this.selectedNodes = new Set();
        this.selectedNodesInOrder = []; // Track node selection order
        this.originalNodeColors = new Map();
        this.originalEdgeColors = new Map();
        this.highlightedEdges = new Set();
        this.lastSelectedNode = null;
        this.reactionChart = null;
        this.modalReactionChart = null;
        this.moleculeModalOpenedFromGrid = false;
        this.graphAdjacency = new Map();
    }

    /**
     * Initialize UI components
     */
    initialize(originalNodeData, originalNetworkEdgeData, graphAdjacency) {
        this.graphAdjacency = graphAdjacency;

        // Store original colors for selection highlighting
        originalNodeData.forEach(node => {
            this.originalNodeColors.set(node.id, node.color);
        });

        originalNetworkEdgeData.forEach(edge => {
            this.originalEdgeColors.set(edge.id, edge.color);
        });

        this.setupEventListeners();
        this.updateSelectionDisplay();
    }

    /**
     * Setup all event listeners
     */
    setupEventListeners() {
        // Modal event listeners
        document.getElementById('moleculeModal').addEventListener('click', (e) => {
            if (e.target === e.currentTarget) this.hideModals();
        });

        document.getElementById('chartModal').addEventListener('click', (e) => {
            if (e.target === e.currentTarget) this.hideModals();
        });

        document.getElementById('gridViewModal').addEventListener('click', (e) => {
            if (e.target === e.currentTarget) this.hideGridModal();
        });

        document.getElementById('helpModal').addEventListener('click', (e) => {
            if (e.target === e.currentTarget) this.hideHelpModal();
        });

        document.getElementById('indicesModal').addEventListener('click', (e) => {
            if (e.target === e.currentTarget) this.hideIndicesModal();
        });

        // Escape key to close modals
        document.addEventListener('keydown', (e) => {
            if (e.key === 'Escape') {
                if (document.getElementById('helpModal').style.display === 'block') {
                    this.hideHelpModal();
                } else if (document.getElementById('gridViewModal').style.display === 'block') {
                    this.hideGridModal();
                } else if (document.getElementById('indicesModal').style.display === 'block') {
                    this.hideIndicesModal();
                } else {
                    this.hideModals();
                }
            }
        });

        // Grid view event listeners
        document.getElementById('gridViewButton').addEventListener('click', () => this.showGridModal());
        document.getElementById('updateGridButton').addEventListener('click', () => this.updateGridView());
        document.getElementById('neutralOnlyCheckbox').addEventListener('change', () => this.updateGridView());
        document.getElementById('helpButton').addEventListener('click', () => this.showHelpModal());

        // Make graph statistics clickable for showing indices
        document.getElementById('nodeCountContainer').addEventListener('click', () => this.showIndicesModal('nodes'));
        document.getElementById('edgeCountContainer').addEventListener('click', () => this.showIndicesModal('edges'));

        // Grid input event listeners
        ['gridMinCount', 'gridMaxBarrier', 'gridEnergyThreshold', 'minMolWeight', 'maxMolWeight'].forEach(id => {
            document.getElementById(id).addEventListener('keypress', (e) => {
                if (e.key === 'Enter') this.updateGridView();
            });
        });

        // Search functionality
        document.getElementById("searchButton").addEventListener("click", () => this.searchNode());
        document.getElementById("nodeSearch").addEventListener("keypress", (e) => {
            if (e.key === "Enter") this.searchNode();
        });
        document.getElementById("nodeSearch").addEventListener("input", () => {
            document.getElementById("searchResult").textContent = "";
        });
    }

    /**
     * Show molecule modal
     */
    showMoleculeModal(nodeId, openedFromGrid = false) {
        const svg = this.networkData.svgMap[nodeId];
        const relEnergy = this.networkData.relEnergyMap[nodeId] || 0;
        const absEnergy = this.networkData.molEnergyMap[nodeId] || 0;

        if (!svg) return;

        this.moleculeModalOpenedFromGrid = openedFromGrid;

        const modal = document.getElementById('moleculeModal');
        const content = document.getElementById('modalMoleculeContent');

        const title = nodeId == this.networkData.substrateId ? `Substrate (Node ${nodeId})` : `Node ${nodeId}`;

        // Show absolute energy (G) for single node selection
        content.innerHTML = `
            <h2 style="margin: 0 0 10px 0; color: #333; font-size: 18px;">${title}</h2>
            <div style="margin-bottom: 10px; color: #666; font-size: 14px;">
                <div style="margin-bottom: 5px;"><strong>Absolute Energy (G):</strong> ${absEnergy.toFixed(6)} hartree</div>
                <div style="color: #888; font-size: 12px;">Relative Energy (∆${this.networkData.energyType}): ${relEnergy.toFixed(2)} kcal/mol</div>
            </div>
            <div style="max-height: 85vh; overflow: auto;">
                ${svg}
            </div>
        `;

        // Scale SVG for modal
        const svgElement = content.querySelector('svg');
        if (svgElement) {
            svgElement.style.cssText = 'max-width: 90vw; height: auto; max-height: 80vh;';
        }

        modal.style.display = 'block';
    }

    /**
     * Show chart modal
     */
    showChartModal() {
        if (!this.reactionChart || this.selectedNodes.size === 0) return;

        const nodeArray = [...this.selectedNodesInOrder];
        const modal = document.getElementById('chartModal');
        const modalTitle = document.getElementById('modalChartTitle');

        const nodeCount = this.selectedNodes.size;
        modalTitle.textContent = `Reaction Profile (${nodeCount} node${nodeCount > 1 ? 's' : ''})`;

        // Initialize modal chart if it doesn't exist
        if (!this.modalReactionChart) {
            this.initializeModalChart();
        }

        // Copy data from main chart
        this.modalReactionChart.data.labels = [...this.reactionChart.data.labels];
        this.modalReactionChart.data.datasets[0].data = [...this.reactionChart.data.datasets[0].data];

        // Update colors based on selection
        const hasSubstrate = nodeArray.includes(this.networkData.substrateId);
        if (hasSubstrate) {
            this.modalReactionChart.data.datasets[0].borderColor = CONFIG.CHART_COLORS.PRIMARY;
            this.modalReactionChart.data.datasets[0].backgroundColor = CONFIG.CHART_COLORS.PRIMARY_BACKGROUND;
            this.modalReactionChart.data.datasets[0].pointBackgroundColor = CONFIG.CHART_COLORS.PRIMARY;
        } else {
            this.modalReactionChart.data.datasets[0].borderColor = CONFIG.CHART_COLORS.SECONDARY;
            this.modalReactionChart.data.datasets[0].backgroundColor = CONFIG.CHART_COLORS.SECONDARY_BACKGROUND;
            this.modalReactionChart.data.datasets[0].pointBackgroundColor = CONFIG.CHART_COLORS.SECONDARY;
        }

        modal.style.display = 'block';

        setTimeout(() => {
            this.modalReactionChart.update('none');
            this.setupModalChartEvents(nodeArray);
        }, 100);
    }

    /**
     * Initialize modal chart
     */
    initializeModalChart() {
        const ctx = document.getElementById('modalReactionChart').getContext('2d');
        this.modalReactionChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: [],
                datasets: [{
                    label: `∆${this.networkData.energyType} (kcal/mol)`,
                    data: [],
                    borderColor: CONFIG.CHART_COLORS.PRIMARY,
                    backgroundColor: CONFIG.CHART_COLORS.PRIMARY_BACKGROUND,
                    pointBackgroundColor: CONFIG.CHART_COLORS.PRIMARY,
                    pointBorderColor: '#fff',
                    pointBorderWidth: 3,
                    pointRadius: 8,
                    pointHoverRadius: 10,
                    fill: true,
                    tension: 0.1,
                    borderWidth: 3
                }]
            },
            options: this.getModalChartOptions()
        });
    }

    /**
     * Get modal chart configuration options
     */
    getModalChartOptions() {
        return {
            responsive: true,
            maintainAspectRatio: false,
            layout: {
                padding: { bottom: 10 }
            },
            scales: {
                y: {
                    beginAtZero: false,
                    title: {
                        display: true,
                        text: `∆${this.networkData.energyType} (kcal/mol)`,
                        font: { size: 16 }
                    },
                    grid: { color: 'rgba(0,0,0,0.1)' },
                    ticks: { font: { size: 14 } }
                },
                x: {
                    title: {
                        display: true,
                        text: 'Node',
                        font: { size: 16 }
                    },
                    grid: { color: 'rgba(0,0,0,0.1)' },
                    ticks: {
                        font: { size: 14 },
                        color: '#4A90E2'
                    }
                }
            },
            plugins: {
                legend: { display: false },
                tooltip: {
                    titleFont: { size: 16 },
                    bodyFont: { size: 14 },
                    callbacks: {
                        title: (context) => {
                            const nodeId = context[0].label;
                            return nodeId == this.networkData.substrateId ? `Substrate (${nodeId})` : `Node ${nodeId}`;
                        },
                        label: (context) => {
                            return `∆${this.networkData.energyType}: ${context.parsed.y.toFixed(2)} kcal/mol`;
                        }
                    }
                }
            },
            interaction: { intersect: false, mode: 'index' },
            onHover: (event, elements) => this.handleModalChartHover(event, elements)
        };
    }

    /**
     * Handle modal chart hover events
     */
    handleModalChartHover(event, elements) {
        try {
            const canvasPosition = Chart.helpers.getRelativePosition(event, this.modalReactionChart);
            const dataX = this.modalReactionChart.scales.x.getValueForPixel(canvasPosition.x);
            const chartArea = this.modalReactionChart.chartArea;

            if (canvasPosition.y > chartArea.bottom) {
                const labelIndex = Math.round(dataX);
                const nodeArray = [...this.selectedNodesInOrder];
                if (labelIndex >= 0 && labelIndex < nodeArray.length) {
                    const nodeId = nodeArray[labelIndex];
                    const mouseX = event.x || event.clientX || 0;
                    const mouseY = event.y || event.clientY || 0;
                    this.showSmallMoleculeModal(nodeId, mouseX, mouseY);
                }
            } else {
                const smallModal = document.getElementById('smallMoleculeModal');
                if (smallModal) smallModal.style.display = 'none';
            }
        } catch (e) {
            console.warn('Error in onHover handler:', e);
            const smallModal = document.getElementById('smallMoleculeModal');
            if (smallModal) smallModal.style.display = 'none';
        }
    }

    /**
     * Setup modal chart event listeners
     */
    setupModalChartEvents(nodeArray) {
        const chartCanvas = document.getElementById('modalReactionChart');

        const handleMouseLeave = () => {
            const smallModal = document.getElementById('smallMoleculeModal');
            if (smallModal) smallModal.style.display = 'none';
        };

        const handleMouseMove = (event) => {
            try {
                const canvasPosition = Chart.helpers.getRelativePosition(event, this.modalReactionChart);
                const chartArea = this.modalReactionChart.chartArea;

                if (canvasPosition.y <= chartArea.bottom) {
                    const smallModal = document.getElementById('smallMoleculeModal');
                    if (smallModal) smallModal.style.display = 'none';
                }
            } catch (e) {
                console.warn('Error in mousemove handler:', e);
            }
        };

        // Remove existing listeners
        chartCanvas.removeEventListener('mouseleave', handleMouseLeave);
        chartCanvas.removeEventListener('mousemove', handleMouseMove);

        // Add new listeners
        chartCanvas.addEventListener('mouseleave', handleMouseLeave);
        chartCanvas.addEventListener('mousemove', handleMouseMove);
    }

    /**
     * Hide modals
     */
    hideModals() {
        if (this.moleculeModalOpenedFromGrid && document.getElementById('moleculeModal').style.display === 'block') {
            document.getElementById('moleculeModal').style.display = 'none';
            this.showGridModal();
            this.moleculeModalOpenedFromGrid = false;
            return;
        }

        document.getElementById('moleculeModal').style.display = 'none';
        document.getElementById('chartModal').style.display = 'none';
        this.moleculeModalOpenedFromGrid = false;
    }

    /**
     * Hide grid modal
     */
    hideGridModal() {
        document.getElementById('gridViewModal').style.display = 'none';
    }

    /**
     * Show help modal
     */
    showHelpModal() {
        document.getElementById('helpModal').style.display = 'block';
    }

    /**
     * Hide help modal
     */
    hideHelpModal() {
        document.getElementById('helpModal').style.display = 'none';
    }

    /**
     * Show grid modal
     */
    showGridModal() {
        const modal = document.getElementById('gridViewModal');
        modal.style.display = 'block';
        this.updateGridView();
    }

    /**
     * Show small molecule modal (tooltip)
     */
    showSmallMoleculeModal(nodeId, mouseX, mouseY) {
        try {
            const svg = this.networkData.svgMap[nodeId];
            const relEnergy = this.networkData.relEnergyMap[nodeId] || 0;
            const absEnergy = this.networkData.molEnergyMap[nodeId] || 0;

            if (!svg) return;

            const modal = document.getElementById('smallMoleculeModal');
            const content = document.getElementById('smallMoleculeContent');

            if (!modal || !content) return;

            const title = nodeId == this.networkData.substrateId ? `Substrate (Node ${nodeId})` : `Node ${nodeId}`;

            const tooltipWidth = window.innerWidth * CONFIG.TOOLTIP_WIDTH_FRACTION;
            const tooltipHeight = window.innerHeight * CONFIG.TOOLTIP_HEIGHT_FRACTION;

            // Show absolute energy (G) for single node tooltip
            content.innerHTML = `
                <div style="font-size: 14px; font-weight: bold; color: #333; margin-bottom: 6px; text-align: center;">${title}</div>
                <div style="font-size: 12px; color: #666; margin-bottom: 8px; text-align: center;">
                    <div style="font-weight: bold;">G: ${absEnergy.toFixed(6)} hartree</div>
                    <div style="color: #888; font-size: 11px;">∆${this.networkData.energyType}: ${relEnergy.toFixed(2)} kcal/mol</div>
                </div>
                <div style="flex: 1; display: flex; align-items: center; justify-content: center; overflow: hidden;">
                    ${svg}
                </div>
            `;

            // Scale SVG to fit the tooltip
            const svgElement = content.querySelector('svg');
            if (svgElement) {
                const availableHeight = tooltipHeight - CONFIG.TOOLTIP_TEXT_HEIGHT;
                const availableWidth = tooltipWidth;

                svgElement.style.cssText = `
                    width: 100%;
                    height: 100%;
                    max-width: ${availableWidth}px;
                    max-height: ${availableHeight}px;
                    object-fit: contain;
                `;
            }

            // Position tooltip
            modal.style.width = tooltipWidth + 'px';
            modal.style.height = tooltipHeight + 'px';

            const safeMouseX = mouseX || 100;
            const safeMouseY = mouseY || 100;

            let left = safeMouseX - tooltipWidth / 2;
            let top = safeMouseY - tooltipHeight / 2;

            // Adjust if tooltip would go off screen
            if (left < CONFIG.TOOLTIP_PADDING) left = CONFIG.TOOLTIP_PADDING;
            if (left + tooltipWidth > window.innerWidth - CONFIG.TOOLTIP_PADDING) {
                left = window.innerWidth - tooltipWidth - CONFIG.TOOLTIP_PADDING;
            }
            if (top < CONFIG.TOOLTIP_PADDING) top = CONFIG.TOOLTIP_PADDING;
            if (top + tooltipHeight > window.innerHeight - CONFIG.TOOLTIP_PADDING) {
                top = window.innerHeight - tooltipHeight - CONFIG.TOOLTIP_PADDING;
            }

            modal.style.left = left + 'px';
            modal.style.top = top + 'px';
            modal.style.display = 'block';
        } catch (e) {
            console.warn('Error showing small molecule modal:', e);
            const modal = document.getElementById('smallMoleculeModal');
            if (modal) modal.style.display = 'none';
        }
    }

    /**
     * Filter pathways by count and barrier criteria (JavaScript implementation of Python logic)
     */
    filterPathways(minCount = 1, maxBarrier = 30) {
        const filteredPathways = {};
        const precomputedPathways = this.networkData.precomputedPathways || {};

        console.log(`🔍 Filtering pathways with min_count=${minCount}, max_barrier=${maxBarrier}`);
        console.log(`📊 Total nodes with precomputed pathways: ${Object.keys(precomputedPathways).length}`);

        let totalPathways = 0;
        let validPathways = 0;

        for (const [nodeIdStr, pathways] of Object.entries(precomputedPathways)) {
            const nodeId = parseInt(nodeIdStr);
            const validPathwaysForNode = [];

            console.log(`🧪 Processing node ${nodeId} with ${pathways.length} pathways`);

            for (let i = 0; i < pathways.length; i++) {
                const pathway = pathways[i];
                totalPathways++;

                // Check pathway barrier
                let pathwayBarrier = 0;
                let pathwayCount = pathway.count || 1; // Default to 1 if count not specified

                // Use barrier_height if available, otherwise calculate from profile
                if (pathway.barrier_height !== undefined) {
                    pathwayBarrier = pathway.barrier_height;
                    console.log(`  Pathway ${i}: using barrier_height = ${pathwayBarrier} kcal/mol`);
                } else if (pathway.profile && Array.isArray(pathway.profile)) {
                    const energies = pathway.profile.map(step => step.energy || 0);
                    if (energies.length > 0) {
                        const minEnergy = Math.min(...energies);
                        const maxEnergy = Math.max(...energies);
                        pathwayBarrier = maxEnergy - minEnergy;
                        console.log(`  Pathway ${i}: calculated barrier = ${pathwayBarrier} kcal/mol (from ${minEnergy} to ${maxEnergy})`);
                    }
                } else {
                    console.log(`  Pathway ${i}: no barrier data available`);
                }

                console.log(`  Pathway ${i}: count=${pathwayCount}, barrier=${pathwayBarrier} kcal/mol`);

                // Check if pathway meets criteria
                const countCheck = pathwayCount >= minCount;
                const barrierCheck = pathwayBarrier <= maxBarrier;

                console.log(`  Pathway ${i}: count check (${pathwayCount} >= ${minCount}) = ${countCheck}`);
                console.log(`  Pathway ${i}: barrier check (${pathwayBarrier} <= ${maxBarrier}) = ${barrierCheck}`);

                if (countCheck && barrierCheck) {
                    validPathwaysForNode.push(pathway);
                    validPathways++;
                    console.log(`  ✅ Pathway ${i}: PASSED all criteria`);
                } else {
                    console.log(`  ❌ Pathway ${i}: FAILED criteria`);
                }
            }

            if (validPathwaysForNode.length > 0) {
                filteredPathways[nodeId] = validPathwaysForNode;
                console.log(`✅ Node ${nodeId}: ${validPathwaysForNode.length} valid pathways`);
            } else {
                console.log(`❌ Node ${nodeId}: no valid pathways`);
            }
        }

        console.log(`✅ Filtered pathways: ${Object.keys(filteredPathways).length} nodes have valid pathways`);
        console.log(`📊 Summary: ${validPathways}/${totalPathways} total pathways passed criteria`);
        return filteredPathways;
    }

    /**
     * Update grid view
     */
    updateGridView() {
        console.log("🔍 DEBUG: Starting degradation product identification...");

        // Get pathway filtering criteria from UI
        const minCount = parseInt(document.getElementById('gridMinCount').value) || 1;
        const maxBarrier = parseFloat(document.getElementById('gridMaxBarrier').value) || 30;

        // Get energy threshold for pathway filtering
        const energyThreshold = parseFloat(document.getElementById('gridEnergyThreshold').value) || 0;

        // Get other filtering criteria
        const neutralOnly = document.getElementById('neutralOnlyCheckbox').checked;
        const minWeight = document.getElementById('minMolWeight').value;
        const maxWeight = document.getElementById('maxMolWeight').value;
        const gridContainer = document.getElementById('gridContainer');

        console.log("🎛️ Current filter settings:");
        console.log(`  - Energy threshold: ${energyThreshold} kcal/mol`);
        console.log(`  - Min count: ${minCount}`);
        console.log(`  - Max barrier: ${maxBarrier} kcal/mol`);
        console.log(`  - Neutral only: ${neutralOnly}`);
        console.log(`  - Min weight: ${minWeight || 'none'}`);
        console.log(`  - Max weight: ${maxWeight || 'none'}`);

        // Get current visible nodes from the network
        const nodes = this.network.body.data.nodes;
        const currentNodes = nodes.get();

        console.log(`� Total visible nodes: ${currentNodes.length}`);
        console.log(`🎯 Sample visible nodes: ${currentNodes.slice(0, 10).map(n => n.id)}`);

        // Check specific nodes from the presentation
        const knownDegradationNodes = [961, 639, 1045, 1311, 986, 710];
        console.log("🧪 Checking known degradation products:");

        knownDegradationNodes.forEach(nodeId => {
            const isVisible = currentNodes.some(node => parseInt(node.id) === nodeId);
            const hasPathways = this.networkData.precomputedPathways && this.networkData.precomputedPathways[nodeId];
            const pathwayCount = hasPathways ? this.networkData.precomputedPathways[nodeId].length : 0;
            const hasCharge = this.networkData.molChargeMap ? this.networkData.molChargeMap[nodeId] : 'unknown';
            const hasSvg = this.networkData.svgMap ? !!this.networkData.svgMap[nodeId] : false;

            console.log(`  Node ${nodeId}:`);
            console.log(`    - Visible: ${isVisible}`);
            console.log(`    - Has pathways: ${!!hasPathways}`);
            console.log(`    - Pathway count: ${pathwayCount}`);
            console.log(`    - Charge: ${hasCharge}`);
            console.log(`    - Has SVG: ${hasSvg}`);

            if (hasPathways && pathwayCount > 0) {
                const firstPathway = this.networkData.precomputedPathways[nodeId][0];
                console.log(`    - First pathway structure:`, Object.keys(firstPathway));

                if (firstPathway.profile && firstPathway.profile.length > 0) {
                    const finalStep = firstPathway.profile[firstPathway.profile.length - 1];
                    const finalEnergy = finalStep.energy;
                    console.log(`    - Final energy: ${finalEnergy} kcal/mol (type: ${typeof finalEnergy})`);
                    console.log(`    - Profile length: ${firstPathway.profile.length}`);
                    console.log(`    - Final step structure:`, Object.keys(finalStep));
                    console.log(`    - Energy threshold check: ${finalEnergy} <= ${energyThreshold} = ${finalEnergy <= energyThreshold}`);

                    if (firstPathway.barrier_height !== undefined) {
                        console.log(`    - Barrier height: ${firstPathway.barrier_height} kcal/mol`);
                    }
                    if (firstPathway.count !== undefined) {
                        console.log(`    - Count: ${firstPathway.count}`);
                    }
                } else {
                    console.log(`    - Profile missing or empty`);
                }
            }
        });

        // Get all pathways that meet count and barrier criteria
        const filteredPathways = this.filterPathways(minCount, maxBarrier);
        console.log(`🔍 Found pathways for ${Object.keys(filteredPathways).length} nodes (count≥${minCount}, barrier≤${maxBarrier})`);

        // Debug pathway filtering for known nodes
        console.log("🔍 Pathway filtering results for known nodes:");
        knownDegradationNodes.forEach(nodeId => {
            const hasFilteredPathways = filteredPathways[nodeId];
            console.log(`  Node ${nodeId}: ${hasFilteredPathways ? hasFilteredPathways.length + ' pathways' : 'no pathways'} after filtering`);
        });

        // Filter visible nodes based on their pathway energy
        const validDegradationProducts = currentNodes.filter(node => {
            const nodeId = parseInt(node.id);

            // Skip substrate
            if (nodeId === parseInt(this.networkData.substrateId)) {
                return false;
            }

            // Must have SVG representation
            if (!this.networkData.svgMap[nodeId]) {
                return false;
            }

            // Must have pathways that meet criteria
            if (!filteredPathways[nodeId]) {
                return false;
            }

            // Check if any pathway has final energy below threshold
            const pathways = filteredPathways[nodeId];
            const hasFavorablePathway = pathways.some(pathway => {
                if (pathway.profile && pathway.profile.length > 0) {
                    const finalEnergy = pathway.profile[pathway.profile.length - 1].energy || 0;
                    return finalEnergy <= energyThreshold;
                }
                return false;
            });

            if (!hasFavorablePathway) {
                return false;
            }

            // Apply additional UI filters
            const charge = this.networkData.molChargeMap[nodeId] || 0;
            const weight = this.networkData.molWeightMap[nodeId] || 0;

            const chargeCheck = !neutralOnly || charge === 0;
            const minWeightCheck = !minWeight || weight >= parseFloat(minWeight);
            const maxWeightCheck = !maxWeight || weight <= parseFloat(maxWeight);

            return chargeCheck && minWeightCheck && maxWeightCheck;
        })
            .sort((a, b) => {
                // Sort by lowest pathway final energy
                const getLowestFinalEnergy = (nodeId) => {
                    const pathways = filteredPathways[nodeId] || [];
                    let minFinalEnergy = Infinity;
                    for (const pathway of pathways) {
                        if (pathway.profile && pathway.profile.length > 0) {
                            const finalEnergy = pathway.profile[pathway.profile.length - 1].energy || 0;
                            minFinalEnergy = Math.min(minFinalEnergy, finalEnergy);
                        }
                    }
                    return minFinalEnergy === Infinity ? 0 : minFinalEnergy;
                };

                return getLowestFinalEnergy(a.id) - getLowestFinalEnergy(b.id);
            });

        console.log(`✅ Final result: ${validDegradationProducts.length} degradation products meet all criteria`);

        // Debug which known nodes made it through
        console.log("🎯 Final results for known degradation products:");
        knownDegradationNodes.forEach(nodeId => {
            const madeIt = validDegradationProducts.some(node => parseInt(node.id) === nodeId);
            console.log(`  Node ${nodeId}: ${madeIt ? 'PASSED' : 'FILTERED OUT'}`);
        });

        // Clear existing content
        gridContainer.innerHTML = '';

        if (validDegradationProducts.length === 0) {
            let filterDesc = `visible nodes with pathways (count≥${minCount}, barrier≤${maxBarrier} kcal/mol, final energy≤${energyThreshold} kcal/mol)`;
            if (neutralOnly) filterDesc += ', neutral only';
            if (minWeight) filterDesc += `, MW≥${minWeight} Da`;
            if (maxWeight) filterDesc += `, MW≤${maxWeight} Da`;

            gridContainer.innerHTML = `
                <div style="grid-column: 1 / -1; text-align: center; color: #999; font-size: 16px; padding: 50px;">
                    No ${filterDesc} found in the current view.
                </div>
            `;
            return;
        }

        // Create grid items for degradation products
        validDegradationProducts.forEach(node => {
            const gridItem = this.createGridItem(node);
            gridContainer.appendChild(gridItem);
        });

        // Update modal title with count
        const modalTitle = document.querySelector('#gridViewModal h2');
        modalTitle.textContent = `Degradation Products (${validDegradationProducts.length} molecule${validDegradationProducts.length !== 1 ? 's' : ''})`;
    }

    /**
     * Create grid item for molecule
     */
    createGridItem(node) {
        const nodeId = node.id;
        const svg = this.networkData.svgMap[nodeId];
        const relEnergy = this.networkData.relEnergyMap[nodeId] || 0;
        const absEnergy = this.networkData.molEnergyMap[nodeId] || 0;
        const charge = this.networkData.molChargeMap[nodeId] || 0;
        const weight = this.networkData.molWeightMap[nodeId] || 0;

        const gridItem = document.createElement('div');
        gridItem.style.cssText = `
            border: 1px solid #ddd;
            border-radius: 8px;
            overflow: hidden;
            background: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            transition: transform 0.2s, box-shadow 0.2s;
            cursor: pointer;
            display: flex;
            flex-direction: column;
            height: 320px;
        `;

        // Add hover effects
        gridItem.addEventListener('mouseenter', function () {
            this.style.transform = 'translateY(-2px)';
            this.style.boxShadow = '0 4px 12px rgba(0,0,0,0.15)';
        });
        gridItem.addEventListener('mouseleave', function () {
            this.style.transform = 'translateY(0)';
            this.style.boxShadow = '0 2px 4px rgba(0,0,0,0.1)';
        });

        // Click to show full modal
        gridItem.addEventListener('click', () => {
            this.hideGridModal();
            this.showMoleculeModal(nodeId, true);
        });

        const headerDiv = document.createElement('div');
        headerDiv.style.cssText = `
            padding: 12px;
            background: #f8f9fa;
            border-bottom: 1px solid #ddd;
            text-align: center;
            flex-shrink: 0;
        `;

        const title = nodeId == this.networkData.substrateId ? `Substrate (Node ${nodeId})` : `Node ${nodeId}`;
        const chargeText = charge === 0 ? '' : ` (charge: ${charge > 0 ? '+' : ''}${charge})`;
        // Show absolute energy (G) for individual molecules in grid view
        headerDiv.innerHTML = `
            <div style="font-size: 14px; font-weight: bold; color: #333; margin-bottom: 4px;">${title}${chargeText}</div>
            <div style="font-size: 12px; color: #666;">
                <strong>G:</strong> ${absEnergy.toFixed(6)} hartree
            </div>
            <div style="font-size: 11px; color: #888;">
                ∆${this.networkData.energyType}: ${relEnergy.toFixed(2)} kcal/mol
            </div>
            <div style="font-size: 11px; color: #999;">
                MW: ${weight} Da
            </div>
        `;

        const svgDiv = document.createElement('div');
        svgDiv.style.cssText = `
            flex: 1;
            display: flex;
            align-items: center;
            justify-content: center;
            padding: 15px;
            background: white;
            overflow: hidden;
        `;
        svgDiv.innerHTML = svg;

        // Scale SVG to fit container
        const svgElement = svgDiv.querySelector('svg');
        if (svgElement) {
            svgElement.style.cssText = `
                max-width: 100%;
                max-height: 100%;
                width: auto;
                height: auto;
                object-fit: contain;
                pointer-events: none;
            `;
        }

        gridItem.appendChild(headerDiv);
        gridItem.appendChild(svgDiv);

        return gridItem;
    }

    /**
     * Search for a node by ID
     */
    searchNode() {
        const searchInput = document.getElementById("nodeSearch");
        const searchResult = document.getElementById("searchResult");
        const nodeId = parseInt(searchInput.value);

        if (isNaN(nodeId)) {
            searchResult.textContent = "Please enter a valid node ID";
            searchResult.style.color = "#e74c3c";
            return;
        }

        const nodes = this.network.body.data.nodes;
        const currentNodes = nodes.get();
        const currentNodeIds = new Set(currentNodes.map(n => n.id));

        if (!currentNodeIds.has(nodeId)) {
            // Check if node exists in original data but is filtered out
            const nodeExistsInOriginal = this.network.originalNodeData &&
                this.network.originalNodeData.some(n => n.id === nodeId);
            if (nodeExistsInOriginal) {
                searchResult.textContent = "Node exists but is currently filtered out";
                searchResult.style.color = "#f39c12";
            } else {
                searchResult.textContent = "Node not found";
                searchResult.style.color = "#e74c3c";
            }
            return;
        }

        // Node found and visible - select it
        this.clearSelection();
        this.selectedNodes.add(nodeId);
        this.updateNodeColor(nodeId, CONFIG.SELECTION_COLOR);
        this.lastSelectedNode = nodeId;
        this.updateSelectedEdges();
        this.updateSelectionDisplay();

        // Focus the network on the selected node
        this.network.focus(nodeId, {
            scale: 1.5,
            animation: {
                duration: CONFIG.ANIMATION.FOCUS_DURATION,
                easingFunction: "easeInOutQuad"
            }
        });

        searchResult.textContent = "";
        searchInput.value = "";
    }

    /**
     * Update selection display
     */
    updateSelectionDisplay() {
        const moleculeContainer = document.getElementById("moleculeContainer");

        if (this.selectedNodes.size === 0) {
            moleculeContainer.innerHTML = `
                <div id="noMoleculesMessage" style="text-align: center; color: #999; margin-top: 50px; font-size: 14px;">
                    No molecules selected.<br>
                    <span style="font-size: 12px;">Click on nodes to view their structures here.</span>
                </div>
            `;
            this.hideReactionChart();
        } else {
            // Use the ordered list instead of sorting by ID
            const nodeArray = [...this.selectedNodesInOrder];
            this.updateMoleculeViewer(nodeArray);
            this.updateReactionChart(nodeArray);
        }
    }

    /**
     * Update molecule viewer
     */
    updateMoleculeViewer(nodeIds) {
        const moleculeContainer = document.getElementById("moleculeContainer");
        moleculeContainer.innerHTML = '';

        nodeIds.forEach((nodeId) => {
            const svg = this.networkData.svgMap[nodeId];
            const relEnergy = this.networkData.relEnergyMap[nodeId] || 0;

            if (svg) {
                const moleculeDiv = this.createMoleculeViewerItem(nodeId, svg, relEnergy);
                moleculeContainer.appendChild(moleculeDiv);
            } else {
                const errorDiv = document.createElement('div');
                errorDiv.style.cssText = `
                    margin-bottom: 10px;
                    padding: 10px;
                    background: #f8d7da;
                    color: #721c24;
                    border: 1px solid #f5c6cb;
                    border-radius: 4px;
                    font-size: 12px;
                `;
                errorDiv.textContent = `No structure available for Node ${nodeId}`;
                moleculeContainer.appendChild(errorDiv);
            }
        });
    }

    /**
     * Create molecule viewer item
     */
    createMoleculeViewerItem(nodeId, svg, relEnergy) {
        const moleculeDiv = document.createElement('div');
        moleculeDiv.style.cssText = `
            margin-bottom: 20px;
            border: 1px solid #ddd;
            border-radius: 8px;
            overflow: hidden;
            background: white;
        `;

        const headerDiv = document.createElement('div');
        headerDiv.style.cssText = `
            padding: 8px 12px;
            background: #f8f9fa;
            border-bottom: 1px solid #ddd;
            font-size: 14px;
            font-weight: bold;
            display: flex;
            justify-content: space-between;
            align-items: center;
        `;

        const titleSpan = document.createElement('span');
        titleSpan.textContent = nodeId == this.networkData.substrateId ? `Substrate (${nodeId})` : `Node ${nodeId}`;

        const energyContainer = document.createElement('div');
        energyContainer.style.cssText = `
            font-size: 12px;
            color: #666;
            font-weight: normal;
            text-align: right;
        `;

        const relEnergySpan = document.createElement('div');
        relEnergySpan.textContent = `∆${this.networkData.energyType}: ${relEnergy.toFixed(2)} kcal/mol`;

        energyContainer.appendChild(relEnergySpan);

        const removeBtn = document.createElement('button');
        removeBtn.textContent = '×';
        removeBtn.style.cssText = `
            background: none;
            border: none;
            font-size: 16px;
            cursor: pointer;
            color: #999;
            padding: 0;
            width: 20px;
            height: 20px;
            display: flex;
            align-items: center;
            justify-content: center;
        `;
        removeBtn.onclick = () => {
            this.toggleNodeSelection(nodeId, true, false);
        };

        headerDiv.appendChild(titleSpan);
        headerDiv.appendChild(energyContainer);
        headerDiv.appendChild(removeBtn);

        const svgDiv = document.createElement('div');
        svgDiv.style.cssText = `
            padding: 10px;
            display: flex;
            justify-content: center;
            align-items: center;
            background: white;
            cursor: pointer;
            transition: background-color 0.2s;
        `;
        svgDiv.innerHTML = svg;

        // Add hover effect and click handler
        svgDiv.addEventListener('mouseenter', function () {
            this.style.backgroundColor = '#f8f9fa';
        });
        svgDiv.addEventListener('mouseleave', function () {
            this.style.backgroundColor = 'white';
        });
        svgDiv.addEventListener('click', () => {
            this.showMoleculeModal(nodeId);
        });

        // Scale SVG to fit container
        const svgElement = svgDiv.querySelector('svg');
        if (svgElement) {
            svgElement.style.cssText = 'max-width: 100%; height: auto; max-height: 300px; pointer-events: none;';
        }

        moleculeDiv.appendChild(headerDiv);
        moleculeDiv.appendChild(svgDiv);

        return moleculeDiv;
    }

    /**
     * Helper function to update node color
     */
    updateNodeColor(nodeId, color) {
        const nodes = this.network.body.data.nodes;
        nodes.update({ id: nodeId, color: color });
    }

    /**
     * Helper function to update edge color
     */
    updateEdgeColor(edgeId, color) {
        const edges = this.network.body.data.edges;
        edges.update({ id: edgeId, color: color });
    }

    /**
     * Clear all selections
     */
    clearSelection() {
        const nodes = this.network.body.data.nodes;
        const currentNodes = nodes.get();
        const currentNodeIds = new Set(currentNodes.map(n => n.id));

        this.selectedNodes.forEach(nodeId => {
            if (currentNodeIds.has(nodeId)) {
                const originalColor = this.originalNodeColors.get(nodeId);
                this.updateNodeColor(nodeId, originalColor);
            }
        });

        // Restore original edge colors
        this.highlightedEdges.forEach(edgeId => {
            const originalColor = this.originalEdgeColors.get(edgeId);
            if (originalColor) {
                this.updateEdgeColor(edgeId, originalColor);
            }
        });
        this.highlightedEdges.clear();

        this.selectedNodes.clear();
        this.selectedNodesInOrder = []; // Clear ordered list
        this.lastSelectedNode = null;
        this.updateSelectionDisplay();

        // Clear pyvis's internal selection
        this.network.unselectAll();
    }

    /**
     * Find shortest path between two nodes
     * Only considers nodes and edges that are currently visible in the network
     */
    findShortestPath(startNode, endNode, visibleNodeIds) {
        if (startNode === endNode) {
            return [startNode];
        }

        // Get the currently visible edges from the network
        const edges = this.network.body.data.edges;
        const visibleEdges = edges.get();

        // Build a graph of only visible connections
        const visibleGraph = new Map();
        visibleNodeIds.forEach(nodeId => {
            visibleGraph.set(nodeId, []);
        });

        // Add edges only if they're visible in the current filtered view
        visibleEdges.forEach(edge => {
            const fromId = edge.from;
            const toId = edge.to;

            // Add bidirectional connections since path can be traversed in either direction
            if (visibleGraph.has(fromId)) {
                visibleGraph.get(fromId).push(toId);
            }
            if (visibleGraph.has(toId)) {
                visibleGraph.get(toId).push(fromId);
            }
        });

        // Standard BFS with queue
        const queue = [[startNode]];
        const visited = new Set([startNode]);

        while (queue.length > 0) {
            const path = queue.shift();
            const currentNode = path[path.length - 1];

            // Get visible neighbors from our filtered graph
            const neighbors = visibleGraph.get(currentNode) || [];

            for (const neighbor of neighbors) {
                if (neighbor === endNode) {
                    return [...path, neighbor];
                }

                if (!visited.has(neighbor)) {
                    visited.add(neighbor);
                    queue.push([...path, neighbor]);
                }
            }
        }

        return [];
    }

    /**
     * Toggle node selection
     */
    toggleNodeSelection(nodeId, ctrlKey, shiftKey) {
        const nodes = this.network.body.data.nodes;
        const currentNodes = nodes.get();
        const currentNodeIds = new Set(currentNodes.map(n => n.id));

        if (!currentNodeIds.has(nodeId)) {
            return;
        }

        if (shiftKey && this.lastSelectedNode !== null) {
            // Range selection mode
            const pathNodes = this.findShortestPath(this.lastSelectedNode, nodeId, currentNodeIds);

            if (pathNodes.length > 0) {
                // Start from index 1 to avoid duplicate of lastSelectedNode
                for (let i = 1; i < pathNodes.length; i++) {
                    const id = pathNodes[i];
                    if (!this.selectedNodes.has(id)) {
                        this.selectedNodes.add(id);
                        this.selectedNodesInOrder.push(id); // Add to ordered list
                        this.updateNodeColor(id, CONFIG.SELECTION_COLOR);
                    }
                }
            } else {
                if (!this.selectedNodes.has(nodeId)) {
                    this.selectedNodes.add(nodeId);
                    this.selectedNodesInOrder.push(nodeId); // Add to ordered list
                    this.updateNodeColor(nodeId, CONFIG.SELECTION_COLOR);
                }
            }

            this.updateSelectedEdges();
            this.lastSelectedNode = nodeId;
        } else if (ctrlKey) {
            // Multi-selection mode
            if (this.selectedNodes.has(nodeId)) {
                // Remove from selection
                this.selectedNodes.delete(nodeId);

                // Remove from ordered list
                const index = this.selectedNodesInOrder.indexOf(nodeId);
                if (index !== -1) {
                    this.selectedNodesInOrder.splice(index, 1);
                }

                const originalColor = this.originalNodeColors.get(nodeId);
                this.updateNodeColor(nodeId, originalColor);

                if (this.lastSelectedNode === nodeId) {
                    this.lastSelectedNode = this.selectedNodesInOrder.length > 0 ?
                        this.selectedNodesInOrder[this.selectedNodesInOrder.length - 1] : null;
                }
            } else {
                // Add to selection
                this.selectedNodes.add(nodeId);
                this.selectedNodesInOrder.push(nodeId);
                this.updateNodeColor(nodeId, CONFIG.SELECTION_COLOR);
                this.lastSelectedNode = nodeId;
            }

            this.updateSelectedEdges();
        } else {
            // Single selection mode
            this.clearSelection();
            this.selectedNodes.add(nodeId);
            this.selectedNodesInOrder = [nodeId]; // Reset ordered list with single node
            this.updateNodeColor(nodeId, CONFIG.SELECTION_COLOR);
            this.lastSelectedNode = nodeId;
        }

        this.updateSelectedEdges();
        this.updateSelectionDisplay();

        // Check if selected nodes form a path and update reaction profile accordingly
        this.checkAndUpdateReactionProfile();
    }

    /**
     * Update edge colors based on selection
     */
    updateSelectedEdges() {
        const edges = this.network.body.data.edges;
        const allEdges = edges.get();

        // Restore all previously highlighted edges
        this.highlightedEdges.forEach(edgeId => {
            const originalColor = this.originalEdgeColors.get(edgeId);
            if (originalColor) {
                this.updateEdgeColor(edgeId, originalColor);
            }
        });
        this.highlightedEdges.clear();

        // Highlight edges that connect selected nodes
        allEdges.forEach(edge => {
            const isSelected = this.selectedNodes.has(edge.from) && this.selectedNodes.has(edge.to);
            if (isSelected) {
                this.updateEdgeColor(edge.id, CONFIG.SELECTION_COLOR);
                this.highlightedEdges.add(edge.id);
            }
        });
    }

    /**
     * Check if selected nodes form a pathway and update reaction profile accordingly
     * This method detects when nodes are selected in a sequence that could represent a reaction path
     */
    checkAndUpdateReactionProfile() {
        // Only proceed if we have multiple nodes selected in order
        if (this.selectedNodesInOrder.length < 2) {
            return;
        }

        // Check if we have a pathway viewer available
        if (!window.pathwayViewer) {
            return;
        }

        // Get the selected path
        const selectedPath = [...this.selectedNodesInOrder];

        console.log(`🔍 Checking reaction profile for selected path: [${selectedPath.join(' → ')}]`);

        // Try to find and display the reaction profile for this path
        try {
            window.pathwayViewer.showReactionProfileForPath(selectedPath);
        } catch (error) {
            console.warn('Error updating reaction profile for selected path:', error);
        }
    }

    /**
     * Hide reaction chart
     */
    hideReactionChart() {
        document.getElementById('reactionProfileContainer').style.display = 'none';
    }

    /**
     * Update reaction chart
     */
    updateReactionChart(nodeArray) {
        if (nodeArray.length === 0) {
            this.hideReactionChart();
            return;
        }

        document.getElementById('reactionProfileContainer').style.display = 'block';

        if (!this.reactionChart) {
            this.initializeReactionChart();
        }

        const labels = [];
        const energies = [];

        nodeArray.forEach(nodeId => {
            const energy = this.networkData.relEnergyMap[nodeId] || 0;
            const label = nodeId == this.networkData.substrateId ? `S` : `${nodeId}`;
            labels.push(label);
            energies.push(energy);
        });

        this.reactionChart.data.labels = labels;
        this.reactionChart.data.datasets[0].data = energies;

        const hasSubstrate = nodeArray.includes(this.networkData.substrateId);
        if (hasSubstrate) {
            this.reactionChart.data.datasets[0].borderColor = CONFIG.CHART_COLORS.PRIMARY;
            this.reactionChart.data.datasets[0].backgroundColor = CONFIG.CHART_COLORS.PRIMARY_BACKGROUND;
            this.reactionChart.data.datasets[0].pointBackgroundColor = CONFIG.CHART_COLORS.PRIMARY;
        } else {
            this.reactionChart.data.datasets[0].borderColor = CONFIG.CHART_COLORS.SECONDARY;
            this.reactionChart.data.datasets[0].backgroundColor = CONFIG.CHART_COLORS.SECONDARY_BACKGROUND;
            this.reactionChart.data.datasets[0].pointBackgroundColor = CONFIG.CHART_COLORS.SECONDARY;
        }

        if (document.getElementById('reactionProfileChart').style.display !== 'none') {
            this.reactionChart.update(CONFIG.ANIMATION.CHART_UPDATE);
        }
    }

    /**
     * Initialize reaction chart
     */
    initializeReactionChart() {
        const ctx = document.getElementById('reactionProfileCanvas').getContext('2d');

        this.reactionChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels: [],
                datasets: [{
                    label: `∆${this.networkData.energyType} (kcal/mol)`,
                    data: [],
                    borderColor: CONFIG.CHART_COLORS.PRIMARY,
                    backgroundColor: CONFIG.CHART_COLORS.PRIMARY_BACKGROUND,
                    pointBackgroundColor: CONFIG.CHART_COLORS.PRIMARY,
                    pointBorderColor: '#fff',
                    pointBorderWidth: 2,
                    pointRadius: 6,
                    pointHoverRadius: 8,
                    fill: true,
                    tension: 0.1
                }]
            },
            options: {
                responsive: true,
                maintainAspectRatio: false,
                scales: {
                    y: {
                        beginAtZero: false,
                        title: {
                            display: true,
                            text: `∆${this.networkData.energyType} (kcal/mol)`,
                            font: { size: 12 }
                        },
                        grid: { color: 'rgba(0,0,0,0.1)' }
                    },
                    x: {
                        title: {
                            display: true,
                            text: 'Node',
                            font: { size: 12 }
                        },
                        grid: { color: 'rgba(0,0,0,0.1)' }
                    }
                },
                plugins: {
                    legend: { display: false },
                    tooltip: {
                        callbacks: {
                            title: (context) => {
                                const nodeId = context[0].label;
                                return nodeId == this.networkData.substrateId ? `Substrate (${nodeId})` : `Node ${nodeId}`;
                            },
                            label: (context) => {
                                return `∆${this.networkData.energyType}: ${context.parsed.y.toFixed(2)} kcal/mol`;
                            }
                        }
                    }
                },
                interaction: { intersect: false, mode: 'index' }
            }
        });
    }

    /**
     * Toggle reaction profile visibility
     */
    toggleReactionProfile() {
        const chartDiv = document.getElementById('reactionProfileChart');
        const toggleSpan = document.getElementById('reactionProfileToggle');

        if (chartDiv.style.display === 'none') {
            chartDiv.style.display = 'block';
            toggleSpan.textContent = '▲';
            if (this.reactionChart && this.selectedNodes.size > 0) {
                this.reactionChart.update('none');
            }
        } else {
            chartDiv.style.display = 'none';
            toggleSpan.textContent = '▼';
        }
    }

    /**
     * Show indices modal with visible nodes and edges
     */
    /**
     * Show indices modal
     * @param {string} [type='nodes'] - The type of indices to focus on, 'nodes' or 'edges'
     * @param {string} [edgeType=null] - When specified, only show edges of this type
     */
    showIndicesModal(type = 'nodes', edgeType = null) {
        const modal = document.getElementById('indicesModal');
        const nodeIndicesList = document.getElementById('nodeIndicesList');
        const edgeIndicesList = document.getElementById('edgeIndicesList');
        const nodeIndicesSection = document.querySelector('.indices-section:nth-child(1)');
        const edgeIndicesSection = document.querySelector('.indices-section:nth-child(2)');

        // Get visible nodes
        const visibleNodes = this.network.body.data.nodes.get({
            filter: function (node) {
                return node.hidden !== true;
            }
        });

        // Sort nodes by ID
        visibleNodes.sort((a, b) => parseInt(a.id) - parseInt(b.id));

        // Create array of node IDs for both display and copying
        const nodeIds = visibleNodes.map(node => node.id);

        // Display node indices
        nodeIndicesList.innerHTML = nodeIds.length > 0 ?
            nodeIds.map(id => `<div class="index-item">${id}</div>`).join('') :
            '<span class="no-indices">No nodes found.</span>';

        // Get visible edges
        let visibleEdges = this.network.body.data.edges.get({
            filter: function (edge) {
                return edge.hidden !== true;
            }
        });

        // Filter edges by type if specified
        if (edgeType) {
            visibleEdges = visibleEdges.filter(edge => {
                const edgeDataInfo = window.edgeData ? window.edgeData.find(e => e.id === edge.id) : null;
                const currentEdgeType = edgeDataInfo ? edgeDataInfo.type : 'unknown';
                return currentEdgeType === edgeType;
            });
        }

        // Create edge indices in format {start_index}-{end_index}
        const edgeIndices = visibleEdges.map(edge => `${edge.from}-${edge.to}`);
        edgeIndices.sort(); // Sort edges alphabetically

        // Display edge indices
        edgeIndicesList.innerHTML = edgeIndices.length > 0 ?
            edgeIndices.map(pair => `<div class="index-item">${pair}</div>`).join('') :
            '<span class="no-indices">No edges found.</span>';        // Set up copy buttons
        const copyNodeIndicesBtn = document.getElementById('copyNodeIndices');
        const copyEdgeIndicesBtn = document.getElementById('copyEdgeIndices');

        // Remove any existing listeners (to prevent duplicates)
        copyNodeIndicesBtn.replaceWith(copyNodeIndicesBtn.cloneNode(true));
        copyEdgeIndicesBtn.replaceWith(copyEdgeIndicesBtn.cloneNode(true));

        // Get the new button references
        const newCopyNodeBtn = document.getElementById('copyNodeIndices');
        const newCopyEdgeBtn = document.getElementById('copyEdgeIndices');

        // Add click event listeners for copying
        newCopyNodeBtn.addEventListener('click', () => {
            this.copyToClipboard(nodeIds.join(' '), newCopyNodeBtn);
        });

        newCopyEdgeBtn.addEventListener('click', () => {
            this.copyToClipboard(edgeIndices.join(' '), newCopyEdgeBtn);
        });

        // Set modal title based on what was clicked
        let modalTitle = type === 'nodes' ? 'Node Indices' : 'Edge Indices';
        if (edgeType && type === 'edges') {
            const capitalizedType = edgeType.charAt(0).toUpperCase() + edgeType.slice(1);
            modalTitle = `${capitalizedType} Edge Indices`;
        }
        document.getElementById('indicesModalTitle').textContent = modalTitle;

        // Show only the relevant section
        if (type === 'nodes') {
            // Show only nodes section
            nodeIndicesSection.style.display = 'block';
            edgeIndicesSection.style.display = 'none';
        } else {
            // Show only edges section
            nodeIndicesSection.style.display = 'none';
            edgeIndicesSection.style.display = 'block';
        }

        // Show the modal
        modal.style.display = 'block';
    }

    /**
     * Hide indices modal
     */
    hideIndicesModal() {
        const modal = document.getElementById('indicesModal');
        modal.style.display = 'none';

        // Reset display properties for next time
        const nodeIndicesSection = document.querySelector('.indices-section:nth-child(1)');
        const edgeIndicesSection = document.querySelector('.indices-section:nth-child(2)');

        // Set both to 'block' as default for next opening
        nodeIndicesSection.style.display = 'block';
        edgeIndicesSection.style.display = 'block';
    }

    /**
     * Hide the stereoisomer modal
     * This is a placeholder implementation for API compatibility
     */
    hideStereoisomerModal() {
        console.log("Stereoisomer modal hide functionality has been disabled for reimplementation");
        // Placeholder for future implementation
    }

    /**
     * Copy text to clipboard and show feedback
     */
    copyToClipboard(text, button) {
        // Create a temporary textarea element to hold the text
        const textarea = document.createElement('textarea');
        textarea.value = text;
        textarea.setAttribute('readonly', ''); // Make it readonly to avoid focus
        textarea.style.position = 'absolute';
        textarea.style.left = '-9999px'; // Move off-screen
        document.body.appendChild(textarea);

        // Select the text and copy it
        textarea.select();
        document.execCommand('copy');

        // Remove the temporary element
        document.body.removeChild(textarea);

        // Show visual feedback
        button.classList.add('copy-success');

        // Revert the button style after a short delay
        setTimeout(() => {
            button.classList.remove('copy-success');
        }, 1500);
    }
}
