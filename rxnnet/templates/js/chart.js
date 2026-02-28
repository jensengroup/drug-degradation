/**
 * Chart.js energy profile charts
 */
const ChartModule = {
    sidebarChart: null,
    modalChart: null,

    createSidebarChart(profile) {
        const ctx = document.getElementById('reactionProfileCanvas')?.getContext('2d');
        if (!ctx) return;

        if (this.sidebarChart) {
            this.sidebarChart.destroy();
        }

        const { labels, data, pointColors } = this.prepareChartData(profile);

        this.sidebarChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels,
                datasets: [{
                    label: '∆E (kcal/mol)',
                    data,
                    borderColor: CONFIG.CHART_COLORS.PRIMARY,
                    backgroundColor: CONFIG.CHART_COLORS.PRIMARY_BACKGROUND,
                    pointBackgroundColor: pointColors,
                    pointBorderColor: '#fff',
                    pointBorderWidth: 2,
                    pointRadius: 5,
                    fill: true,
                    tension: 0.1,
                    borderWidth: 2
                }]
            },
            options: this.getChartOptions(false)
        });
    },

    createModalChart(profile) {
        const ctx = document.getElementById('modalReactionChart')?.getContext('2d');
        if (!ctx) return;

        if (this.modalChart) {
            this.modalChart.destroy();
        }

        const { labels, data, pointColors } = this.prepareChartData(profile);

        this.modalChart = new Chart(ctx, {
            type: 'line',
            data: {
                labels,
                datasets: [{
                    label: '∆E (kcal/mol)',
                    data,
                    borderColor: CONFIG.CHART_COLORS.PRIMARY,
                    backgroundColor: CONFIG.CHART_COLORS.PRIMARY_BACKGROUND,
                    pointBackgroundColor: pointColors,
                    pointBorderColor: '#fff',
                    pointBorderWidth: 3,
                    pointRadius: 8,
                    pointHoverRadius: 10,
                    fill: true,
                    tension: 0.1,
                    borderWidth: 3
                }]
            },
            options: this.getChartOptions(true)
        });
    },

    prepareChartData(profile) {
        const labels = profile.map(p => {
            if (p.type === 'transition_state') return 'TS';
            return p.node_id.toString();
        });

        const data = profile.map(p => p.energy);

        const pointColors = profile.map(p =>
            p.type === 'transition_state' ? CONFIG.CHART_COLORS.TS_POINT : CONFIG.CHART_COLORS.PRIMARY
        );

        return { labels, data, pointColors };
    },

    getChartOptions(isModal) {
        return {
            responsive: true,
            maintainAspectRatio: false,
            layout: {
                padding: isModal ? { bottom: 10 } : 5
            },
            scales: {
                y: {
                    beginAtZero: false,
                    title: {
                        display: isModal,
                        text: '∆E (kcal/mol)',
                        font: { size: isModal ? 14 : 10 }
                    },
                    grid: { color: 'rgba(0,0,0,0.1)' },
                    ticks: { font: { size: isModal ? 12 : 9 } }
                },
                x: {
                    title: {
                        display: isModal,
                        text: 'Node',
                        font: { size: isModal ? 14 : 10 }
                    },
                    grid: { color: 'rgba(0,0,0,0.1)' },
                    ticks: {
                        font: { size: isModal ? 12 : 9 },
                        color: '#4A90E2'
                    }
                }
            },
            plugins: {
                legend: { display: false },
                tooltip: {
                    titleFont: { size: isModal ? 14 : 11 },
                    bodyFont: { size: isModal ? 12 : 10 },
                    callbacks: {
                        title: (ctx) => {
                            const label = ctx[0].label;
                            if (label === 'TS') return 'Transition State';
                            const nodeId = parseInt(label);
                            return nodeId === Data.raw.substrateId ? `Substrate (${nodeId})` : `Node ${nodeId}`;
                        },
                        label: (ctx) => `∆E: ${ctx.parsed.y.toFixed(2)} kcal/mol`
                    }
                }
            },
            interaction: { intersect: false, mode: 'index' }
        };
    },

    destroy(chartType) {
        if (chartType === 'sidebar' && this.sidebarChart) {
            this.sidebarChart.destroy();
            this.sidebarChart = null;
        }
        if (chartType === 'modal' && this.modalChart) {
            this.modalChart.destroy();
            this.modalChart = null;
        }
    }
};
