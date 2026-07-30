/**
 * Configuration constants for the network visualization
 */
const CONFIG = {
    // Default filter values
    DEFAULT_COUNT_THRESHOLD: 1,
    DEFAULT_PATH_ENERGY_THRESHOLD: 50,

    // Selection colors
    SELECTION_COLOR: '#DC143C',
    SELECTED_HIGHLIGHT_COLOR: '#FFD700',

    // Node colors
    NODE_COLORS: {
        default: '#C2C2C2',
        substrate: '#4E4E4E',
        selected: '#FFD700',
        positive: '#27ae60',
        negative: '#e74c3c'
    },

    // Edge colors by reaction type
    EDGE_COLORS: {
        'tautomerization': '#55A69A',
        'tautomer': '#55A69A',
        'protonation': '#E2A334',
        'protonate': '#E2A334',
        'deprotonation': '#3C56E9',
        'deprotonate': '#3C56E9',
        'reaction': '#808080',
        'mtd-reaction': '#808080',
        'add': '#808080',
        'default': '#00000084'
    },

    // Chart colors
    CHART_COLORS: {
        PRIMARY: '#4A90E2',
        PRIMARY_BACKGROUND: 'rgba(74, 144, 226, 0.1)',
        SECONDARY: '#8C4AE2',
        SECONDARY_BACKGROUND: 'rgba(140, 74, 226, 0.1)',
        TS_POINT: '#e94560'
    },

    // Animation
    ANIMATION: {
        FOCUS_DURATION: 500
    }
};
