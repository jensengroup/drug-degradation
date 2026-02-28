/**
 * Configuration constants for the network visualization
 */
const CONFIG = {
    // Default values
    DEFAULT_MAX_BARRIER_THRESHOLD: 50,

    // UI constants
    TOOLTIP_WIDTH_FRACTION: 0.35, // 35% of viewport width
    TOOLTIP_HEIGHT_FRACTION: 0.35, // 35% of viewport height
    TOOLTIP_PADDING: 10,
    TOOLTIP_TEXT_HEIGHT: 60,

    // Selection colors
    SELECTION_COLOR: "#DC143C",
    SELECTED_HIGHLIGHT_COLOR: "#FFD700",

    // Edge colors by reaction type
    EDGE_COLORS: {
        'tautomerization': '#55A69A',
        'tautomer': '#55A69A',
        'protonation': '#E2A334',
        'deprotonation': '#3C56E9',
        'reaction': '#808080',
        'default': '#00000084'
    },

    // Chart colors
    CHART_COLORS: {
        PRIMARY: '#4A90E2',
        PRIMARY_BACKGROUND: 'rgba(74, 144, 226, 0.1)',
        SECONDARY: '#8C4AE2',
        SECONDARY_BACKGROUND: 'rgba(140, 74, 226, 0.1)'
    },

    // Animation durations
    ANIMATION: {
        FOCUS_DURATION: 500,
        CHART_UPDATE: 'none' // No animation for responsiveness
    }
};
