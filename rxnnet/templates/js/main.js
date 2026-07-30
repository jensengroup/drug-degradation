/**
 * Main entry point - Initialize the application
 */
(function() {
    // Initialize data
    Data.init(NETWORK_DATA);

    // Initialize network graph
    NetworkGraph.init('network-container');

    // Initialize UI
    UI.init();

    // Apply the sliders' default values to the initial view instead of
    // showing the raw, unfiltered network.
    UI.applyFilters();

    console.log('Reaction Network Visualization initialized');
    console.log(`Loaded ${Data.nodes.length} nodes, ${Data.edges.length} edges`);
    if (Data.junctionNodes.length > 0) {
        console.log(`  Including ${Data.junctionNodes.length} hyperedge junctions`);
    }
})();
