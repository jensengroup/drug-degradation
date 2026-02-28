/**
 * Main entry point - Initialize the application
 */
(function() {
    // Initialize data
    Data.init(NETWORK_DATA);

    // Initialize network graph
    NetworkGraph.init('network-container');
    NetworkGraph.update(Data.nodes, Data.edges, Data.junctionNodes);

    // Initialize UI
    UI.init();

    console.log('Reaction Network Visualization initialized');
    console.log(`Loaded ${Data.nodes.length} nodes, ${Data.edges.length} edges`);
    if (Data.junctionNodes.length > 0) {
        console.log(`  Including ${Data.junctionNodes.length} hyperedge junctions`);
    }
})();
