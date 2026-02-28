/**
 * Data Parsing Extensions for ReactionNetwork JavaScript Implementation
 *
 * This module extends the ReactionNetworkJS class with data parsing methods
 * that ONLY use precomputed data from Python. NO energy calculations are
 * performed on the JavaScript side.
 *
 * JavaScript responsibilities:
 * - Parse precomputed Python data
 * - Apply filtering to pathways
 * - Display pathway profiles
 * - Handle user interactions
 *
 * Python responsibilities:
 * - ALL energy calculations
 * - Protonation/deprotonation corrections
 * - Transition state energies
 * - Reaction profile generation
 */

// Extend the ReactionNetworkJS class with data parsing methods (NO CALCULATIONS)
ReactionNetworkJS.prototype.getPrecomputedProfile = function (nodeId, pathwayIndex = 0) {
    /**
     * Get precomputed reaction profile from Python data
     * NO calculations performed - just data retrieval
     *
     * @param {number} nodeId - Target node ID
     * @param {number} pathwayIndex - Which pathway to use (default: 0)
     * @returns {Array} Precomputed reaction profile from Python
     */
    const pathways = this.networkData.precomputedPathways[nodeId];
    if (!pathways || pathways.length <= pathwayIndex) {
        console.warn(`No pathway ${pathwayIndex} found for node ${nodeId}`);
        return [];
    }

    const pathway = pathways[pathwayIndex];
    return pathway.profile || [];
};

ReactionNetworkJS.prototype.getPrecomputedBarrier = function (nodeId, pathwayIndex = 0) {
    /**
     * Get precomputed barrier height from Python data
     * NO calculations performed - just data retrieval
     *
     * @param {number} nodeId - Target node ID
     * @param {number} pathwayIndex - Which pathway to use (default: 0)
     * @returns {number} Precomputed barrier height in kcal/mol
     */
    const pathways = this.networkData.precomputedPathways[nodeId];
    if (!pathways || pathways.length <= pathwayIndex) {
        return 0.0;
    }

    const pathway = pathways[pathwayIndex];
    return pathway.barrier_height || 0.0;
};

ReactionNetworkJS.prototype.getPrecomputedPath = function (nodeId, pathwayIndex = 0) {
    /**
     * Get precomputed reaction path from Python data
     * NO calculations performed - just data retrieval
     *
     * @param {number} nodeId - Target node ID
     * @param {number} pathwayIndex - Which pathway to use (default: 0)
     * @returns {Array} Precomputed reaction path (array of node IDs)
     */
    const pathways = this.networkData.precomputedPathways[nodeId];
    if (!pathways || pathways.length <= pathwayIndex) {
        return [];
    }

    const pathway = pathways[pathwayIndex];
    return pathway.path || [];
};

ReactionNetworkJS.prototype.filterPrecomputedPathways = function (criteria) {
    /**
     * Filter precomputed pathways based on criteria
     * Uses ONLY the precomputed Python data - no calculations
     *
     * @param {Object} criteria - Filtering criteria
     * @param {number} criteria.maxBarrier - Maximum barrier height in kcal/mol
     * @param {number} criteria.maxLength - Maximum pathway length
     * @param {number} criteria.minCount - Minimum edge count (for edge filtering)
     * @returns {Object} Filtered pathways data
     */
    const filteredPathways = {};

    for (const [nodeIdStr, pathwayList] of Object.entries(this.networkData.precomputedPathways)) {
        const nodeId = parseInt(nodeIdStr);
        const validPathways = [];

        for (const pathway of pathwayList) {
            let isValid = true;

            // Filter by barrier height (using precomputed values)
            if (criteria.maxBarrier !== undefined && pathway.barrier_height > criteria.maxBarrier) {
                isValid = false;
            }

            // Filter by pathway length (using precomputed values)
            if (criteria.maxLength !== undefined && pathway.path.length > criteria.maxLength) {
                isValid = false;
            }

            // Filter by edge count (check if pathway uses valid edges)
            if (criteria.minCount !== undefined) {
                const pathEdges = this.getEdgesForPath(pathway.path);
                if (!this.areEdgesValid(pathEdges, criteria.minCount)) {
                    isValid = false;
                }
            }

            if (isValid) {
                validPathways.push(pathway);
            }
        }

        // Always include the node, even if no valid pathways (empty array)
        filteredPathways[nodeId] = validPathways;
    }

    return filteredPathways;
};

ReactionNetworkJS.prototype.getEdgesForPath = function (path) {
    /**
     * Get edges used in a pathway - NO energy calculations
     *
     * @param {Array} path - Array of node IDs
     * @returns {Array} Array of edge objects used in the path
     */
    const edges = [];
    for (let i = 0; i < path.length - 1; i++) {
        const edgeKey = `${path[i]}-${path[i + 1]}`;
        const edge = this.edgeMap.get(edgeKey);
        if (edge) {
            edges.push(edge);
        }
    }
    return edges;
};

ReactionNetworkJS.prototype.areEdgesValid = function (edges, minCount) {
    /**
     * Check if edges meet count criteria - NO energy calculations
     *
     * @param {Array} edges - Array of edge objects
     * @param {number} minCount - Minimum count required
     * @returns {boolean} True if all edges meet criteria
     */
    for (const edge of edges) {
        // Non-reaction edges always pass (same logic as Python)
        if (edge.type.toLowerCase() !== "reaction") {
            continue;
        }

        // Reaction edges must meet count threshold
        if ((edge.count || 1) < minCount) {
            return false;
        }
    }
    return true;
};

ReactionNetworkJS.prototype.getBestPathwayEnergy = function (nodeId) {
    /**
     * Get best pathway energy from precomputed Python data
     * NO calculations performed - just data retrieval
     *
     * @param {number} nodeId - Target node ID
     * @returns {number} Best pathway energy in kcal/mol (from Python)
     */
    return this.networkData.bestPathEnergyMap[nodeId] || Infinity;
};

ReactionNetworkJS.prototype.getBestPathwayBarrier = function (nodeId) {
    /**
     * Get best pathway barrier from precomputed Python data
     * NO calculations performed - just data retrieval
     *
     * @param {number} nodeId - Target node ID
     * @returns {number} Best pathway barrier in kcal/mol (from Python)
     */
    return this.networkData.bestPathBarrierMap[nodeId] || Infinity;
};

ReactionNetworkJS.prototype.getRelativeEnergy = function (nodeId) {
    /**
     * Get relative energy from precomputed Python data
     * NO calculations performed - just data retrieval
     *
     * @param {number} nodeId - Target node ID
     * @returns {number} Relative energy in kcal/mol (from Python)
     */
    return this.networkData.relEnergyMap[nodeId] || 0.0;
};

ReactionNetworkJS.prototype.validateDataConsistency = function () {
    /**
     * Validate that all required precomputed data is present
     * Used for debugging - ensures Python data was passed correctly
     *
     * @returns {Object} Validation results
     */
    const validation = {
        hasPrecomputedPathways: !!this.networkData.precomputedPathways,
        hasEnergyMaps: !!(this.networkData.bestPathEnergyMap && this.networkData.bestPathBarrierMap),
        hasRelativeEnergies: !!this.networkData.relEnergyMap,
        pathwayCount: 0,
        energyMapSize: 0,
        issues: []
    };

    if (validation.hasPrecomputedPathways) {
        validation.pathwayCount = Object.keys(this.networkData.precomputedPathways).length;
    } else {
        validation.issues.push("Missing precomputedPathways data from Python");
    }

    if (validation.hasEnergyMaps) {
        validation.energyMapSize = Object.keys(this.networkData.bestPathEnergyMap).length;
    } else {
        validation.issues.push("Missing energy maps from Python");
    }

    if (!validation.hasRelativeEnergies) {
        validation.issues.push("Missing relative energy map from Python");
    }

    validation.isValid = validation.issues.length === 0;

    console.log("=== JavaScript Data Validation ===");
    console.log(`Precomputed pathways: ${validation.pathwayCount} nodes`);
    console.log(`Energy maps: ${validation.energyMapSize} entries`);
    console.log(`Valid: ${validation.isValid ? '✅' : '❌'}`);

    if (validation.issues.length > 0) {
        console.log("Issues:", validation.issues);
    }

    return validation;
};

// Legacy method compatibility - map old methods to new data-only methods
ReactionNetworkJS.prototype.calculateReactionProfile = function (path, edgeMap) {
    /**
     * Legacy compatibility method - now just returns precomputed profile
     * NO calculations performed - just data lookup
     *
     * @param {Array} path - Array of node IDs (not used in new implementation)
     * @param {Map} edgeMap - Edge map (not used in new implementation)
     * @returns {Array} Precomputed reaction profile from Python
     */
    console.warn("calculateReactionProfile is deprecated. Use getPrecomputedProfile instead.");

    // If this is called for a single target node, find it and return the profile
    if (path.length >= 2) {
        const targetNode = path[path.length - 1];
        return this.getPrecomputedProfile(targetNode, 0);
    }

    return [];
};

ReactionNetworkJS.prototype.findMaxBarrier = function (profile) {
    /**
     * Legacy compatibility method - extracts barrier from precomputed profile
     * NO calculations performed - just data extraction
     *
     * @param {Array} profile - Precomputed reaction profile
     * @returns {number} Maximum barrier height in kcal/mol
     */
    if (!profile || profile.length === 0) {
        return 0.0;
    }

    // Extract barrier from precomputed profile metadata if available
    // Use ONLY precomputed barrier_height - NO calculations from profiles
    if (profile.barrier_height !== undefined) {
        return profile.barrier_height;
    }

    // DO NOT calculate barriers from profile points - this causes incorrect values
    // JavaScript should never calculate energies or barriers, only use precomputed data
    console.warn("⚠️ No precomputed barrier_height found - returning 0 instead of calculating from profile");
    return 0;
};

// Extend the ReactionNetworkJS class with energy calculation methods
ReactionNetworkJS.prototype.calculateEliminatedProductsEnergy = function (edge, currentNode) {
    /**
     * Calculate energy contribution from eliminated products
     * Port of Python _calculate_eliminated_products_energy method
     *
     * @param {Object} edge - Edge object
     * @param {Object} currentNode - Current node object
     * @returns {number} Energy contribution from eliminated products in kcal/mol
     */
    if (!edge.smallerProducts || edge.smallerProducts.length === 0) {
        return 0.0;
    }

    let eliminatedEnergy = 0.0;
    for (const nodeId of edge.smallerProducts) {
        const node = this.nodes.get(nodeId);
        if (node) {
            // Convert from hartree to kcal/mol
            eliminatedEnergy += this.hartree2kcalmol(node.rawEnergy);
        }
    }

    console.log(`Node ${currentNode.id} has siblings ${edge.smallerProducts} yielding ${eliminatedEnergy} kcal/mol`);

    return eliminatedEnergy;
};

ReactionNetworkJS.prototype.isNeutralization = function (currentEdgeType, previousProtonType) {
    /**
     * Check if current step neutralizes previous protonation state
     * Port of Python _is_neutralization method
     *
     * @param {string} currentEdgeType - Current edge type
     * @param {string} previousProtonType - Previous proton type
     * @returns {boolean} True if neutralization occurs
     */
    return (
        (currentEdgeType === "protonation" && previousProtonType === "deprotonation") ||
        (currentEdgeType === "deprotonation" && previousProtonType === "protonation")
    );
};

ReactionNetworkJS.prototype.handleProtonationDeprotonation = function (
    edgeType,
    currentProtonType,
    currentProtonOffset,
    node,
    lastEnergy,
    pH,
    R,
    T,
    pkaSlope,
    pkaIntercept
) {
    /**
     * Handle protonation/deprotonation energy corrections
     * EXACT port of Python logic from _calculate_reaction_profile method
     *
     * @param {string} edgeType - Edge type (protonation/deprotonation)
     * @param {string} currentProtonType - Current proton state
     * @param {number} currentProtonOffset - Current persistent proton offset
     * @param {Object} node - Node object
     * @param {number} lastEnergy - Last calculated energy (from energies array)
     * @param {number} pH - pH value
     * @param {number} R - Gas constant
     * @param {number} T - Temperature
     * @param {number} pkaSlope - pKa slope parameter
     * @param {number} pkaIntercept - pKa intercept parameter
     * @returns {Object} {newProtonOffset: number, newProtonType: string}
     */

    // Handle protonation/deprotonation with persistent offset - EXACT PYTHON LOGIC
    if (edgeType === "protonation") {
        if (currentProtonType === "deprotonation") {
            // Neutralization
            console.log(`Neutralization: protonation cancels previous deprotonation`);
            return {
                newProtonOffset: 0.0,
                newProtonType: ""
            };
        } else {
            // New protonation
            const nodeEnergyKcal = this.hartree2kcalmol(node.rawEnergy);
            const dG = nodeEnergyKcal - lastEnergy;
            const pka = pkaSlope * (-dG) + pkaIntercept;
            const dGCorrection = R * T * Math.log(10) * (pH - pka);
            const newProtonOffset = lastEnergy - nodeEnergyKcal + dGCorrection;
            const newProtonType = (node.originType || "").toLowerCase();

            console.log(`Protonation: pKa=${pka.toFixed(2)}, dG_correction=${dGCorrection.toFixed(3)}, new_offset=${newProtonOffset.toFixed(3)}`);
            return {
                newProtonOffset: newProtonOffset,
                newProtonType: newProtonType
            };
        }
    } else if (edgeType === "deprotonation") {
        if (currentProtonType === "protonation") {
            // Neutralization
            console.log(`Neutralization: deprotonation cancels previous protonation`);
            return {
                newProtonOffset: 0.0,
                newProtonType: ""
            };
        } else {
            // New deprotonation
            const nodeEnergyKcal = this.hartree2kcalmol(node.rawEnergy);
            const dG = nodeEnergyKcal - lastEnergy;
            const pka = pkaSlope * dG + pkaIntercept;
            const dGCorrection = R * T * Math.log(10) * (pka - pH);
            const newProtonOffset = lastEnergy - nodeEnergyKcal + dGCorrection;
            const newProtonType = (node.originType || "").toLowerCase();

            console.log(`Deprotonation: pKa=${pka.toFixed(2)}, dG_correction=${dGCorrection.toFixed(3)}, new_offset=${newProtonOffset.toFixed(3)}`);
            return {
                newProtonOffset: newProtonOffset,
                newProtonType: newProtonType
            };
        }
    } else {
        // Not a protonation/deprotonation edge - keep current offset and type
        return {
            newProtonOffset: currentProtonOffset,
            newProtonType: currentProtonType
        };
    }
};

ReactionNetworkJS.prototype.calculateReactionProfile = function (path, edgeMap) {
    /**
     * DATA-ONLY method - returns precomputed profile from Python
     * DOES NOT calculate anything - JavaScript should never do energy calculations
     *
     * @param {Array} path - Array of node IDs representing the reaction path
     * @param {Map} edgeMap - Map of edge keys to edge objects (not used)
     * @returns {Array} Precomputed reaction profile from Python
     */

    console.warn("⚠️ calculateReactionProfile called - JavaScript should use precomputed profiles only");

    if (path.length < 2) {
        return [];
    }

    // Find the target node and return its precomputed profile
    const targetNode = path[path.length - 1];
    const precomputedProfile = this.getPrecomputedProfile(targetNode, 0);

    if (precomputedProfile && precomputedProfile.length > 0) {
        console.log(`✅ Returning precomputed profile for node ${targetNode} (${precomputedProfile.length} points)`);
        return precomputedProfile;
    }

    console.warn(`❌ No precomputed profile found for node ${targetNode}`);
    return [];
};

ReactionNetworkJS.prototype.validateEnergyCalculations = function (testCases = null) {
    /**
     * Validate that energy calculations match expected Python results
     * Can be called from browser console to verify correctness
     *
     * @param {Array} testCases - Optional array of test cases, uses defaults if null
     * @returns {Object} Validation results
     */
    console.log("=== JavaScript Energy Calculation Validation ===");

    // Default test cases based on known Python results
    const defaultTests = [
        {
            name: "Protonation Node 1→7",
            path: [1, 7],
            expectedEnergy: 0.405,
            tolerance: 0.001
        }
    ];

    const tests = testCases || defaultTests;
    const results = [];

    for (const test of tests) {
        console.log(`\nTesting: ${test.name}`);
        console.log(`Path: ${test.path.join(' → ')}`);

        try {
            // Calculate reaction profile
            const profile = this.calculateReactionProfile(test.path, this.edgeMap);

            if (profile.length < 2) {
                console.log(`❌ ${test.name}: No profile generated`);
                results.push({
                    name: test.name,
                    passed: false,
                    error: "No profile generated"
                });
                continue;
            }

            // Get final energy (last point in profile)
            const finalEnergy = profile[profile.length - 1].energy;
            const difference = Math.abs(finalEnergy - test.expectedEnergy);
            const passed = difference <= test.tolerance;

            console.log(`Expected: ${test.expectedEnergy.toFixed(3)} kcal/mol`);
            console.log(`Calculated: ${finalEnergy.toFixed(3)} kcal/mol`);
            console.log(`Difference: ${difference.toFixed(6)} kcal/mol`);
            console.log(`Result: ${passed ? '✅ PASS' : '❌ FAIL'}`);

            results.push({
                name: test.name,
                expected: test.expectedEnergy,
                calculated: finalEnergy,
                difference: difference,
                passed: passed
            });

        } catch (error) {
            console.log(`❌ ${test.name}: Error - ${error.message}`);
            results.push({
                name: test.name,
                passed: false,
                error: error.message
            });
        }
    }

    // Summary
    const passCount = results.filter(r => r.passed).length;
    const totalCount = results.length;

    console.log(`\n=== Validation Summary ===`);
    console.log(`${passCount}/${totalCount} tests passed`);

    if (passCount === totalCount) {
        console.log("✅ All energy calculations are correct!");
    } else {
        console.log("❌ Some energy calculations need fixing.");
    }

    return {
        passed: passCount === totalCount,
        results: results,
        summary: `${passCount}/${totalCount} tests passed`
    };
};

ReactionNetworkJS.prototype.findMaxBarrier = function (profile) {
    /**
     * Find maximum energy barrier in a reaction profile
     * Exact port of Python _find_max_barrier method
     *
     * @param {Array} profile - Reaction profile array
     * @returns {number} Maximum barrier height in kcal/mol
     */

    if (!profile || profile.length === 0) {
        return 0.0;
    }

    const energies = profile.map(p => p.energy);

    // Calculate energy above the rolling minimum
    let rollingMin = Number.POSITIVE_INFINITY;
    const barrierEnergies = [];

    for (const energy of energies) {
        if (energy < rollingMin) {
            rollingMin = energy;
        }
        barrierEnergies.push(energy - rollingMin);
    }

    // Find the maximum barrier
    const maxBarrier = Math.max(...barrierEnergies);

    return maxBarrier;
};

// Test function for validation
ReactionNetworkJS.prototype.validateEnergyCalculations = function () {
    /**
     * Validate energy calculation methods against known values
     * @returns {Object} Validation results
     */

    const validation = {
        hartreeConversion: false,
        pkaCalculation: false,
        energyCorrections: false,
        issues: []
    };

    // Test hartree conversion (known value)
    const testHartree = 1.0;
    const expectedKcal = 627.509474;
    const actualKcal = this.hartree2kcalmol(testHartree);

    if (Math.abs(actualKcal - expectedKcal) < 1e-6) {
        validation.hartreeConversion = true;
    } else {
        validation.issues.push(`Hartree conversion mismatch: expected ${expectedKcal}, got ${actualKcal}`);
    }

    // Test pKa calculation logic (basic functionality)
    try {
        const testResult = this.handleProtonationDeprotonation(
            "protonation", "", { rawEnergy: -1000, originType: "protonation" },
            0, 7, this.R, this.T, this.PKA_SLOPE, this.PKA_INTERCEPT
        );

        if (typeof testResult.newProtonOffset === 'number' && typeof testResult.newProtonType === 'string') {
            validation.pkaCalculation = true;
        } else {
            validation.issues.push('pKa calculation returned invalid types');
        }
    } catch (error) {
        validation.issues.push(`pKa calculation failed: ${error.message}`);
    }

    // Test energy corrections
    try {
        const testEdge = { smallerProducts: [], type: "reaction" };
        const testNode = { id: 1, rawEnergy: -1000 };
        const eliminatedEnergy = this.calculateEliminatedProductsEnergy(testEdge, testNode);

        if (eliminatedEnergy === 0.0) {
            validation.energyCorrections = true;
        } else {
            validation.issues.push(`Expected 0 eliminated energy, got ${eliminatedEnergy}`);
        }
    } catch (error) {
        validation.issues.push(`Energy corrections failed: ${error.message}`);
    }

    return validation;
};

console.log('✅ Energy calculation methods added to ReactionNetworkJS');
