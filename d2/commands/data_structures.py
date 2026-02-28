"""Data structures for the network visualization refactoring."""

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional


@dataclass
class VisualizationData:
    """Container for all visualization data that gets passed to JavaScript."""

    svg_map: Dict[int, str]
    mol_energy_map: Dict[int, float]
    charge_map: Dict[int, int]
    weight_map: Dict[int, int]
    original_edge_data: List[Dict[str, Any]]
    substrate_id: int
    energy_type: str
    temperature: float
    fit_params: Dict[str, Any]
    stereoisomer_groups: Optional[Dict[int, List[int]]] = None
    node_labels: Optional[Dict[int, str]] = None

    def to_json(self) -> str:
        """Convert to JSON for JavaScript consumption."""
        return json.dumps(
            {
                "svgMap": self.svg_map,
                "molEnergyMap": self.mol_energy_map,
                "chargeMap": self.charge_map,
                "weightMap": self.weight_map,
                "originalEdgeData": self.original_edge_data,
                "substrateId": self.substrate_id,
                "energyType": self.energy_type,
                "temperature": self.temperature,
                "fitParams": self.fit_params,
                "stereoisomerGroups": self.stereoisomer_groups or {},
                "nodeLabels": self.node_labels or {},
            }
        )


class NetworkDataProcessor:
    """Handles the processing of molecular and reaction data."""

    def __init__(self, config, substrate_id: int, prop_name: str):
        self.config = config
        self.substrate_id = substrate_id
        self.prop_name = prop_name
        self.energy_type = "G" if "gibbs" in prop_name else "E"

    def process_molecules(self, visualizer) -> Dict[str, Any]:
        """Process molecular data and return structured information."""
        return {
            "svg_map": visualizer.mol_svg_map,
            "energy_map": visualizer.mol_energy_map,
            "charge_map": visualizer.mol_charge_map,
            "weight_map": visualizer.mol_weight_map,
        }

    def process_reactions(self, edges) -> List[Dict[str, Any]]:
        """Process reaction data and return edge information."""
        return [
            {
                "id": i,
                "begin": edge["begin"],
                "end": edge["end"],
                "type": edge["type"],
                "deltaE": edge["deltaE"],
                "count": edge.get("count", 1),
                "smaller_products": edge.get("smaller_products", []),
            }
            for i, edge in enumerate(edges)
        ]


class HTMLTemplateRenderer:
    """Handles HTML template generation and rendering."""

    def __init__(self, template_dir: Optional[Path] = None):
        self.template_dir = template_dir or Path(__file__).parent / "templates"

    def render_main_template(self, visualization_data: VisualizationData) -> str:
        """Render the main HTML template with data."""
        template_path = self.template_dir / "network_visualization.html"

        if not template_path.exists():
            return self._generate_inline_template(visualization_data)

        with open(template_path, "r") as f:
            template = f.read()

        # Replace the placeholder with the actual data
        html_content = template.replace(
            "{{ VISUALIZATION_DATA }}", visualization_data.to_json()
        )

        # Update relative paths to account for the output location
        commands_dir = str(Path(__file__).parent)
        html_content = html_content.replace(
            'href="static/', f'href="{commands_dir}/static/'
        )
        html_content = html_content.replace(
            'src="static/', f'src="{commands_dir}/static/'
        )

        return html_content

    def _generate_inline_template(self, data: VisualizationData) -> str:
        """Fallback inline template generation."""
        commands_dir = str(Path(__file__).parent)
        return f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Reaction Network Visualization</title>
    <link rel="stylesheet" href="{commands_dir}/static/css/network-visualization.css">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
</head>
<body>
    <div id="network-container"></div>
    <div id="controls-panel"></div>

    <script>
        window.VISUALIZATION_DATA = {data.to_json()};
    </script>
    <script src="{commands_dir}/static/js/network-config.js"></script>
    <script src="{commands_dir}/static/js/network-data.js"></script>
    <script src="{commands_dir}/static/js/network-filters.js"></script>
    <script src="{commands_dir}/static/js/network-ui.js"></script>
    <script src="{commands_dir}/static/js/network-utils.js"></script>
    <script src="{commands_dir}/static/js/network-main.js"></script>
</body>
</html>
        """
