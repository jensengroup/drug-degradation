"""Network visualization module."""

import json
import webbrowser
from pathlib import Path
from typing import Any, Dict, Optional

import typer
from rdkit.Chem import rdDepictor

from d2.commands.rn import ReactionNetwork  # Use unified ReactionNetwork
from d2.config import NetworkConfig

# Ensure consistent 2D depiction
rdDepictor.SetPreferCoordGen(True)

app = typer.Typer()


class HTMLTemplateRenderer:
    """Handles HTML template generation and rendering."""

    def __init__(
        self, template_dir: Optional[Path] = None, static_dir: Optional[Path] = None
    ):
        self.template_dir = template_dir or Path(__file__).parent / "templates"
        self.static_dir = static_dir or Path(__file__).parent / "static"

    def render_template(self, visualization_data: Dict[str, Any]) -> str:
        """Render the HTML template with visualization data as a self-contained HTML file."""
        template_path = self.template_dir / "network_visualization.html"

        if not template_path.exists():
            return self._generate_inline_template(visualization_data)

        with open(template_path, "r") as f:
            template = f.read()

        # Replace the placeholder with the actual data
        html_content = template.replace(
            "{{ VISUALIZATION_DATA }}", json.dumps(visualization_data)
        )

        # Inline all CSS files
        html_content = self._inline_css_files(html_content)

        # Inline all JavaScript files (except external CDN links)
        html_content = self._inline_js_files(html_content)

        return html_content

    def _inline_css_files(self, html_content: str) -> str:
        """Replace CSS link tags with inline styles."""
        import re

        # Find all CSS link tags that reference local files - handle version parameters
        css_pattern = (
            r'<link[^>]*href=["\']static/([^"\']+\.css)(?:\?[^"\']*)?["\'][^>]*>'
        )

        def replace_css_link(match):
            css_path = match.group(1)
            css_file_path = self.static_dir / css_path

            if css_file_path.exists():
                try:
                    with open(css_file_path, "r", encoding="utf-8") as f:
                        css_content = f.read()
                    return f"<style>\n{css_content}\n</style>"
                except Exception as e:
                    print(f"Could not read CSS file {css_file_path}: {e}")
                    return ""
            else:
                print(f"CSS file not found: {css_file_path}")
                return ""

        return re.sub(css_pattern, replace_css_link, html_content)

    def _inline_js_files(self, html_content: str) -> str:
        """Replace JavaScript script tags with inline scripts (except external CDN links)."""
        import re

        # Find all script tags that reference local files (not CDN) - handle version parameters
        js_pattern = r'<script[^>]*src=["\']static/([^"\']+\.js)(?:\?[^"\']*)?["\'][^>]*></script>'

        def replace_js_script(match):
            js_path = match.group(1)
            js_file_path = self.static_dir / js_path

            # Skip debug files
            if any(debug_name in js_path for debug_name in ["debug-", "test-"]):
                return ""

            if js_file_path.exists():
                try:
                    with open(js_file_path, "r", encoding="utf-8") as f:
                        js_content = f.read()
                    return f"<script>\n{js_content}\n</script>"
                except Exception as e:
                    print(f"Could not read JS file {js_file_path}: {e}")
                    return ""
            else:
                print(f"JS file not found: {js_file_path}")
                return ""

        return re.sub(js_pattern, replace_js_script, html_content)

    def _generate_inline_template(self, data: Dict[str, Any]) -> str:
        """Fallback inline template generation."""
        css_content = ""
        js_files_content = []

        try:
            css_path = self.static_dir / "css" / "network-visualization.css"
            if css_path.exists():
                with open(css_path, "r", encoding="utf-8") as f:
                    css_content = f.read()

            # Read JavaScript files in order
            js_files = [
                "network-config.js",
                "network-data.js",
                "network-path-analysis.js",
                "network-filters.js",
                "network-stereoisomers.js",
                "network-ui.js",
                "network-utils.js",
                "network-data-viewer.js",
                "network-main.js",
            ]

            for js_file in js_files:
                js_path = self.static_dir / "js" / js_file
                if js_path.exists():
                    with open(js_path, "r", encoding="utf-8") as f:
                        js_files_content.append(f.read())

        except Exception as e:
            print(f"Error reading static files in fallback: {e}")

        # Combine all JavaScript
        combined_js = "\n\n".join(js_files_content)

        return f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>Reaction Network Visualization</title>
    <script src="https://unpkg.com/vis-network/standalone/umd/vis-network.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
{css_content}
    </style>
</head>
<body>
    <div id="mynetwork"></div>
    <div id="moleculeViewer" class="molecule-viewer">
        <!-- Controls and UI will be dynamically generated by JavaScript -->
    </div>

    <script>
        window.VISUALIZATION_DATA = {json.dumps(data)};
    </script>
    <script>
{combined_js}
    </script>
</body>
</html>"""


@app.command("vis")
def vis_command(
    substrate_id: int = 1,
    lot: str = "gfn2-xtb",
    prop_name: str = "l2l1_gibbs-energy",
    pH: float = 7.0,
):
    """Visualize the reaction network."""
    config = NetworkConfig()

    if not config.is_initialized():
        print("Network not initialized. Run setup command first.")
        return

    print(f"Creating reaction network visualization (pH = {pH})...")

    # Create network object and process data
    network = ReactionNetwork(config, substrate_id, prop_name, pH)

    print("Generating self-contained HTML visualization...")

    visualization_data = network.to_visualization_data()

    template_renderer = HTMLTemplateRenderer()
    html_content = template_renderer.render_template(visualization_data)

    output_path = Path(".reaction_network_vis.html")
    output_path.write_text(html_content, encoding="utf-8")

    print(f"Saved visualization to: {output_path.resolve()}")
    webbrowser.open(f"file://{output_path.resolve()}")
