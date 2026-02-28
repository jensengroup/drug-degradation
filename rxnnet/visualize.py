import json
from pathlib import Path
from typing import Optional

from rxnnet.network import ReactionNetwork


class Visualizer:
    """Template-based visualization generator."""

    def __init__(self, output_dir: Optional[str | Path] = None):
        self.template_dir = Path(__file__).parent / "templates"
        self.output_dir = Path(output_dir) if output_dir else Path.cwd()

    def render(
        self, network: ReactionNetwork, output_name: str = "network.html"
    ) -> Path:
        """Render network visualization to HTML file."""
        self.output_dir.mkdir(parents=True, exist_ok=True)
        output_path = self.output_dir / output_name

        data = network.to_visualization_data()

        base = self._load_template("base.html")
        sidebar = self._load_template("components/sidebar.html")
        modals = self._load_template("components/modals.html")
        css = self._load_template("css/main.css")

        js_modules = [
            "config",
            "data",
            "network-graph",
            "chart",
            "ui",
            "main",
        ]
        js_code = "\n".join([self._load_template(f"js/{m}.js") for m in js_modules])

        html = base.replace("{{SIDEBAR}}", sidebar)
        html = html.replace("{{MODALS}}", modals)
        html = html.replace("{{STYLES}}", css)
        html = html.replace("{{SCRIPTS}}", js_code)
        html = html.replace("{{NETWORK_DATA}}", json.dumps(data, default=str))

        output_path.write_text(html)
        print(f"Visualization saved to {output_path}")
        return output_path

    def _load_template(self, name: str) -> str:
        path = self.template_dir / name
        if not path.exists():
            print(f"Warning: Template not found: {path}")
            return ""
        return path.read_text()


def visualize(
    network_dir: str | Path,
    output_path: Optional[str | Path] = None,
    substrate_id: int = 1,
    pH: float = 7.0,
    temperature: float = 313.15,
) -> Path:
    """Generate visualization HTML for a reaction network."""
    network = ReactionNetwork(
        network_dir=network_dir,
        substrate_id=substrate_id,
        pH=pH,
        temperature=temperature,
    )

    network.compute_pathways()

    output_dir = Path(output_path).parent if output_path else Path(network_dir)
    output_name = Path(output_path).name if output_path else "network.html"

    viz = Visualizer(output_dir=output_dir)
    return viz.render(network, output_name=output_name)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Visualize reaction network")
    parser.add_argument(
        "-d",
        "--network-dir",
        default=".",
        help="Network directory path (default: current directory)",
    )
    parser.add_argument("-o", "--output", help="Output HTML path")
    parser.add_argument("-s", "--substrate", default=1, type=int, help="Substrate ID")
    parser.add_argument("--pH", type=float, default=7.0, help="pH")
    parser.add_argument("--temp", type=float, default=313.15, help="Temperature (K)")

    args = parser.parse_args()

    visualize(
        network_dir=args.network_dir,
        output_path=args.output,
        substrate_id=args.substrate,
        pH=args.pH,
        temperature=args.temp,
    )
