import json
import os
import shutil

try:
    import tomllib
except ModuleNotFoundError:
    import pip._vendor.tomli as tomllib

from pathlib import Path

import rxnnet


def load_toml(path: str | Path) -> dict:
    """Load settings from a TOML file."""
    if path and Path(path).is_file():
        with open(path, "rb") as f:
            return tomllib.load(f)
    raise FileNotFoundError(f"Config file not found: {path}")


class Config:
    """Network directory structure and settings."""

    def __init__(self, network_dir: str = "."):
        self.network_dir = Path(network_dir).resolve()

        # Directory structure
        self.network_file = self.network_dir / "network-info.json"
        self.node_data = self.network_dir / "nodes"
        self.product_data = self.network_dir / "products"
        self.new_nodes = self.network_dir / "new-nodes"
        self.reaction_data = self.network_dir / "reactions"
        self.new_reactions = self.network_dir / "new-reactions"
        self.qm_data = self.network_dir / "qm-data"
        self.tmp_dir = self.network_dir / ".tmp"
        self.label_db = self.network_dir / ".label-db"

        config_path = self.network_dir / "config.toml"
        if not config_path.exists():
            template = Path(rxnnet.__path__[0]) / "config_template.toml"
            shutil.copy(template, self.network_dir / "config.toml")

        self.settings = load_toml(config_path)

        # Network state (from network-info.json)
        self.state = self._load_network_state()

        # add executables to path
        for _, path in self.settings["executables"].items():
            os.environ["PATH"] = str(Path(path).parent) + ":" + os.environ["PATH"]

    def _load_network_state(self) -> dict:
        """Load network state from network-info.json."""
        if self.network_file.is_file():
            with open(self.network_file, "r") as f:
                return json.load(f)
        return {}

    # Network state (from network-info.json)
    @property
    def node_count(self) -> int:
        return self.state.get("node_count", 0)

    def is_initialized(self) -> bool:
        return self.network_file.is_file()

    def setup(self) -> None:
        """Initialize the network: create directories and network-info.json."""
        self._create_directories()

        # Create initial network state
        initial_state = {
            "node_count": 0,
            "tmp_count": 1000000,
        }

        with open(self.network_file, "w") as f:
            json.dump(initial_state, f, indent=2)

        # Reload state
        self.state = self._load_network_state()

    def _create_directories(self) -> None:
        """Create necessary directories for the network."""
        directories = [
            self.node_data,
            self.reaction_data,
            self.product_data,
            self.new_nodes,
            self.qm_data,
            self.tmp_dir,
            self.new_reactions,
        ]

        for directory in directories:
            os.makedirs(directory.resolve(), exist_ok=True)


class IdxCounter:
    """Manages index counters for nodes."""

    def __init__(self, parent_file: str, name: str):
        self.parent_file = parent_file
        self.name = name

    @property
    def idx(self) -> int:
        """Get current index value."""
        with open(self.parent_file, "r") as f:
            global_settings = json.load(f)
        return global_settings.get(self.name, 0)

    def increment(self) -> int:
        """Increment counter and return new value."""
        with open(self.parent_file, "r") as f:
            global_settings = json.load(f)
        global_settings[self.name] += 1
        with open(self.parent_file, "w") as f:
            json.dump(global_settings, f, indent=2)
        return global_settings[self.name]


def main():
    """CLI entry point: Initialize a reaction network."""
    import argparse

    parser = argparse.ArgumentParser(description="Initialize a reaction network")
    parser.add_argument(
        "network_dir",
        nargs="?",
        default=".",
        help="Network directory path (default: current directory)",
    )
    args = parser.parse_args()

    config = Config(args.network_dir)
    if not config.is_initialized():
        config.setup()
        print(f"Network initialized in {config.network_dir}")
    else:
        print(f"Network already initialized in {config.network_dir}")


if __name__ == "__main__":
    main()
