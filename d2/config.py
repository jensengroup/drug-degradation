"""Configuration management for the reaction network application."""

import json
import os
from pathlib import Path
from typing import Optional

import submitit

try:
    import typer
    from rich.prompt import Prompt

    HAS_RICH_TYPER = True
except ImportError:
    HAS_RICH_TYPER = False

    # Fallback for basic input
    class Prompt:
        @staticmethod
        def ask(question, default=None, choices=None, show_default=False):
            if choices:
                question += f" [{'/'.join(choices)}]"
            if default is not None:
                question += f" (default: {default})"
            response = input(f"{question}: ").strip()
            return response if response else default

    class typer:
        @staticmethod
        def confirm(question, default=True):
            response = (
                input(f"{question} [{'Y/n' if default else 'y/N'}]: ").strip().lower()
            )
            if not response:
                return default
            return response.startswith("y")

        class BadParameter(Exception):
            pass


class NetworkConfig:
    """Manages network configuration and global settings."""

    def __init__(self, network_dir: str = "."):
        self.network_dir = Path(network_dir).resolve()
        self.network_file = self.network_dir / "network-info.json"
        self.node_data = self.network_dir / "nodes"
        self.product_data = self.network_dir / "products"
        self.new_nodes = self.network_dir / "new-nodes"
        self.reaction_data = self.network_dir / "reactions"
        self.new_reactions = self.network_dir / "new-reactions"
        self.qm_data = self.network_dir / "qm-data"
        self.tmp_dir = self.network_dir / ".tmp"
        self.label_db = self.network_dir / ".label-db"

        # Energy offsets for different solvents
        self.offsets = {
            "water": {
                "deprotonation": -0.4059596304589945,
                "protonation": 0.4059596304589945,
            }
        }

        self.fit_params = {"a": 0.43496246, "b": -114.66798804}

        # Solvent choices
        self.solvent_choices = [
            "None",
            "Acetone",
            "Acetonitrile",
            "CH2Cl2",
            "CHCl3",
            "DMF",
            "DMSO",
            "Hexane",
            "Methanol",
            "Octanol",
            "THF",
            "Water",
        ]

        # Load settings
        self._load_settings()

    def _load_settings(self) -> None:
        """Load global settings from network file."""
        if self.network_file.is_file():
            with open(self.network_file, "r") as f:
                self.global_settings = json.load(f)
            self.temperature = self.global_settings["temperature"]
            self.pressure = self.global_settings["pressure"]
            self.solvent = self.global_settings.get("solvent", None)
            self.node_count = self.global_settings.get("node_count", 0)
            self.tmp_count = self.global_settings.get("tmp_count", 1000000)
            self.lot = self.global_settings.get("lot", "gfn2-xtb")

            # Load executor settings
            self.executor_settings = self.global_settings.get(
                "executor",
                {
                    "slurm_partition": "kemi1",
                    "timeout_min": 1200,
                    "slurm_mem_per_cpu": "2000",
                    "cpus_per_task": 1,
                    "slurm_array_parallelism": 100,
                    "slurm_job_name": "d2",
                },
            )
        else:
            print("No `network-info.json` file found, run `control.py setup`")
            self.node_count = 0
            self.temperature = None
            self.pressure = None
            self.solvent = None
            self.lot = "gfn2-xtb"
            self.executor_settings = {
                "slurm_partition": "kemi1",
                "timeout_min": 1200,
                "slurm_mem_per_cpu": "2000",
                "cpus_per_task": 1,
                "slurm_array_parallelism": 100,
                "slurm_job_name": "d2",
            }

    @property
    def executor(self):
        return self.get_executor()

    def is_initialized(self) -> bool:
        """Check if network is initialized."""
        return self.network_file.is_file()

    def setup_network(
        self,
        temperature: Optional[float] = None,
        pressure: Optional[float] = None,
        solvent: Optional[str] = None,
        lot: Optional[str] = None,
    ) -> None:
        """Initialize network with given parameters."""
        if self.is_initialized():
            print(f"Network is already initialized in {self.network_dir}.")
            return

        print(f"No Network initialized in {self.network_dir}.")
        if not typer.confirm("Do you want to initialize a new Network?", default=True):
            return

        # Get missing parameters
        if lot is None:
            lot = Prompt.ask(
                "Enter default level of theory (LOT)",
                default="r2scan-3c/gfn2-xtb",
                show_default=True,
            ).lower()

        if temperature is None:
            temperature = Prompt.ask("Enter temperature (K)", default=298.15)
        if pressure is None:
            pressure = Prompt.ask("Enter pressure (atm)", default=1.0)
        if solvent is None:
            solvent = Prompt.ask(
                "Enter solvent",
                choices=self.solvent_choices,
                default="None",
            )

        # Validate solvent
        if solvent is not None and solvent not in self.solvent_choices:
            raise typer.BadParameter(
                f"Invalid solvent choice. Must be one of {self.solvent_choices}"
            )

        print(
            f"Initializing network with Temperature={temperature}, "
            f"Pressure={pressure}, Solvent={solvent}"
        )

        # Create global settings
        global_settings = {
            "temperature": temperature,
            "pressure": pressure,
            "node_count": 0,
            "tmp_count": 1000000,
            "lot": lot,
            "executor": {
                "slurm_partition": "kemi1",
                "timeout_min": 1200,
                "slurm_mem_per_cpu": "2000",
                "cpus_per_task": 1,
                "slurm_array_parallelism": 100,
                "slurm_job_name": "d2",
            },
        }

        if solvent != "None":
            global_settings["solvent"] = solvent.lower()

        # Save settings and create directories
        with open(self.network_file, "w") as f:
            json.dump(global_settings, f, indent=2)

        self._create_directories()
        self._load_settings()  # Reload settings

        print(f"Network initialized at {self.network_file}")

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

    def get_executor(self, **kwargs) -> submitit.Executor:
        """Get submitit AutoExecutor for job submission."""
        # Always use AutoExecutor - automatically chooses between local and cluster
        executor = submitit.AutoExecutor(folder=self.tmp_dir)

        # Configure with settings, allowing kwargs to override defaults
        settings = {**self.executor_settings, **kwargs}

        executor.update_parameters(
            timeout_min=settings.get("timeout_min", 1200),
            slurm_mem_per_cpu=settings.get("slurm_mem_per_cpu", "2000"),
            slurm_array_parallelism=settings.get("slurm_array_parallelism", 100),
            cpus_per_task=settings.get("cpus_per_task", 1),
            slurm_job_name=settings.get("slurm_job_name", "d2"),
            slurm_partition=settings.get("slurm_partition", "kemi1"),
        )

        return executor

    def configure_executor(
        self,
        **kwargs,
    ) -> submitit.Executor:
        """Configure and return a submitit AutoExecutor with custom
        parameters."""
        # Always use AutoExecutor - automatically chooses between local and cluster
        executor = self.get_executor()
        executor.update_parameters(
            **kwargs,
        )

        return executor


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
