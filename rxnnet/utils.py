import hashlib
import os
from pathlib import Path
from typing import Tuple
import shutil
import random
import subprocess
from typing import Generator
import signal
import warnings
import string

STANDARD_PROPERTIES = {"xtb": {"total energy": "electronic_energy"}, "orca": {}}

alphabet = string.ascii_lowercase + string.digits


def get_random_str() -> str:
    """Generate a random string of 6 characters."""
    return hashlib.md5(os.urandom(16)).hexdigest()[:6]


def fragment_charge(
    total_charge: int, n_fragments: int, min_max_frag_charge: Tuple[int, int] = (-2, 2)
):
    """Generate possible charge combinations for molecular fragments."""
    from itertools import product

    min_max_frag_charge = list(min_max_frag_charge)
    min_max_frag_charge[1] += 1
    min_max_frag_charge = range(*tuple(min_max_frag_charge))
    return [
        combo
        for combo in product(min_max_frag_charge, repeat=n_fragments)
        if sum(combo) == total_charge
    ]


def stream(
    cmd: str, cwd: None | Path = None, shell: bool = True
) -> Generator[str, None, None]:
    """Execute a command and stream stdout and stderr concurrently."""
    with subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        shell=shell,
        cwd=cwd,
        preexec_fn=os.setsid,
        bufsize=1,
    ) as process:
        try:
            for line in iter(process.stdout.readline, ""):
                yield line
            for line in iter(process.stderr.readline, ""):
                yield line
        except KeyboardInterrupt:
            print("\nCtrl+C pressed. Terminating the process...")
            os.killpg(os.getpgid(process.pid), signal.SIGTERM)
            process.wait()
            print("Process terminated.")
        finally:
            process.stdout.close()
            process.stderr.close()
            process.wait()


def check_executable(executable: str):
    """Check if executable is in PATH."""
    results = stream(f"which {executable}")
    result = next(results)
    if result.startswith("which: no"):
        warnings.warn(f"Executable {executable} not found in PATH")


class WorkingDir:
    def __init__(self, root: str = ".", name: str = None) -> None:
        self.root = Path(root)
        self.name = name if name else self._random_str()
        self.dir = self.root / self.name
        self.create()

    def __str__(self) -> str:
        return str(self.dir.resolve())

    def __repr__(self) -> str:
        return self.__str__()

    def __truediv__(self, name: str) -> str:
        return self.dir / name

    def _random_str(self) -> str:
        name = "_" + "".join(random.choices(alphabet, k=6))
        while (self.root / name).exists():
            name = "".join(random.choices(alphabet, k=6))
        return name

    def create(self) -> None:
        self.dir.mkdir(parents=True, exist_ok=True)

    def cleanup(self) -> None:
        try:
            # print("removing ", self.dir.absolute())
            shutil.rmtree(self.dir.absolute())
        except FileNotFoundError:
            pass


class WorkingFile:
    def __init__(self, root: str = ".", filename: str = None, mode="w") -> None:
        self.root = Path(root)
        self.filename = filename if filename else self._random_str()
        self.mode = mode
        self.path = self.root / self.filename

    def _random_str(self) -> str:
        name = "".join(random.choices(alphabet, k=6)) + ".ttxt"
        while (self.root / name).exists():
            name = "".join(random.choices(alphabet, k=6)) + ".ttxt"
        return name

    def __str__(self) -> str:
        return str(self.path.resolve())

    def __repr__(self) -> str:
        return self.__str__()

    def create(self) -> None:
        with open(str(self), self.mode) as _:
            pass

    def cleanup(self) -> None:
        try:
            shutil.rmtree(self.path)
        except FileNotFoundError:
            pass

    @property
    def stem(self):
        return str(self.path.stem)
