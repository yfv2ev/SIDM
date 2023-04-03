"""Module to define classes and methods that are helpful for scaleout"""

from dask.distributed import Client, Worker, WorkerPlugin
from typing import List
import os

class DependencyInstaller(WorkerPlugin):
    """Class to install dependencies on dask workers"""
    def __init__(self, dependencies: List[str]):
        self._dependencies = " ".join(f"'{dep}'" for dep in dependencies)

    def setup(self, worker: Worker):
        os.system(f"pip install {self._dependencies}")


def make_dask_client(address):
    """Create dask client that includes dependency installer"""
    dependency_installer = DependencyInstaller([
        "git+https://github.com/btcardwell/SIDM.git@scaleout",
    ])
    client = Client(address)
    client.register_worker_plugin(dependency_installer)
    return client
