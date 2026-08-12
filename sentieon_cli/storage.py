"""Storage staging for job inputs and outputs."""

import pathlib
from abc import ABC, abstractmethod


class StorageProvider(ABC):
    """Stage a job's declared files in before it runs and out after.

    On a shared filesystem the files are already local, so staging
    degenerates to verifying presence (see :class:`LocalStorageProvider`).
    A cloud provider would instead pull inputs from and push outputs to
    object storage, deriving each object's URI from a configured prefix and
    the local path. The methods are async so they compose with the async
    executors.
    """

    @abstractmethod
    async def stage_in(self, local: pathlib.Path) -> None:
        """Make ``local`` available before a job runs.

        Pull it from storage, or verify it is present. Raise if it cannot be
        made available.
        """

    @abstractmethod
    async def stage_out(self, local: pathlib.Path) -> None:
        """Persist ``local`` after a job runs.

        Push it to storage, or verify it was produced. Raise if it is missing.
        """


class LocalStorageProvider(StorageProvider):
    """Shared-filesystem provider: no data movement, just verify presence.

    A file counts as present when it exists, even if it is empty (zero bytes).
    """

    async def stage_in(self, local: pathlib.Path) -> None:
        if not local.exists():
            raise FileNotFoundError(f"missing input: {local}")

    async def stage_out(self, local: pathlib.Path) -> None:
        if not local.exists():
            raise FileNotFoundError(f"missing output: {local}")
