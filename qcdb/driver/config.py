from enum import Enum
from typing import Any, Dict

from pydantic import BaseModel, Field


class FrameEnum(str, Enum):

    native = "native"
    unified = "unified"


class ProgramEnum(str, Enum):

    cfour = "cfour"
    gamess = "gamess"
    nwchem = "nwchem"
    psi4 = "psi4"
    qcdb = "qcdb"


class ModeConfig(BaseModel):
    """Description of the configuration used to formulate and interpret a task."""

    # Specifications
    #    result_frame: FrameEnum = Field(None, description="Orientation/frame/alignment for returned results")
    implicit_program: ProgramEnum = Field(
        ProgramEnum.qcdb, description="Program to use without prefixing method or keywords"
    )
    module_fallback: bool = Field(None, description="")

    #    ncores: int = pydantic.Field(None, description="Number cores per task on each node")
    #    nnodes: int = pydantic.Field(None, description="Number of nodes per task")
    #    memory: float = pydantic.Field(
    #        None, description="Amount of memory in GiB (2^30 bytes; not GB = 10^9 bytes) per node."
    #    )
    #    scratch_directory: Optional[str]  # What location to use as scratch
    #    retries: int  # Number of retries on random failures
    #    mpiexec_command: Optional[str]  # Command used to launch MPI tasks, see NodeDescriptor
    #    use_mpiexec: bool = False  # Whether it is necessary to use MPI to run an executable
    #    cores_per_rank: int = pydantic.Field(1, description="Number of cores per MPI rank")
    #    scratch_messy: bool = pydantic.Field(
    #        False, description="Leave scratch directory and contents on disk after completion."
    #    )

    class Config:
        extra = "forbid"


# def get_mode_config(key: str = None):
#    if key is None:
#        return mode_toggles
#    else:
#        return mode_toggles[key]


def get_mode_config(*, mode_options: Dict[str, Any] = None) -> ModeConfig:
    """
    Returns the mode configuration set for qcdb.
    """

    if mode_options is None:
        mode_options = {}

    config = {}

    # Override any settings
    if mode_options:
        config.update(mode_options)

    return ModeConfig(**config)


# def set_mode_config(key: str, val):
#    mode_toggles[key] = val
#
#
# mode_toggles = {
#    "unfixed_output_orientation": None,  # "dsl", "unified"
# }
