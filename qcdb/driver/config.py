from enum import Enum
from typing import Any, Dict

from pydantic import BaseModel, Field, root_validator


class FrameEnum(str, Enum):

    native = "native"
    unified = "unified"


class ProgramEnum(str, Enum):

    cfour = "cfour"
    gamess = "gamess"
    nwchem = "nwchem"
    psi4 = "psi4"
    qcdb = "qcdb"


class ModeEnum(str, Enum):

    sandwich = "sandwich"
    unified = "unified"


class ModeConfig(BaseModel):
    """Description of the configuration used to formulate and interpret a task."""

    # Specifications
    #    result_frame: FrameEnum = Field(None, description="Orientation/frame/alignment for returned results")
    mode: ModeEnum = Field(None, description="Broad mode of operation specifying degree of interoperability. `sandwich` aims for unified interface. `unified` aims for unified results. Usually, set only this as it sets further fields, but other fields can be tweaked to override broad mode.")
    implicit_program: ProgramEnum = Field(
        ProgramEnum.qcdb, description="Program to use without prefixing method or keywords"
    )
    module_fallback: bool = Field(None, description="")
    translate_method_algorithm: bool = Field(None, description="Suggest the QCDB value for QC method algorithm, e.g., conventional vs. density-fitted")
    translate_orbital_space: bool = Field(None, description="Suggest the QCDB value for frozen core")

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

    @root_validator(pre=True)
    def setup_modes_aggregate(cls, values):
        mode = values.get("mode")
        if mode == "unified":
            values["module_fallback"] = values.get("module_fallback", True)
            values["translate_method_algorithm"] = values.get("translate_method_algorithm", True)
            values["translate_orbital_space"] = values.get("translate_orbital_space", True)
        elif mode == "sandwich":
            values["module_fallback"] = values.get("module_fallback", False)
            values["translate_method_algorithm"] = values.get("translate_method_algorithm", False)
            values["translate_orbital_space"] = values.get("translate_orbital_space", False)
        return values


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
