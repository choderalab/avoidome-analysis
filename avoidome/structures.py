import numpy as np
from pydantic import BaseModel, Field

from avoidome.schema import ProteinEntity


class StructureEntry(BaseModel):
    """
    A class to represent a structure entry
    """

    components: list[ProteinEntity | Ligand] = Field(
        ..., title="The components of the structure"
    )


class ExperimentalStructure(StructureEntry):
    """
    A class to represent an experimental structure entry
    """

    pdb_id: str = Field(..., title="The PDB ID of the structure")
    method: str = Field(..., title="The method used to determine the structure")
    resolution: float = Field(..., title="The resolution (Å) of the structure")

    @property
    def sequence_coverage(self):
        return np.mean([c.length for c in self.components])


class PredictedStructure(StructureEntry):
    """
    A class to represent a predicted structure entry
    """

    af_id: str = Field(..., title="The AlphaFold DB ID of the structure")
    confidence: float = Field(..., title="The confidence of the structure")
    model_url: str = Field(..., title="The URL of the PDB file of the structure")
