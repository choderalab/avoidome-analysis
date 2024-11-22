import numpy as np
from pydantic import Field, BaseModel

from avoidome.structures import (
    StructureEntry,
    ExperimentalStructure,
    PredictedStructure,
)
from avoidome.uniprot import UniprotEntry


class TargetStructureData(BaseModel):
    """
    A class to represent a protein target
    """

    uniprot_id: str = Field(..., title="The UniProt ID of the protein")
    target_name: str = Field(..., title="The name of the protein target")
    structures: list[StructureEntry] = Field(..., title="A list of structure entries")
    sequence: str = Field(..., title="The sequence of the protein")

    @property
    def sequence_length(self) -> int:
        return len(self.sequence)

    @property
    def experimental_structures(self):
        return [s for s in self.structures if isinstance(s, ExperimentalStructure)]

    @property
    def predicted_structures(self):
        return [s for s in self.structures if isinstance(s, PredictedStructure)]

    @property
    def average_coverage(self):
        return np.mean(
            np.array([s.sequence_coverage for s in self.experimental_structures])
            / self.sequence_length
        )

    @property
    def average_confidence(self):
        return np.mean(np.array([s.confidence for s in self.predicted_structures]))

    @property
    def n_experimental_structures(self):
        return len([s for s in self.structures if isinstance(s, ExperimentalStructure)])

    @property
    def n_predicted_structures(self):
        return len([s for s in self.structures if isinstance(s, PredictedStructure)])

    @classmethod
    def from_uniprot_entry(cls, uniprot: UniprotEntry):
        return cls(
            uniprot_id=uniprot.uniprot_id,
            target_name=uniprot.name,
            sequence=uniprot.sequence,
            structures=uniprot.get_experimental_structures()
            + uniprot.get_alphafold_structures(),
        )

    @classmethod
    def from_dict(cls, data: dict):
        return cls(
            uniprot_id=data["uniprot_id"],
            target_name=data["target_name"],
            sequence=data["sequence"],
            structures=[
                ExperimentalStructure(**s) if "pdb_id" in s else PredictedStructure(**s)
                for s in data["structures"]
            ],
        )
