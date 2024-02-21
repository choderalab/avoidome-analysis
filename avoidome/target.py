import numpy as np
from pydantic import Field

from avoidome.structures import StructureEntry, ExperimentalStructure, PredictedStructure


class TargetStructureData(Target):
    """
    A class to represent a protein target
    """

    name: str = Field(..., title="The name of the protein")
    structures: list[StructureEntry] = Field(..., title="A list of structure entries")

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
    def from_uniprot_entry
