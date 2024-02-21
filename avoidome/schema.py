from pydantic import BaseModel, Field
import numpy as np
from asapdiscovery.data.schema_v2.target import Target
from asapdiscovery.data.schema_v2.ligand import Ligand


class ProteinEntity(BaseModel):
    """
    A class to represent a single protein entity using the UniProt ID as a system of record
    Since the UniProt ID isn't always 1:1 with a single protein, this class is used to represent a single protein entity
    as a subset of the UniProt ID if needed.
    """

    name: str = Field(
        ...,
        title="The name of the entity. Should be unique but not guaranteed to be so.",
    )
    uniprot_id: str = Field(..., title="The UniProt ID of the protein")
    sequence: str = Field(..., title="The protein sequence in FASTA format")
    start: int = Field(
        ...,
        title="The start position of the protein of interest in the UniProt sequence",
    )
    end: int = Field(
        ..., title="The end position of the protein of interest in the UniProt sequence"
    )

    @property
    def length(self) -> int:
        return self.end - self.start


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
