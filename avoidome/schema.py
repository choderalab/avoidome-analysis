from pydantic import BaseModel, Field
import numpy as np
from asapdiscovery.data.schema_v2.target import Target

class UniprotID(BaseModel):
    """
    A class to represent a uniprot entry
    """
    uniprot_id: str = Field(..., title="The UniProt ID of the protein")


class MolecularEntity(BaseModel):
    """
    A class to represent a single protein entity using the UniProt ID as a system of record
    Since the UniProt ID isn't always 1:1 with a single protein, this class is used to represent a single protein entity
    as a subset of the UniProt ID if needed.
    """
    name: str = Field(..., title="The name of the protein")
    uniprot_id: UniprotID = Field(..., title="The UniProt ID of the protein")
    sequence: str = Field(..., title="The protein sequence in FASTA format")
    uniprot_start: int = Field(..., title="The start position of the structure")
    uniprot_end: int = Field(..., title="The end position of the structure")

class StructureEntry(BaseModel):
    pass

class ResolvedChain(BaseModel):
    """
    A class to represent a resolved chain
    """
    chain_id: str = Field(..., title="The chain ID")
    start: int = Field(..., title="The start position of the chain")
    end: int = Field(..., title="The end position of the chain")

    @property
    def length(self) -> int:
        return self.end - self.start


class ExperimentalStructure(UniprotID):
    """
    A class to represent an experimental structure entry
    """
    pdb_id: str = Field(..., title="The PDB ID of the structure")
    method: str = Field(..., title="The method used to determine the structure")
    resolution: float = Field(..., title="The resolution (å) of the structure")
    resolved_chains: list[ResolvedChain] = Field(..., title="A list of resolved chains")

    @property
    def sequence_coverage(self):
        return np.mean([c.length for c in self.resolved_chains])


class PredictedStructure(UniprotID):
    """
    A class to represent a predicted structure entry
    """
    af_id: str = Field(..., title="The AlphaFold DB ID of the structure")
    uniprot_start: int = Field(..., title="The start position of the structure")
    uniprot_end: int = Field(..., title="The end position of the structure")
    confidence: float = Field(..., title="The confidence of the structure")


class LabeledTarget(Target):
    """
    A class to contain the info for a particular protein target, containing potentially multiple molecular entities
    """
    components: list[MolecularEntity] = Field(..., title="A list of molecular entities")

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
        return np.mean(np.array([s.sequence_coverage for s in self.experimental_structures]) / self.sequence_length)

    @property
    def average_confidence(self):
        return np.mean(np.array([s.confidence for s in self.predicted_structures]))

    @property
    def n_experimental_structures(self):
        return len([s for s in self.structures if isinstance(s, ExperimentalStructure)])

    @property
    def n_predicted_structures(self):
        return len([s for s in self.structures if isinstance(s, PredictedStructure)])