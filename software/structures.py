from pydantic import BaseModel, Field
import numpy as np


class StructureEntry(BaseModel):
    """
    A class to represent a structure entry
    """
    uniprot_id: str = Field(..., title="The UniProt ID of the protein")


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


class ExperimentalStructure(StructureEntry):
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

class PredictedStructure(StructureEntry):
    """
    A class to represent a predicted structure entry
    """
    af_id: str = Field(..., title="The AlphaFold DB ID of the structure")
    uniprot_start: int = Field(..., title="The start position of the structure")
    uniprot_end: int = Field(..., title="The end position of the structure")
    confidence: float = Field(..., title="The confidence of the structure")