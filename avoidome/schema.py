from pydantic import BaseModel, Field


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
