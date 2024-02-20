
class Target(BaseModel):
    """
    A class to represent a protein target
    """
    name: str = Field(..., title="The name of the protein")
    uniprot_id: str = Field(..., title="The UniProt ID of the protein")
    sequence: str = Field(..., title="The protein sequence")
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