"""
This module contains functions to request and parse data from the UniProt API
"""

import requests, sys
from pydantic import BaseModel, Field
import asapdiscovery.data.schema_v2.schema_base as schema_base
from avoidome.schema import (
    ProteinEntity,
)
from avoidome.structures import (
    StructureEntry,
    ExperimentalStructure,
    PredictedStructure,
)


def request_uniprot(uniprot_id):
    """
    A function to request a protein entry from the UniProt API
    """
    requestURL = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    r = requests.get(requestURL, headers={"Accept": "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r.json()


def request_alphafold(uniprot_id):
    """
    A function to request a protein entry from the UniProt API
    """
    requestURL = f"https://alphafold.ebi.ac.uk/api/uniprot/summary/{uniprot_id}.json"
    r = requests.get(requestURL)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r.json()


def download_alphafold_structure(model_url):
    """
    A function to download a predicted structure from the AlphaFold API
    """
    r = requests.get(model_url)
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r.content


class UniprotEntry(schema_base.DataModelAbstractBase):
    """
    A class to represent a uniprot entry
    """

    uniprot_id: str = Field(..., title="The UniProt ID of the protein")
    data: dict = Field(..., title="The data associated with the UniProt ID")

    def __str__(self):
        return self.uniprot_id

    @property
    def sequence(self):
        return self.data["sequence"]["sequence"]

    @property
    def sequence_length(self):
        return len(self.sequence)

    @property
    def name(self):
        return self.data["id"]

    @classmethod
    def from_uniprot_id(cls, uniprot_id: str):
        return cls(uniprot_id=uniprot_id, data=request_uniprot(uniprot_id))

    def get_experimental_structures(self) -> list[ExperimentalStructure]:
        """
        A function to get the experimental structures from the UniProt entry
        :return:
        """
        experimental_structures = []
        for ref in self.data["dbReferences"]:
            if ref["type"] == "PDB" and ref["properties"].get("resolution"):
                properties = ref["properties"]
                chain_str = properties["chains"]

                components = []
                for entry in chain_str.split(", "):
                    chain_letters = entry.split("=")[0].split("/")
                    start, end = entry.split("=")[1].split("-")
                    for chain in chain_letters:
                        components.append(
                            ProteinEntity(
                                name=self.name,
                                uniprot_id=self.uniprot_id,
                                sequence=self.sequence,
                                start=int(start),
                                end=int(end),
                            )
                        )
                try:
                    experimental_structures.append(
                        ExperimentalStructure(
                            pdb_id=ref["id"],
                            method=properties["method"],
                            resolution=properties["resolution"].split(" ")[0],
                            components=components,
                        )
                    )
                except KeyError:
                    print(f"Error parsing {ref}")
        return experimental_structures

    def get_alphafold_structures(self):
        af_ids = [
            ref for ref in self.data["dbReferences"] if ref["type"] == "AlphaFoldDB"
        ]
        r = request_alphafold(self.uniprot_id)
        summary_data = r["structures"][0]["summary"]
        return [
            PredictedStructure(
                af_id=summary_data["model_identifier"],
                confidence=summary_data["confidence_avg_local_score"],
                model_url=summary_data["model_url"],
                components=[
                    ProteinEntity(
                        name=self.name,
                        uniprot_id=self.uniprot_id,
                        sequence=self.sequence,
                        start=summary_data["uniprot_start"],
                        end=summary_data["uniprot_end"],
                    )
                ],
            )
        ]
