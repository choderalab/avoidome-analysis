"""
This module contains functions to request and parse data from the UniProt API
"""
import requests, sys

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


def request_uniprot(uniprot_id):
    """
    A function to request a protein entry from the UniProt API
    """
    requestURL = f"https://www.ebi.ac.uk/proteins/api/proteins/{uniprot_id}"
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    return r.json()



def parse_resolved_chains(chain_str):
    """
    A function to parse resolved chains for a particular uniprot id from a Uniprot string
    """
    chain_letters = chain_str.split('=')[0].split('/')
    start, end = chain_str.split('=')[1].split('-')

    resolved_chains = []
    for chain in chain_letters:
        resolved_chains.append(ResolvedChain(chain_id=chain, start=int(start), end=int(end)))
    return resolved_chains


def parse_uniprot_accession(target_name: str, uniprot_id: str, uniprot_dict: dict) -> Target:
    """
    A function to parse the UniProt accession from a UniProt API response
    """
    pdb_ids = [ref for ref in uniprot_dict['dbReferences'] if ref['type'] == 'PDB']
    af_ids = [ref for ref in uniprot_dict['dbReferences'] if ref['type'] == 'AlphaFoldDB']

    refs = []
    for ref in pdb_ids:
        properties = ref['properties']
        refs.append(ExperimentalStructure(
            uniprot_id=uniprot_id,
            pdb_id=ref['id'],
            method=properties['method'],
            resolution=properties['resolution'].split(' ')[0],
            resolved_chains=parse_resolved_chains(properties['chains'])
        ))
    for ref in af_ids:
        r = request_alphafold(uniprot_id)
        summary_data = r['structures'][0]['summary']
        refs.append(PredictedStructure(
            uniprot_id=uniprot_id,
            af_id=summary_data['model_identifier'],
            uniprot_start=summary_data['uniprot_start'],
            uniprot_end=summary_data['uniprot_end'],
            confidence=summary_data['confidence_avg_local_score']

        ))
    return Target(
        name=target_name,
        uniprot_id=uniprot_id,
        sequence=uniprot_dict['sequence']['sequence'],
        structures=refs
    )