import requests
import itertools
import os
import glob
import pubchempy as pcp

from rcsbapi.search import TextQuery


# Module-level print flag - set from modrag.py
print_flag = False

def smiles_node(names_list: list[str]) -> (str):
  '''
    Queries Pubchem for the smiles string of the molecule based on the name.
      Args:
        names_list: the list of molecule names
      Returns:  
        smiles_string: a string of the tool results
  '''
  print("smiles tool")
  print('===================================================')

  smiles_string = ''
  for name in names_list:
    try:
        res = pcp.get_compounds(name, "name")
        smiles = res[0].smiles
        #smiles = smiles.replace('#','~')
        smiles_string += f'{name}: The SMILES string for the molecule is: {smiles}\n'
    except:
        smiles = "unknown"
        smiles_string += f'{name}: Fail\n'

  return smiles_string

def get_protein_from_pdb(pdb_id: str, protein_name: str) -> str:
  '''
    Helper function to get the protein information from the PDB database.
    Args:
      pdb_id: the PDB ID of the protein
      protein_name: the name of the protein
    Returns:
      r.text: the PDB information as a string
  '''
  print('PDB retrieval tool')
  print('===================================================')

  # Check whether a .pdb file for this PDB ID is already present in pdb_files/.
  # The leading part of the filename may differ (e.g. "sult1a3" vs "sulfotransferase"),
  # so only the PDB ID is used as the match test.
  pdb_id_upper = pdb_id.upper()
  for existing in glob.glob('pdb_files/*.pdb'):
    stem = os.path.basename(existing).rsplit('.', 1)[0]
    if stem.upper().endswith(f'_{pdb_id_upper}'):
      print(f'PDB file for PDB ID {pdb_id} already present at {existing}; reusing it.')
      return (f'The PDB file for protein {protein_name} with PDB ID {pdb_id} was already '
              f'present in pdb_files/ as {existing} and has been reused.')

  url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
  r = requests.get(url)
  with open(f"pdb_files/{protein_name}_{pdb_id}.pdb", 'w') as f:
    f.write(r.text)

  return f'The PDB file for protein {protein_name} with PDB ID {pdb_id} has been retrieved and saved as pdb_files/{protein_name}_{pdb_id}.pdb'

def find_PDBID_node(test_protein_list: list[str]) -> str:
  '''
    Accepts a protein name and searches the protein databank for PDB IDs that match along with the entry titles.
      Args:
        test_protein_list: the protein names to query
      Returns:
        pdb_string: a string containing the results of the PDB search.
  '''

  print(f"PDB search tool")
  print('===================================================')

  pdb_string = ''
  which_pdbs = 0

  for test_protein in test_protein_list:
    try:
      query = TextQuery(value=test_protein)
      results = query()

      def pdb_gen():
        for rid in results:
          yield(rid)

      take10 = itertools.islice(pdb_gen(), which_pdbs, which_pdbs+10, 1)

      pdb_string += f'10 PDBs that match the protein {test_protein} are: \n'
      for pdb in take10:
        data = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}").json()
        title = data['struct']['title']
        pdb_string += f'PDB ID: {pdb}, with title: {title} \n'
    except:
      pdb_string += f'Failed to get PDB IDs for protein {test_protein}\n'

  return pdb_string