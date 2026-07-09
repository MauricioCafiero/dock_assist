import requests
import itertools
import os
import glob
import math
import pubchempy as pcp

from rcsbapi.search import TextQuery
from rdkit import Chem


# Module-level print flag - set from modrag.py
print_flag = False

# Distance (Angstrom) under which a co-crystallized molecule is reported as
# "close" to the docked ligand (minimum atom-to-atom distance).
NEARBY_DISTANCE_CUTOFF = 5.0

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

      take10 = itertools.islice(pdb_gen(), which_pdbs, which_pdbs+20, 1)

      pdb_string += f'20 PDBs that match the protein {test_protein} are: \n'
      for pdb in take10:
        data = requests.get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}").json()
        title = data['struct']['title']
        pdb_string += f'PDB ID: {pdb}, with title: {title} \n'
    except:
      pdb_string += f'Failed to get PDB IDs for protein {test_protein}\n'

  return pdb_string

def check_nearby_molecules(pdb_filepath: str, ligand_filepath: str) -> str:
  '''
    Checks for nearby molecules in the PDB file to asses if the docking has 
    located the correct binding site. Should be called to verify blind docking
    in the case where a ligand is present in the crystal structure. 

    Args:
      pdb_filepath: the path to the PDB file
      ligand_filepath: the path to the ligand file
    Returns:
      nearby_molecules: a string containing the results of the check.
  '''
  print(f"Nearby molecules check tool")
  print('===================================================')

  with open(pdb_filepath, 'r') as pdb_file:
      pdb_content = pdb_file.readlines()

  ligand_names = {}
  for line in pdb_content:
    if line.startswith('HETNAM'):
      molecule_symbol = line[11:15].strip()
      molecule_name = line[15:70].strip()
      if not molecule_symbol:
        continue
      if molecule_symbol not in ligand_names:
        ligand_names[molecule_symbol] = molecule_name
      elif molecule_name:
        ligand_names[molecule_symbol] += ' ' + molecule_name

  molecule_dict = {}
  for line in pdb_content:
    if line.startswith('HETATM'):
      molecule_symbol = line[17:20].strip()
      if molecule_symbol not in ligand_names:
        continue
      chain_id = line[21].strip()
      occ_key = (molecule_symbol, chain_id)
      if occ_key not in molecule_dict:
        molecule_dict[occ_key] = {'name': ligand_names[molecule_symbol] or molecule_symbol,
                                  'coords': []}
      x_coord = float(line[30:38])
      y_coord = float(line[38:46])
      z_coord = float(line[46:54])
      molecule_dict[occ_key]['coords'].append((x_coord, y_coord, z_coord))

  ligand_atom_coords = []
  try:
    supplier = Chem.SDMolSupplier(ligand_filepath, removeHs=True)
    poses = [m for m in supplier if m is not None]
  except Exception:
    poses = []
  for mol in poses:
    if not ligand_atom_coords:
      conf = mol.GetConformer()
      ligand_atom_coords = [(conf.GetAtomPosition(j).x,
                             conf.GetAtomPosition(j).y,
                             conf.GetAtomPosition(j).z) for j in range(mol.GetNumAtoms())]

  nearby_molecules = f'Checked for nearby molecules in {pdb_filepath}.\n'
  if not ligand_atom_coords:
    nearby_molecules += 'No ligand could be read from the SDF.\n'
  else:
    distances = []
    for (molecule_symbol, chain_id), info in molecule_dict.items():
      mol_coords = info['coords']
      if not mol_coords or not ligand_atom_coords:
        continue
      best = None
      for lx, ly, lz in ligand_atom_coords:
        for mx, my, mz in mol_coords:
          d = math.sqrt((lx - mx) ** 2 + (ly - my) ** 2 + (lz - mz) ** 2)
          if best is None or d < best:
            best = d
      distances.append((best, molecule_symbol, chain_id, info))
    distances.sort(key=lambda d: d[0])

    nearby_molecules += (f'Molecules within {NEARBY_DISTANCE_CUTOFF:.1f} A '
                         f'of the docked ligand (min atom-to-atom distance):\n')
    any_close = False
    for dist, molecule_symbol, chain_id, info in distances:
      molecule_name = info['name']
      line = (f'  {molecule_name} ({molecule_symbol}, chain {chain_id}): '
             f'min atom-to-atom {dist:.2f} A')
      if dist <= NEARBY_DISTANCE_CUTOFF:
        any_close = True
        nearby_molecules += line + ' -- CLOSE\n'
    if not any_close:
      nearby_molecules += '  (none found)\n'

  return nearby_molecules