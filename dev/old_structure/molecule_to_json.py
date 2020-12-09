#
# # task: convert molecule to JSON and recreate molecule from JSON
#
# from pprint import pprint
# from src.molecule_classes_and_functions import Molecule, MoleculeList
#
# m1 = Molecule.from_smiles('CCOC(=O)[C@@H]1CCCN(C(=O)c2nc(-c3ccc(C)cc3)n3c2CCCCC3)C1')
# # m1.serialize_as_json()
# # pprint(m1.__dict__, width=600)
# df = m1.bond_distances_as_dataframe()
# # print(df)
#
#
# # for bond in df['atoms']:
# #     a1, a2 = bond
# #     print((a1, a1.element, a1.total_index), (a2, a2.element, a2.total_index))
# #     print(a1.charge, a1.radical, a1.stereochem)
#
#
# # create a class which generates a molecule from its json-repr:
# import json
#
#
# inchi = m1.inchi
# canonical_SMILES = m1.smiles
# chemical_formula = m1.chem_formula
# mol_weight = m1.get_molweight
# tot_ele = m1.get_total_electrons
# # atom_coords_origin = m1.coordinate_metadata[0]
#
#
# from rdkit import Chem
#
# def neutralize_atoms(mol):
#     pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
#     at_matches = mol.GetSubstructMatches(pattern)
#     at_matches_list = [y[0] for y in at_matches]
#     if len(at_matches_list) > 0:
#         for at_idx in at_matches_list:
#             atom = mol.GetAtomWithIdx(at_idx)
#             chg = atom.GetFormalCharge()
#             hcount = atom.GetTotalNumHs()
#             atom.SetFormalCharge(0)
#             atom.SetNumExplicitHs(hcount - chg)
#             atom.UpdatePropertyCache()
#     return mol


# smiList = ['CC(CNC[O-])[N+]([O-])=O',
#            'C[N+](C)(C)CCC([O-])=O',
#            '[O-]C1=CC=[N+]([O-])C=C1',
#            '[O-]CCCN=[N+]=[N-]',
#            'C[NH+](C)CC[S-]',
#            'CP([O-])(=O)OC[NH3+]']
# for m in smiList:
#     print(m)

# mols = [Chem.MolFromSmiles(m) for m in smiList]
# print(mols)

# neutralized_mols = []
# # Neutralize molecules by atom
# for mol in mols:
#     neutralize_atoms(mol)
#     print(Chem.MolToSmiles(mol))
#     neutralized_mols.append(mol)
#
# print(neutralized_mols)


# with open('../file_import/150_smiles-strings_example.smi') as fr:
#     smiles = fr.read()
#     for smiles in smiles.split('\n'):
#         rdkit_mol = Chem.MolFromSmiles(smiles)
#         rdkit_mol = neutralize_atoms(rdkit_mol)
#         neutral_smiles = Chem.MolToSmiles(rdkit_mol)
#         with open('../file_import/150_smiles_neutralized.smi', 'a') as fw:
#             fw.write(neutral_smiles + '\n')


# mol_lst = MoleculeList.from_smiles_file('../file_import/150_smiles_neutralized.smi')
# print(mol_lst)



import pandas as pd

# # molhashs:
# from rdkit.Chem import rdMolHash
# import rdkit
# molhashf = rdMolHash.HashFunction.names
# pprint(molhashf)
#
# # Generate MolHashes for molecule 's' with all defined hash functions.
# s = Chem.MolFromSmiles('CC(C(C1=CC(=C(C=C1)O)O)O)N(C)C(=O)OCC2=CC=CC=C2')
#
# for i, j in molhashf.items():
#     print(i, rdMolHash.MolHash(s, j))
#


# # inserting a document (=row)
# import datetime
# post = {"author": "Lois",
#         "text": "where am I?",
#         "tags": ["pandas", "scikit-learn"],
#         "date": datetime.datetime.utcnow()}
# print(datetime.datetime.utcnow())
#
# post_id = posts.insert_one(post).inserted_id
# print(post_id)  # 5f96ce45ed27ea266dfcc2d9

# print(db.list_collection_names())  # ['posts_collection']
#
#
# print(posts.count_documents({"author": "Mike"}))  # 2
# if more than 1 document is found: -> cursor object
# cursor = posts.find({"author": "Mike"})
# print(cursor)
# for post in cursor:
#     pprint(post)



#
# import pymongo
# from pprint import pprint
#
# # connects to default host/port:
# client = pymongo.MongoClient()
#
# # select a database (test)
# db = client['test-database']
#
# # select a collection (eq. to sql-table)
# posts = db['posts_collection']
# print(posts)
#
# # iterate over every document in the
# for post in posts.find():
#     pprint(post)
# print()
#
# # not_existing = posts.find_one({"author": "Steven"})['tags'].get('other_DFT')
# # print(not_existing)
# #
# #
# #
# # # short (adding a key-value pair to a nested dictionary):
# old_entry = posts.find_one({"author": "Steven"})['tags']
# pprint(old_entry)
# old_entry['DFT22'] = [2, 's', 0]
# result2 = posts.update_one({"author": "Lois"}, {'$set': {'tags': old_entry}})
# print(result2.matched_count)
# print(result2.modified_count)

# look at the newly generated document:
# found_post = posts.find_one({"author": "Lois"})
# pprint(found_post)

# m = Molecule.from_smiles('CCOC(=O)[C@@H]1CCCN(C(=O)c2nc(-c3ccc(C)cc3)n3c2CCCCC3)C1')
# # pprint(m.__dict__)
# pprint(m.atom_coordinates_matrix)





# import datetime
# molecule = {"InChI": "inchi",
#             "canonical_SMILES": "smiles",
#             "Molecular_Weight": 219.3321,
#             "Chemical_Formula": "C20H19",
#             "Total_Count_Electrons": 240,
#             "bonds": "where am I?",
#             "atoms": '',
#             "atom_coordinates": {'EKTG': {'meta_data': (),
#                                           'atomic coordinates': []},
#                                  'DFT': 'atomic coordinates'},
#             "date": datetime.datetime.utcnow()}
#
#
# m = Molecule.from_smiles('CCOC(=O)[C@@H]1CCCN(C(=O)c2nc(-c3ccc(C)cc3)n3c2CCCCC3)C1')
# # print(m.save_in_mongodb())
#
# atom_dict = dict(rdkit=[1, 2, 3])
# print(atom_dict)
# atom_dict['DFT'] = [2, 3, 4]
# print(atom_dict)

# ['total_index', 'Element', 'x', 'y', 'z']
# atom_coordinates = []
# for atom in m.contained_atoms:
#     x = [atom.total_index, atom.element, atom.x, atom.y, atom.z]
#     atom_coordinates.append(x)
# pprint(atom_coordinates)


# x = m.bond_properties_list
# pprint(x)
# m.visualize_as_3d_plot()

# testo = {'_id': 2,
#          'author': 'Steven',
#          'date': 3,
#          'tags': {'DFT': [2, 5, 1], 'EKTG': [1, 2, 3]},
#          'text': 'I send you a request!'}
#
# pprint(testo)
#
# testo['tags']['another-DFT'] = [8, 7, 6, 5, 4]
#
# pprint(testo)
#


# m = Molecule.from_smiles('CC')
#
# molblock = ''
# for atom in m.contained_atoms:
#     # print(atom.element, atom.x, atom.y, atom.z)
#     prr = f'{atom.element}    {atom.x}    {atom.y}    {atom.z}\n'
#     molblock += prr


#
# from rdkit.Chem import AllChem
# molecule = Chem.MolFromSmiles('CC')
# molecule = Chem.AddHs(molecule)
# AllChem.EmbedMolecule(molecule, randomSeed=0xf00d)
# molblock = Chem.MolToXYZBlock(molecule)
# print(molblock)
# bo = molblock.split('\n')[2:]
# print(bo)
# molblock = '\n'.join(bo)
#
#
# print(molblock)
#
#
# import uuid
#
# # add command to the molblock:
# command = '! B3LYP def2-SVP Opt\n'
#
# modfirst = '\n* xyz 0 1\n'
# modlast = '*'
#
# inputt = command + modfirst + molblock + modlast
#
# print(inputt)
#
# mol_name = str(uuid.uuid1())[:18]
# with open(f'../ORCA_data/{mol_name}.inp', 'w') as mol:
#     mol.write(inputt)
#
#
# from subprocess import Popen, PIPE
# import subprocess
# from pprint import pprint
#
# import time, os
# start_time1 = time.perf_counter()
#
#
# path = '../ORCA_data/'
#
# input_output_list = [f'{mol_name}']
# cmds_list = [['/home/niklas/orca/orca', f'{input_output_name}.inp', '>', f'{input_output_name}.out'] for input_output_name in input_output_list]
# procs_list = [Popen(cmd, stdout=PIPE, stderr=PIPE, cwd=path, shell=False) for cmd in cmds_list]  # cwd=path, , text=True
#
# for proc in procs_list:
#     print(f'args of this process (subprocess PID={proc.pid}, os PID={os.getpid()}): {proc.args}')
#     stdout, stderr = proc.communicate()
#
#     file_name = proc.args[-1][:-4]
#
#     # proc.wait()
#
#     with open(f'../ORCA_data/{file_name}.xyz', 'r') as read_file:
#         kaboom = read_file.read().split('\n')
#
#     subprocess.call([f'rm {file_name}*'], shell=True, cwd=path)
#     proc.wait()
#
# # extract data from xyz-files:
#
# def is_float(s):
#     try:
#         float(s)
#         return True
#     except ValueError:
#         return False
#
#
# print(kaboom)
#
# atom_enumerator, xyz_lst = 0, []
#
# for line in kaboom:
#     line = line.split()
#     if len(line) == 4 and len(line[0]) <= 2 and is_float(line[1]):
#         atom_enumerator += 1
#         element, x, y, z = line
#         print(atom_enumerator, element, x, y, z)
#         xyz_lst.append([atom_enumerator, element, float(x), float(y), float(z)])
#
# pprint(xyz_lst)  # use matrix to overwrite old coordinates with new ones
#
# # # --------------
# # save in mongodb!?
#
# print(f"Script execution finished after {(round(time.perf_counter()-start_time1, 3))} s")


# extract xyz-block and


# ----------------------------------------------------------------------------------------------------------
# complete:

# # calculate a single molecule:
# import uuid
# from subprocess import Popen, PIPE
# import subprocess
# from pprint import pprint
# import time
#
# m = Molecule.from_smiles('CC')
#
# molblock = ''
# for atom in m.contained_atoms:
#     molblock += f'{atom.element}    {atom.x}    {atom.y}    {atom.z}\n'
# print(molblock)
#
# # add command to the molblock:
# command = '! B3LYP def2-SVP Opt\n'
# modfirst = '\n* xyz 0 1\n'
# modlast = '*'
# inputt = command + modfirst + molblock + modlast
#
#
# mol_name = str(uuid.uuid1())[:18]
# with open(f'../ORCA_data/{mol_name}.inp', 'w') as mol:
#     mol.write(inputt)
#
# start_time1 = time.perf_counter()
# # calculate only a single molecule (maybe built into the Molecule class)
# path = '../ORCA_data/'
# cmd = ['/home/niklas/orca/orca', f'{mol_name}.inp', '>', f'{mol_name}.out']
# process = Popen(cmd, stdout=PIPE, stderr=PIPE, cwd=path, shell=False)
# stdout, stderr = process.communicate()
# file_name = process.args[-1][:-4]
# with open(f'../ORCA_data/{file_name}.xyz', 'r') as read_file:
#     kaboom = read_file.read().split('\n')
# subprocess.call([f'rm {file_name}*'], shell=True, cwd=path)
# process.wait()
#
#
# def is_float(s):
#     try:
#         float(s)
#         return True
#     except ValueError:
#         return False
#
#
# def xyz_extractor(xyz_string):
#     """ Extracts atomic coordinates of a xyz-string """
#     atom_enumerator, xyz_lst = 0, []
#     for line in xyz_string:
#         line = line.split()
#         if len(line) == 4 and len(line[0]) <= 2 and is_float(line[1]):
#             atom_enumerator += 1
#             element, x, y, z = line
#             print(atom_enumerator, element, x, y, z)
#             xyz_lst.append([atom_enumerator, element, float(x), float(y), float(z)])
#     return xyz_lst
#
#
# pprint(xyz_extractor(kaboom))
# print(f"Script execution finished after {(round(time.perf_counter()-start_time1, 3))} s")


# --------------------------------------------------------------------------------------------------------

# calculate a few molecules in parallel:
import uuid
from subprocess import Popen, PIPE
import subprocess
from pprint import pprint
import time, os
from rdkit.Chem import AllChem
from rdkit import Chem
#
#
# def add_command_to_molblock(molblock):
#     command = '! B3LYP def2-SVP Opt\n'
#     return command + '\n* xyz 0 1\n' + molblock + '*'
#
#
# smiles_lst = ['C', 'O', 'B']
#
# mol_lst = MoleculeList.from_smiles_list(smiles_lst)
#
# for num, molecule in enumerate(mol_lst):
#     print(molecule)
#     molblock = molecule.xyz_molblock
#     final_molblock = add_command_to_molblock(molblock)
#     print(final_molblock)
#
#     mol_name = str(uuid.uuid1())[:18]
#     with open(f'../ORCA_data/{mol_name}.inp', 'w') as mol:
#         mol.write(final_molblock)
#     pair = (num, mol_name, final_molblock)



# # Very small molecules used as testset for DFT Calculations in parallel:
# with open('../file_import/orca_canon_testset.smi', 'r') as f:
#     # filtering invalid SMILES strings out by canonicalizing it using RDKit:
#     canon_smiles_lst = [Chem.MolToSmiles(Chem.MolFromSmiles(mol)) for mol in f.read().split('\n') if Chem.MolFromSmiles(mol) is not None]
#     print(canon_smiles_lst)
#





# def smiles_to_preprocessed_molblock(smiles):
#     molecule = Chem.MolFromSmiles(smiles)
#     molecule = Chem.AddHs(molecule)
#     AllChem.EmbedMolecule(molecule, randomSeed=0xf00d)
#     molblock = Chem.MolToXYZBlock(molecule)
#     bo = molblock.split('\n')[2:]
#     molblock = '\n'.join(bo)
#     return '\n* xyz 0 1\n' + molblock + '*'


# for num, smiles in enumerate(smiles_lst):
#     molblock = smiles_to_preprocessed_molblock(smiles)
#     final_molblock = add_command_to_molblock(molblock)
#     print(final_molblock)
#     print(f'\n\n\n\n')
#     mol_name = str(uuid.uuid1())[:18]
#     with open(f'../ORCA_data/{mol_name}.inp', 'w') as mol:
#         mol.write(final_molblock)
#     pair = (num, mol_name, final_molblock)

# # add command to the molblock:
# command = '! B3LYP def2-SVP Opt\n'
# modfirst = '\n* xyz 0 1\n'
# modlast = '*'
# inputt = command + modfirst + molblock + modlast
# print(inputt)

# mol_name = str(uuid.uuid1())[:18]
# with open(f'../ORCA_data/{mol_name}.inp', 'w') as mol:
#     mol.write(inputt)

# start_time1 = time.perf_counter()
# path = '../ORCA_data/'
# input_output_list = [f'{mol_name}']
# cmds_list = [['/home/niklas/orca/orca', f'{input_output_name}.inp', '>', f'{input_output_name}.out'] for input_output_name in input_output_list]
# procs_list = [Popen(cmd, stdout=PIPE, stderr=PIPE, cwd=path, shell=False) for cmd in cmds_list]  # cwd=path, , text=True
#
# xyz_data = []
# for proc in procs_list:
#     print(f'args of this process (subprocess PID={proc.pid}, os PID={os.getpid()}): {proc.args}')
#
#     stdout, stderr = proc.communicate()
#     file_name = proc.args[-1][:-4]
#
#     with open(f'../ORCA_data/{file_name}.xyz', 'r') as read_file:
#         kaboom = read_file.read()
#     xyz_data.append(kaboom)
#     subprocess.call([f'rm {file_name}*'], shell=True, cwd=path)
#     proc.wait()
#
# print(xyz_data)
#
#
# def is_float(s):
#     try:
#         float(s)
#         return True
#     except ValueError:
#         return False
#
#
# def xyz_extractor(xyz_string):
#     """ Extracts atomic coordinates of a xyz-string """
#     atom_enumerator, xyz_lst = 0, []
#     for line in xyz_string.split('\n'):
#         line = line.split()
#         if len(line) == 4 and len(line[0]) <= 2 and is_float(line[1]):
#             atom_enumerator += 1
#             element, x, y, z = line
#             print(atom_enumerator, element, x, y, z)
#             xyz_lst.append([atom_enumerator, element, float(x), float(y), float(z)])
#     return xyz_lst
#
#
# for xyz_string in xyz_data:
#     extracted_lst = xyz_extractor(xyz_string)
#     pprint(extracted_lst)
#
# print(f"Script execution finished after {(round(time.perf_counter()-start_time1, 3))} s")
#

# # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# # # use multiprocessing:
# import multiprocessing
# from tqdm import tqdm
#
# # Very small molecules used as testset for DFT Calculations in parallel (invalid strings are removed, what about charg/proton changed stuff?):
# with open('../file_import/orca_canon_testset.smi', 'r') as f:  # ../file_import/orca_canon_testset.smi
#     # filtering invalid SMILES strings out by canonicalizing it using RDKit:
#     canon_smiles_lst = [Chem.MolToSmiles(Chem.MolFromSmiles(mol)) for mol in f.read().split('\n') if Chem.MolFromSmiles(mol) is not None and mol != '']
#     print(canon_smiles_lst)
#
# # mol_lst = MoleculeList.from_smiles_file('../file_import/orca_canon_testset.smi')
# # print(mol_lst)
#
# pool = multiprocessing.Pool()
#
# # idea: enumerate molecule with code, then list it and make multiprocessing
# result = pool.imap(Molecule.from_smiles, canon_smiles_lst)
#
# return_lst = []
# with tqdm(total=len(canon_smiles_lst)) as pbar:
#     for molecule in result:
#         return_lst.append(molecule)
#         pbar.update(1)
# print(return_lst)

# for num, mol in enumerate(return_lst):
#     print(mol)
#     print(mol.bond_properties_list)
#     print(mol.coordinate_metadata)
#     if num >= 2:
#         break






# mol_lst = MoleculeList.from_smiles_file('../file_import/orca_canon_testset.smi')
# print(mol_lst)
# mol_lst.save_in_mongodb()

# mol = Molecule.from_smiles('CC')
# print(mol.inchi)
# mol.save_in_mongodb()
# for


# # check all entries:
# import pymongo
# from pprint import pprint
# # connects to default host/port:
# client = pymongo.MongoClient()
# # select a database (test)
# db = client['test_database']
# # select a collection (eq. to sql-table)
# mols = db['molecules']
# print(mols)
# # iterate over every document in the
# for mol in mols.find():
#     pprint(mol)

# for mol in mols.find({'InChI': 'InChI=1S/CH5N/c1-2/h2H2,1H3'}):
#     pprint(mol)


# class MongoDB:
#     def __init__(self):
#         pass
#
#     def query_all(self):
#         pass
#
#     def


# # # iterate over every document in the
# for mol in mols.find({"InChI": "InChI=1S/H2O/h1H2"}):
#     pprint(mol)

# import pandas as pd
#
# mol = Molecule.from_mongo_db('InChI=1S/C3H6O/c1-3(2)4/h1-2H3')
# pprint(mol)
# print(type(mol))
#
# # molecule_entry['atoms']
# # molecule_entry['bonds']
# # molecule_entry['canonical_SMILES']
# # molecule_entry['InChI']
#
# empty_mol = Molecule()
#
# from src.molecule_classes_and_functions import Atom
# a2 = Atom('H', index=1, total_index=2)
# a3 = Atom('H', index=2, total_index=3)
# a4 = Atom('H', index=3, total_index=4)
# a1 = Atom('C', index=1, total_index=1)
# a5 = Atom('H', index=4, total_index=5)
# empty_mol.mol_graph[a1] = 'no bonds'
# empty_mol.mol_graph[a2] = 'no bonds'
# empty_mol.mol_graph[a3] = 'no bonds'
# empty_mol.mol_graph[a4] = 'no bonds'
# empty_mol.mol_graph[a5] = 'no bonds'

# pprint(empty_mol.mol_graph)
#
# for atom_index in empty_mol.mol_graph:
#     print(atom_index, atom_index.total_index)
#
# rd_mol = Chem.MolFromSmiles('CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC')
# # rd_mol = Chem.MolFromSmiles('CCCC')
# rd_mol = Chem.AddHs(rd_mol)
# AllChem.EmbedMolecule(rd_mol, randomSeed=0xf00d)
# rd_molbl = Chem.MolToMolBlock(rd_mol)
# print(rd_molbl)
# with open('superlong_chain.smi', 'w') as fw:
#     fw.write(rd_molbl)
# mol = Chem.MolFromMolBlock(rd_molbl)
# mol = Chem.MolToSmiles(mol)
# print(mol)
# test_molblock = """
#      RDKit          3D
#
#  10  9  0  0  0  0  0  0  0  0999 V2000
#    -1.2993   -0.0224   -0.1785 C   0  0  0  0  0  0  0  0  0  0  0  0
#     0.0124   -0.0562    0.5044 C   0  0  0  0  0  0  0  0  0  0  0  0
#     1.3154    0.0365   -0.1875 C   0  0  0  0  0  0  0  0  0  0  0  0
#     0.0085   -0.1681    1.7195 O   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.8782    0.8636    0.2104 H   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.2233    0.0609   -1.2707 H   0  0  0  0  0  0  0  0  0  0  0  0
#    -1.9194   -0.9067    0.0979 H   0  0  0  0  0  0  0  0  0  0  0  0
#     1.9956   -0.7005    0.3197 H   0  0  0  0  0  0  0  0  0  0  0  0
#     1.7423    1.0330    0.0523 H   0  0  0  0  0  0  0  0  0  0  0  0
#     1.2460   -0.1400   -1.2675 H   0  0  0  0  0  0  0  0  0  0  0  0
#   1  2  1  0
#   2  3  1  0
#   2  4  2  0
#   1  5  1  0
#   1  6  1  0
#   1  7  1  0
#   3  8  1  0
#   3  9  1  0
#   3 10  1  0
# M  END
#
# """
# mol = Chem.MolFromMolBlock(test_molblock)
# mol = Chem.MolToSmiles(mol)
# print(mol)

# with open('MolToXYZBlock', 'w') as fw:
#     fw.write(Chem.MolToXYZBlock(rd_mol))


# with open('MolToFASTA', 'w') as fw:
#     fw.write(Chem.MolToFASTA(rd_mol))


# import json
# with open('MolToJSON', 'r') as f:
#     kaboom = f.read()
#     loaded_json = json.loads(kaboom)
#     with open('MolToJSON_new', 'w') as fw:
#         fw.write(json.dumps(loaded_json, indent=4, sort_keys=True))

# import json
# loaded_json = json.loads(Chem.MolToJSON(rd_mol))
# with open('MolToJSON_newwwww', 'w') as fw:
#     fw.write(json.dumps(loaded_json, indent=4, sort_keys=True))
#
# with open('MolToPDBBlock', 'w') as fw:
#     fw.write(Chem.MolToPDBBlock(rd_mol))
#
# with open('MolToV3KMolBlock', 'w') as fw:
#     fw.write(Chem.MolToV3KMolBlock(rd_mol))

# filter()

#
# dict1 = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5}
# # Double each value in the dictionary
# double_dict1 = {k*2: v*2 for (k, v) in dict1.items()}
# print(double_dict1)
# # Check for items greater than 2 (conditionals)
# dict1_cond = {k: v for (k, v) in dict1.items() if v > 2}
# print(dict1_cond)
# dict1_doubleCond = {k: v for (k, v) in dict1.items() if v > 2 if v % 2 == 0}
# print(dict1_doubleCond)
#
#
# # Initialize `fahrenheit` dictionary
# fahrenheit = {'t1': -30, 't2': -20, 't3': -10, 't4': 0}
# # easier: use comp. dict.:
# cels_dict = {k: round(((v-32) * 5 / 9), 2) for (k, v) in fahrenheit.items()}
# print('fahr dict conv to cels dict:', cels_dict)
#
#
#
# dict1 = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5, 'f':6}
#
# # Identify odd and even entries
# dict1_tripleCond = {k: ('even' if v % 2 == 0 else 'odd') for (k, v) in dict1.items()}
# print(dict1_tripleCond)
#
# # inner comp dicts --> can be used for working with graphs !
# nested_dict = {'first': {'a': 1}, 'second': {'b': 2}}
# float_dict = {outer_k: {float(inner_v) for (inner_k, inner_v) in outer_v.items()} for (outer_k, outer_v) in nested_dict.items()}
# print(float_dict)
#
# str_nums = ["1", "2", "3", "4", "5"]
# str_conv = [int(x) for x in str_nums]
# print(str_conv)
#
# numbers = [-2, -1, 0, 1, 2, 3]
# abs_nums = [abs(x) for x in numbers]
# print(abs_nums)
# pom = [(x**2, x**3) for x in numbers]
# print(pom)
#
# some_text = ['I', 'am', 'Writing']
# lower_text = [x.lower() for x in some_text]
# print(lower_text)
#
# str_nums = ["1,3", "2,1", "3,", "4,97", "5,0"]
# str_conv = [float(x.replace(',', '.')) for x in str_nums]
# print(str_conv)
#
# # incredible strong module !
# import re
#
# more_text = ['hello;', 'my friend!', 'how', 'are you?']
# proc_text = [re.sub(r'[!?.:;,"()-]', "", word) for word in more_text]
# print(proc_text)
#
# # temp conversion:
# temps_in_fahrenheit = [75, 79, 64]
# temps_in_celsius = [round(((x-32) * 5 / 9), 2) for x in temps_in_fahrenheit]
# print('temps in °C:', temps_in_celsius)
#
# # inner comp dicts --> can be used for working with graphs !
# nested_dict = {'first': {'a': 1}, 'second': {'b': 2}}
# float_dict = {outer_k: {float(inner_v) for (inner_k, inner_v) in outer_v.items()} for (outer_k, outer_v) in nested_dict.items()}
# print(float_dict)
#
# import multiprocessing
# pool = multiprocessing.Pool()
# def square(n):
#     return n * n
#
#
# nums = [1, 2, 3]
# result = pool.imap(lambda x: square(x), nums)
# for num in result:
#     print(num)
#




# import multiprocessing
#
# try:
#     cpus = multiprocessing.cpu_count()
#     print(cpus)
# except NotImplementedError:
#     cpus = 2   # arbitrary default
#
#
# def square(n):
#     return n * n
#
# start_time = time.perf_counter()
# pool = multiprocessing.Pool(processes=cpus)  # nice: number of cpu cores/threads can be defined!
# pool.map(square, range(100))  # .imap creates an iterator, 1.209 s
#
# # normal .map() function -> lambda-function can be used (not with multiproc. map!)
# # list(map(lambda x: x * x, list(range(def_range))))   # 1.073
# print('elapsed time:', round(time.perf_counter()-start_time, 3))



# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# from queue import PriorityQueue
# from multiprocessing import Queue as MultiprocQueue
# from queue import Queue
# q = Queue()
# q.put("1")
# q.put("2")
# q.put("3")
# q.put("4")
# print(q)
# # print(q.get())
# # print(q.get())
# # print(q.get())
# # print(q.get())
# while not q.empty():
#     next_item = q.get()
#     print(next_item)


# # Queue + scheduler + generator + multiprocessing??

# from sched import scheduler
# from time import sleep
# pprint(dir(scheduler))
#
#
# def job_function(a_string):
#     sleep(1)
#     return len(a_string)
#
# import sched
#
#
# test_list = ('dgvdv', 'dsfg', 'f', 'sds', '34rf', 'dbgd', 'f', 'ss', 'sef')
# ma_gen = (e for e in test_list)
# for e in ma_gen:
#     print(e)

# import sched
# import time
# # instance is created
# scheduler = sched.scheduler(time.time, time.sleep)
# # function to print time
# # and name of the event
# def print_event(name):
#     print('EVENT:', time.time(), name)
# # printing starting time
# print('START:', time.time())
# # first event with delay of
# # 1 second
# e1 = scheduler.enter(1, 1, print_event, ('1 st',))
# # second event with delay of
# # 2 seconds
# e1 = scheduler.enter(2, 1, print_event, (' 2nd',))
# # executing the events
# scheduler.run()
#
# import datetime
# print(datetime.datetime.now())   # 2020-11-02 15:27:14.159020
#
# # library imported
# import sched
# import time
# # instance is created
# scheduler = sched.scheduler(time.time, time.sleep)
# # function to print time and
# # name of the event
# def print_event(name):
#     print('EVENT:', time.time(), name)
# # printing starting time
# print('START:', time.time())
# # event x with delay of 1 second
# # enters queue using enterabs method
# e_x = scheduler.enterabs(time.time() + 1, 1, print_event, argument=("Event X",))
# # executing the events
# scheduler.run()



# # retrieving a json containing data about bitcoing value:
# url = 'https://api.coindesk.com/v1/bpi/currentprice.json'
# import requests
# import json
# while True:
#     url_request = requests.get(url)
#     bitcoin_data = json.loads(url_request.text)
#     bitcoin_euro_value = bitcoin_data['bpi']['EUR']['rate_float']
#     print(round(float(bitcoin_euro_value), 2))
#     sleep(30)  # value is updated every 30s



# from apscheduler.schedulers.background import BackgroundScheduler
#
#
# scheduler = BackgroundScheduler()
# scheduler.add_job()

# # convert timestamp to datetime
# from datetime import datetime
# timestamp = 1545730073
# dt_object = datetime.fromtimestamp(timestamp)
# print("dt_object =", dt_object)
# print("type(dt_object) =", type(dt_object))



# # clean?
# with open('../file_import/2000_smi.smi', 'r') as f:  # ../file_import/orca_canon_testset.smi
#
#     canon_smiles_lst = [Chem.MolToSmiles(Chem.MolFromInchi(Chem.MolToInchi(neutralize_atoms((Chem.MolFromSmiles(mol)))))) for mol in f.read().split('\n') if Chem.MolFromSmiles(mol) is not None and mol != '']
#
#     #canon_smiles_lst = [Chem.MolToSmiles(Chem.MolFromSmiles(mol)) for mol in f.read().split('\n') if Chem.MolFromSmiles(mol) is not None and mol != '']
#     print(canon_smiles_lst)
#     canon_smiles_lst
#     joined_smiles = '\n'.join(canon_smiles_lst)
#     print(joined_smiles)
#     with open('../file_import/cleaned_2000_smi.smi', 'w') as fw:
#         fw.write(joined_smiles)

# -----------

# testing behaviour of dictionaries (1st: __contains__ dunder-method)
test_dict = {'a_out1': {'a_in1': 1, 'a_in2': 2}, 'a_out2': {'a_in3': 3, 'a_in4': 4}}

# no iteration necessary for outter dict:
if 'a_out1' in test_dict.values():
    print('a_out1 was found using')

# iteration through inner dict for 'contains':
[print('kaboom') for x in test_dict.values() if 'a_in1' in x]

# ------
from src.molecule_classes_and_functions import Atom, Molecule
from collections import OrderedDict
# a1 = Atom('C', index=1, total_index=1)
# a2 = Atom('H', index=1, total_index=2)
# a3 = Atom('H', index=2, total_index=3)
#
# m1 = Molecule()
# m1.add_atom(a1)
# m1.add_atom(a3)
# m1.add_atom(a2)
# m1.add_bond(a1, a2, 's')
# m1.add_bond(a1, a3, 's')
# # print(m1.mol_graph)
# pprint(m1.__dict__)
# # hydro = m1.select_atom_by_total_index(2)
# # print(hydro)

import networkx
networkx.Graph()

methane = Molecule.from_smiles('CC')
print(methane)
pprint(methane.__dict__)



# hydro = methane.select_atom_by_total_index(3)
# print(hydro)



# class TestClass:
#     def __init__(self, name):
#         self.name = name
#
# p = TestClass('boy')
# print(p.name)
# # delattr(p, 'name')
# del p.name
# # print(p.name)








