
import multiprocessing
import pickle
from collections import UserList

from .molecule import Molecule


class MoleculeList(UserList):
    """An especially for molecules constructed list. The goal is to add exceptional well sorting
    and filtering algorithms for molecules. the class should contain methods which are capable of creating
    more complex data structures than just a linear set!"""
    def __init__(self, data):
        super().__init__(data)
        # if data is list:
        #     print(f'initialized MoleculeList was filled with another list!')
        #     self.data.append(data)
        #     print(self.data)
        # # elif data is not None:
        # #     self.data = list(data)
        # else:
        #     print('Empty MoleculeList initialized.')
        #     self.data = list()

    # def __repr__(self):
    #     return self.data

    def sort_function1(self):
        new_lst = []
        for x in reversed(self.data):
            new_lst.append(x)
        self.data = new_lst

    @classmethod
    def from_smiles_file(cls, smiles_path):
        with open(f"{smiles_path}", "r") as read_file:
            smiles_list = read_file.read().splitlines()
            pool = multiprocessing.Pool()
            result = pool.imap(Molecule.from_smiles, smiles_list)

            total = len(smiles_list)
            dummy_lst = []
            print("SMILES-file is converted into a List of Molecules...")
            for molecule in result:
                dummy_lst.append(molecule)

                # [12:09:56] WARNING: Omitted undefined stereo occurs
                # [12:09:56] WARNING: Charges were rearranged

            # with tqdm(total=total) as pbar:
            #     for molecule in result:
            #         dummy_lst.append(molecule)
            #         pbar.update(1)
        print("...MoleculeList is created!")
        return cls(data=dummy_lst)

    @classmethod
    def from_smiles_list(cls, smiles_lst):
        """ generates a list of molecues using a list of smiles-strings """
        pool = multiprocessing.Pool()
        result = pool.imap(Molecule.from_smiles, smiles_lst)
        total = len(smiles_lst)
        dummy_lst = []
        print("SMILES-file is converted into a List of Molecules...")
        for molecule in result:
            dummy_lst.append(molecule)
        print("...MoleculeList is created!")
        return cls(data=dummy_lst)

    # @classmethod
    # def from_smiles_file(cls, smiles_path):
    #     with open(f"{smiles_path}", 'r') as read_file:
    #         smiles_list = read_file.read().splitlines()
    #
    #         mol_lst = []
    #         for smiles in smiles_list:
    #             print(smiles)
    #             mol_lst.append(Molecule.from_smiles(smiles))
    #             print()
    #             # [Molecule.from_smiles(x) for x in smiles_list]
    #     return cls(data=mol_lst)

    @classmethod
    def from_sdf_file(self):
        # goes to the Molecule.from_molblock and
        # applies multiprocessing
        pass

    @classmethod
    def load_from_pickle(cls, pickle_file):
        # if 'pickle_storage' not in os.getcwd():
        #     os.chdir('../pickle_storage')
        with open(f"{pickle_file}", "rb") as molecule_file:
            molecule_list = pickle.load(molecule_file)
        return cls(molecule_list)

    def save_as_pickle(self, pickle_name):
        # if 'pickle_storage' not in os.getcwd():
        #     os.chdir('../pickle_storage')
        with open(f"{pickle_name}", "wb") as molecule_file:
            pickle.dump(self.data, molecule_file)

    def save_in_mongodb(self):
        """ saves all molecules in the database """
        for mol in self.data:
            mol.save_in_mongodb()

    @staticmethod
    def convert_to_dataframe(self):
        """convert contained information into a pandas dataframe
        for further data exploration"""
        print(f"MoleculeList is converted into a Dataframe"
              f"structure: SMILES, Atomweight, electroncount")
        # DataFrame()
        pass


