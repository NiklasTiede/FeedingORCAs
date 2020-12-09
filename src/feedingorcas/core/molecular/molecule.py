"""
Molecule
========

"""

import datetime
import time
from bisect import insort
from math import sqrt

import matplotlib.pyplot as plt
import numpy as np
import pymongo
import rdkit
from pandas import DataFrame
from rdkit import Chem
from rdkit.Chem import AllChem

from ..utilities import AtomNotFoundError, annotate3d, neutralize_atoms
from .atom import Atom


class Molecule(object):
    """The Molecule class let's you create molecules. A Molecule-instance connects atoms in a graph-structure with
    each other. Each atoms cartesian coordinates, stereochemical information, charge and properties derived from its
    element (atomic number, weight, valence electrons) are stored within the atom-objects (see Atom-class). The graph
    representation of each molecule is unique and unambiguously. Individual Molecules are created by extracting their
    data from chemical format files (.smi, .mol, .pdb etc) or a string representation (SMILES, InChI etc.)."""

    # element symbol order for molecule representation in chemical formula (according to Hill notation)
    chem_form_seq = ("C", "H", "B", "Br", "Cl", "F", "I", "N", "O", "P", "S",
                     "Si")

    def __init__(self, smiles="", inchi=""):
        """The attribute fields of an initialized molecule are empty. They are filled by adding instantiated atoms
        and their connectivity information to the molecule. Thereby the molecular graph, the chemical formula
        (according to Hill notation) and a list of contained atoms are generated.

        Arguments:
            - smiles ('str'):  SMILES-representation of the molecule                (example: 'C=CC')

        """

        self.inchi = inchi
        self.smiles = smiles
        self.contained_atoms = list()
        self.mol_graph = dict()
        self.element_counter_index = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        self.chem_formula = ""

        # 5-tuple: [software, version, method, calc.-time, datetime]
        self.coordinate_metadata = [None, None, None, None, None]

    def __len__(self):
        """ The len-Method returns the number of atoms within the molecule. """
        return len(self.mol_graph)

    def __str__(self):
        """A single molecule is represented as its chemical formula in this format: 'mol(<num atoms>):chem_formula'
        print(m1) -> m(41):C19H17F2N1O2"""
        return f"m({len(self.mol_graph)}):{self.chem_formula}"

    def __repr__(self):
        """Representation of molecules within a collection (see __str__).
        print(mol_lst) -> [mol(34):C13H19N1O1, mol(26):C9H14N3, mol(44):C15H22O5S2]"""
        return f"m({len(self.mol_graph)}):{self.chem_formula}"

    def __iter__(self):
        """ The atoms within the molecule are iterable """
        return iter(self.mol_graph)

    def __next__(self):
        try:
            return next(self.mol_graph)
        except AtomNotFoundError as anfe:
            print(anfe)

    def add_atom(self, atom_instance):
        """The .add_atom-method appends atoms to the contained_atoms-list. Thereby each atom is indexed and each element
        counted to generate a chemical formula which is used for proper representation (__str__/__repr__)."""

        # a list of all atoms contained within the molecules (order equates to .sdf import reading)
        if atom_instance not in self.contained_atoms:
            self.contained_atoms.append(atom_instance)

        # to calculate chemical formula and count each atom-type to represent atoms as C1 C2 H1 H2 ...
        self.element_counter_index[atom_instance.atoms_internal_number] += 1

        # add index to atom, for their __str__/__repr__-methods:
        atom_instance.index = self.element_counter_index[
            atom_instance.atoms_internal_number]

        # update chemical formula (order of repr according to Hill notation)
        self.chem_formula = ""
        chem_formula_seq = (
            self.element_counter_index[2],
            self.element_counter_index[0],
            self.element_counter_index[1],
            self.element_counter_index[10],
            self.element_counter_index[9],
            self.element_counter_index[5],
            self.element_counter_index[11],
            self.element_counter_index[3],
            self.element_counter_index[4],
            self.element_counter_index[7],
            self.element_counter_index[8],
            self.element_counter_index[6],
        )

        for num, element_symbol in enumerate(chem_formula_seq):
            # exclude all non-existing elements
            if element_symbol > 0:
                self.chem_formula += Molecule.chem_form_seq[num] + str(
                    element_symbol)

    def add_bond(self, atom1, atom2, bond_order):
        """The .add_bond-method creates the molecular graph of the molecule (data structure: dict in a dict)
        bond order: s = single, d = double, t = triple. mol_graph is created independently from self.atom_contained."""

        # inserts atom-instances into the "first level" of the mol_graph (which is a dict in a dict)
        if (atom1 not in self.mol_graph
                and atom1 in self.contained_atoms):  # self.mol_graph
            self.mol_graph[atom1] = dict()
        if atom2 not in self.mol_graph and atom2 in self.contained_atoms:
            self.mol_graph[atom2] = dict()

        # inserts the connected atoms (and their bond order) into the "second level" of the mol_graph
        self.mol_graph[atom2][atom1] = bond_order
        self.mol_graph[atom1][atom2] = bond_order

    def molblock_data_extraction(self, molblock_input, rdkit_molblock=False):
        """A helper method which extracts molecular information (atom types, coordinates, connectivity, bond-type
        stereochemistry, charge, radical) from a RDKit generated molblock. By default it is assumed that the given
        molblock is not created by rdkit, so it is converted to a molecule and then back to a molblock, to stanardize
        the input for the molecule-data-handler"""

        # molblock is standardized using RDKit:
        if not rdkit_molblock:
            rdkit_mol = Chem.MolFromMolBlock(molblock_input)
            molblock = Chem.MolToMolBlock(rdkit_mol)
        else:
            molblock = molblock_input

        # extracting the molblock information:
        splitted_molblock = molblock.split("\n")
        atom_enumerator = 0
        for line in splitted_molblock:
            line_list = line.split()

            # Atom-block information is extracted (stereoatom_info represents an atoms charge):
            if len(line_list) == 16:
                atom_enumerator += 1  # atoms are enumerated starting at 1
                charge, radical = 0, False
                x, y, z = float(line_list[0]), float(line_list[1]), float(
                    line_list[2])
                element_symbol, stereoatom_info = line_list[3], int(
                    line_list[6])
                if stereoatom_info == 3:  # positive charge
                    charge = 1
                if stereoatom_info == 4:  # radical
                    radical = True
                if stereoatom_info == 5:  # negative charge
                    charge = -1
                new_atom = Atom(
                    element_symbol=element_symbol,
                    x=x,
                    y=y,
                    z=z,
                    index=atom_enumerator,
                    charge=charge,
                    radical=radical,
                    total_index=atom_enumerator,
                )
                self.add_atom(new_atom)

            # Bond-block information is extracted (stereobond_info contains an atoms chiral information):
            if len(line_list) == 4 and line_list[0] != "M":
                atom1_index, atom2_index, bond_order_num, stereobond_info = (
                    int(line_list[0]),
                    int(line_list[1]),
                    int(line_list[2]),
                    int(line_list[3]),
                )
                atom1 = self.contained_atoms[atom1_index - 1]
                atom2 = self.contained_atoms[atom2_index - 1]
                bond_order = ["s", "d", "t", "a"][bond_order_num - 1]
                self.add_bond(atom1, atom2, bond_order)
                if stereobond_info == 1:
                    atom1.stereochem = "S"  # CIP rule
                if stereobond_info == 6:
                    atom1.stereochem = "R"

    @classmethod
    def from_smiles(cls, smiles_input, neutralize_molecule=True):
        """ A SMILES-string is used to generate the Molecules InChI-string and it's graph. """
        # molecule instantiation and adding the smiles-string:
        start_time = time.perf_counter()
        mdh_mol = cls()

        # adding InChI-string to the molecule (for finding entries in database) and canonicalization of SMILES:
        molecule = Chem.MolFromSmiles(smiles_input)

        # neutralize atoms
        if neutralize_molecule:
            molecule = neutralize_atoms(molecule)
        mdh_mol.inchi = AllChem.MolToInchi(
            molecule)  # TODO: producing warning, generate inchi's separately?
        mdh_mol.smiles = Chem.MolToSmiles(molecule)

        # generate 3D atomic coordinates based on RDKit's EKTD-method:
        molecule = Chem.AddHs(molecule)
        AllChem.EmbedMolecule(molecule, randomSeed=0xF00D)
        molblock = Chem.MolToMolBlock(molecule)

        # add metadata:
        mdh_mol.coordinate_metadata[0] = rdkit.__name__
        mdh_mol.coordinate_metadata[1] = rdkit.__version__
        mdh_mol.coordinate_metadata[2] = "EKTG"  # TODO: change to ETKDG
        mdh_mol.coordinate_metadata[4] = datetime.datetime.utcnow()

        # 5-tuple: [software, version, method, calc.-time, datetime]
        # self.coordinate_metadata = [None, None, None, None, None]

        # generate/extract cartesian coordinates and add them to the molecule:
        mdh_mol.molblock_data_extraction(molblock, rdkit_molblock=True)

        mdh_mol.coordinate_metadata[3] = round(
            time.perf_counter() - start_time, 3)

        return mdh_mol

    @classmethod
    def from_inchi(cls, inchi_input):
        """ A InChI-string is used to generate the Molecules SMILES-string and it's graph. """
        # molecule instantiation and adding the inchi-string:
        mdh_mol = cls(inchi=inchi_input)

        # adding InChI-string to the molecule (for finding entries in database)
        molecule = Chem.MolFromInchi(inchi_input)
        mdh_mol.smiles = Chem.MolToSmiles(molecule)
        mdh_mol.inchi = Chem.MolToInchi(molecule)

        # generate 3D-coords based on RDKit's EKTD-method:
        molecule = Chem.AddHs(molecule)
        AllChem.EmbedMolecule(molecule, randomSeed=0xF00D)
        molblock = Chem.MolToMolBlock(molecule)

        # generate/extract cartesian coordinates and add them to the molecule:
        mdh_mol.molblock_data_extraction(molblock, rdkit_molblock=True)
        return mdh_mol

    @classmethod
    def from_selfies(cls, selfies_input):  # TODO: fillout
        """ Description. """
        mol = cls()

        # add selfies decoder, add also the @property 'get_selfies'

        # generate/extract cartesian coordinates and add them to the molecule:
        mol.molblock_data_extraction(string_type="selfies",
                                     string=selfies_input)
        return mol

    @classmethod
    def from_molfile(cls, molfile_path):
        """ a molfile is used to generate a molecule. The path of the molfile has to be given. """
        pass

    @classmethod
    def load_molecule_from_mongodb(
        cls,
        InChI,
        atomic_coords='EKTG',  # TODO: change to ETKDG
        db_name="test_database",
        collection_name="molecules",
    ):
        """ a molecule is created from it's MongoDB document.  """
        mdh_mol = cls()
        client = pymongo.MongoClient()
        # select a database (test)
        db = client[f"{db_name}"]
        # select a collection (eq. to sql-table)
        molecules_doc = db[f"{collection_name}"]

        # specified molecule:
        document = molecules_doc.find_one({"InChI": InChI})

        if document is None:
            print(
                f'Search for {InChI!r} returned no matching document from the collection {collection_name!r}.'
            )
        else:
            return Molecule.mongodb_data_extraction(document, atomic_coords)

    @classmethod
    def mongodb_data_extraction(cls,
                                mongodb_dict,
                                atomic_coords='EKTG'):  # TODO: change to ETKDG
        """ Method which extracts data from the MongoDB document (a dictionary) and creates a Molecule-instance.
            By default EKTG-generated atomic coordinates are loaded, but it's also possible to load other
            coordinates. """

        # data are pulled into the Molecule:
        molecule = cls()

        molecule.inchi = mongodb_dict['InChI']
        molecule.smiles = mongodb_dict['canonical_SMILES']

        # extracting atom-data from atoms-key:
        for total_index, element, charge, stereochem, radical in mongodb_dict[
                'atoms']:
            atom = Atom(total_index=total_index,
                        element_symbol=element,
                        charge=charge,
                        stereochem=stereochem,
                        radical=radical)
            molecule.add_atom(atom)

        # bond information are extracted:
        for atom1, atom2, bondorder in mongodb_dict['bonds']:
            atom1 = molecule.contained_atoms[atom1 - 1]
            atom2 = molecule.contained_atoms[atom2 - 1]
            molecule.add_bond(atom1, atom2, bondorder)

        # atomic coordinates are selected and pulled into the molecule:
        for total_index, element, x, y, z in mongodb_dict['xyz_coordinates'][
                atomic_coords]['atomic coordinates']:
            atom = molecule.contained_atoms[total_index - 1]
            if total_index == atom.total_index and element == atom.element:
                atom.x, atom.y, atom.z = x, y, z

        # extract metadata
        metadata = mongodb_dict['xyz_coordinates'][atomic_coords]
        molecule.coordinate_metadata[0] = metadata['software']
        molecule.coordinate_metadata[1] = metadata['version']
        molecule.coordinate_metadata[2] = metadata['method']
        molecule.coordinate_metadata[3] = metadata['calc-time']
        molecule.coordinate_metadata[4] = metadata['datetime']

        return molecule

    def calc_geometry_opt(self, method):
        """perform geometry optimization, using different quantum chemical
        methods. atomic coordinates are saved within db """
        methods_list = ["method1", "method2"]  # use subprocess module
        # 5-tuple: (used software+version, method, datetime, calc-time, machine (OS + CPU))
        self.coordinate_metadata = ()
        pass

    @property
    def atom_coordinates_matrix(self):
        """ returns the atoms coordinates as matrix """
        xyz_lst = []
        for atom in self.contained_atoms:
            xyz_lst.append(
                [atom.total_index, atom.element, atom.x, atom.y, atom.z])
        return xyz_lst

    @property
    def xyz_molblock(self):
        """ generates the molecules XYZ molblock. """
        xyz_molblock = ""
        for atom in self.contained_atoms:
            xyz_molblock += f"{atom.element}    {atom.x}    {atom.y}    {atom.z}\n"
        return xyz_molblock

    @property
    def atom_coordinates_numpy_matrix(self):
        """ returns the atoms coordinates as matrix """
        xyz_lst = []
        for atom in self.contained_atoms:
            xyz_lst.append([atom.x, atom.y, atom.z])
        return np.array(xyz_lst)

    @property
    def atom_properties_list(self):
        """ used for saving atomic properties within the MongoDB. """
        atom_prop_list = []
        for atom in self.contained_atoms:
            atom_prop_list.append([
                atom.total_index,
                atom.element,
                atom.charge,
                atom.stereochem,
                atom.radical,
            ])
        return atom_prop_list

    @property
    def bond_properties_list(self):
        data = []
        molecular_graph = self.mol_graph
        for atom1 in molecular_graph:
            for atom2 in molecular_graph[atom1]:
                # prevent repetitions and sort atoms1 vs atom2:
                if (
                    (atom1.total_index, atom2.total_index),
                        molecular_graph[atom1][atom2],
                ) not in data and (
                    (atom2.total_index, atom1.total_index),
                        molecular_graph[atom1][atom2],
                ) not in data:
                    data.append((
                        atom2.total_index,
                        atom1.total_index,
                        molecular_graph[atom1][atom2],
                    ))
        return data

    def save_in_mongodb(self,
                        db_name="test_database",
                        collection_name="molecules"):
        """By default RDKit-generated coordinates and other molecular features are saved
        within a molecules-collection (MongoDB). If the molecules coordinates are generated using
        quantum chemical geometry optimization (ab initio methods etc.) then only
        the atomic cartesian coordinates are added to the 'atom_coordinates'-dictionary."""

        # connects to default host/port:
        client = pymongo.MongoClient()
        # select a database (test)
        db = client[f"{db_name}"]
        # select a collection (eq. to sql-table)
        molecules_doc = db[f"{collection_name}"]

        # if molecule has yet no document within the Mono database at all:
        if molecules_doc.find_one({"InChI": self.inchi}) is None:
            molecule = {
                "InChI": self.inchi,
                "canonical_SMILES": self.smiles,
                "Molecular_Weight": self.get_molweight,
                "Chemical_Formula": self.chem_formula,
                "Total_Count_Electrons": self.get_total_electrons,
                "bonds": self.bond_properties_list,
                "atoms": self.atom_properties_list,
                "xyz_coordinates": {
                    self.coordinate_metadata[2]: {
                        "software": self.coordinate_metadata[0],
                        "version": self.coordinate_metadata[1],
                        "method": self.coordinate_metadata[2],
                        "calc-time": self.coordinate_metadata[3],
                        "datetime": self.coordinate_metadata[4],
                        "atomic coordinates": self.atom_coordinates_matrix,
                    },
                },
            }
            molecules_doc.insert_one(molecule)
            print(f"{self.inchi} was added to the database.")

        # if atomic coordinates are generated by another method (not RDKit's EKTG-method) add them also to MongoDB:
        elif (molecules_doc.find_one({
                "InChI": self.inchi
        })["xyz_coordinates"].get(self.coordinate_metadata[2]) is None):
            # add atomic coordinates of the calculated method:
            xyz_coordinates = {
                "software": self.coordinate_metadata[0],
                "version": self.coordinate_metadata[1],
                "method": self.coordinate_metadata[2],
                "calc-time": self.coordinate_metadata[3],
                "datetime": self.coordinate_metadata[4],
                "atomic coordinates": self.atom_coordinates_matrix,
            }
            updated_entry = molecules_doc.find_one(
                {"InChI": self.inchi})["xyz_coordinates"][
                    self.coordinate_metadata[2]] = xyz_coordinates
            molecules_doc.update_one(
                {"InChI": self.inchi},
                {"$set": {
                    "xyz_coordinates": updated_entry
                }})

    @staticmethod
    def query_inchi(inchi,
                    db_name="test_database",
                    collection_name="molecules"):
        """ querying for one molecule in the MongoDB using it's InChI """
        # connects to default host/port:
        client = pymongo.MongoClient()
        # select a database (test)
        db = client[f"{db_name}"]
        # select a collection (eq. to sql-table)
        molecules_doc = db[f"{collection_name}"]
        molecule_entry = molecules_doc.find_one({"InChI": inchi})
        return molecule_entry

    # def select_atom_by_total_index(self, total_index):
    #     """ An specific atom-instance from a molecule is selected using the atoms index (int) and
    #     element-symbol (str) by iterating over the molecules atoms.
    #
    #     example: atom = Molecule.select_atom_by_ele_idx(element='C', index='1')
    #              print(atom)           # -> C1
    #              print(atom.__dict__)  # -> displays all stored attributes of the atom-object
    #     """
    #     for num, atom in enumerate(self.mol_graph):
    #         if atom.total_index == total_index:
    #             return atom
    #         elif (num + 1) >= len(self.mol_graph):
    #             raise AtomNotFoundError

    def select_atom_by_ele_idx(self, element, index):
        """An specific atom-instance from a molecule is selected using the atoms index (int) and
        element-symbol (str) by iterating over the molecules atoms.

        example: atom = Molecule.select_atom_by_ele_idx(element='C', index='1')
                 print(atom)           # -> C1
                 print(atom.__dict__)  # -> displays all stored attributes of the atom-object
        """
        for num, atom in enumerate(self.contained_atoms):
            if atom.element == element and atom.index == index:
                return atom
            elif (num + 1) >= len(self.contained_atoms):
                raise AtomNotFoundError(element, index, self.contained_atoms)

    def get_info_of_selected_atom(self, element, index):  # example: 'C3'
        """The specified atom is selected and information about the atom are returned.
        information: radical/charge/stereochem, covalent bonding partners, hybridization,
        cartesian coordinates and their origin
        """

        # select_atom_by_ele_idx-method is used to select the atom:
        selected_atom = self.select_atom_by_ele_idx(element, index)

        # an information-string is concatenated:
        info = ""
        info += f"{selected_atom} has the following properties:\n"

        if selected_atom.radical is True:
            info += f"It is a radical\n"
        if selected_atom.charge > 0:
            info += f"It has a positive charge ({selected_atom.charge})\n"
        if selected_atom.charge < 0:
            info += f"It has a negative charge ({selected_atom.charge})\n"
        if selected_atom.stereochem is not None:
            info += f"It is chiral: {repr(selected_atom.stereochem)}-configured in accordance to CIP rules\n"

        # neighboured atoms -> bonds
        info += f"Covalent bonds: {self.mol_graph[selected_atom]} (s=single, d=double, t=triple bond)\n"

        # (octet rule (8 electrons)) - (elements num valence-electrons) - (number of bonds) = hybridization-type
        # 8 - 4 - 4 = 0 Csp3
        # 8 - 5 - 3 = 0 Nsp3
        # 8 - 5 - 2 = 1 Nsp2
        # hybridization_category: 0 -> sp3, 1 -> sp2, 2 -> sp

        atoms_bonds_lst = self.mol_graph[selected_atom].values()
        atoms_valence_electrons = (
            Atom.valence_electrons[selected_atom.atoms_internal_number] -
            selected_atom.charge)
        hybridization_category = 8 - atoms_valence_electrons - len(
            atoms_bonds_lst)

        if hybridization_category == 0:
            info += f"Hybridization: {selected_atom.element}sp3\n"
        if hybridization_category == 1:
            info += f"Hybridization: {selected_atom.element}sp2\n"
        if hybridization_category == 2:
            info += f"Hybridization: {selected_atom.element}sp\n"

        info += f"Cartesian coordinates: {selected_atom.x, selected_atom.y, selected_atom.z}\n"
        # info += f'origin: {self.coordinate_origin}'
        return info

    @property
    def get_implicit_graph(self):
        """ Returns an implicit graph-representation of the molecule (without hydrogen). """
        # delete all atom with atom.element = 'H'
        implicit_graph = self.mol_graph.copy()

        for atom in self.mol_graph:

            # H-atoms to be removed are saved within this list (while iteration elements of a dict cannot be removed)
            atom_cache = []
            for inner_atom in self.mol_graph[atom]:
                if inner_atom.element == "H":
                    atom_cache.append(inner_atom)

            # now the H-atoms are removed:
            for atom_next in atom_cache:
                del implicit_graph[atom][atom_next]

            # remove all H-atoms from outer dict:
            if atom.element == "H":
                del implicit_graph[atom]

        return implicit_graph

    @property
    def get_molweight(self):
        """ Returns the molecules molecular weight in unit (1 uni == 1g/mol). """
        mol_weight = 0
        for atom in self.contained_atoms:
            mol_weight += Atom.atom_weights[atom.atoms_internal_number]
        return round(mol_weight, 4)

    @property
    def get_total_electrons(self):
        """ Returns the molecules total number of electrons. """
        total_electrons = 0
        for atom in self.contained_atoms:
            total_electrons += Atom.atomic_numbers[atom.atoms_internal_number]
        return total_electrons

    def bond_distances_as_dataframe(self):
        """The distances of each valence bond of a molecule is calculated and returned as dataframe.
        Atom distances are shown in pm.  # TODO: add more than just C1, H1...

        Returns:
               atoms  distance bondorder
        0  (C1, C2)    154.22         s
        1  (C1, H1)    110.49         s
        2  (C1, H2)    111.57         s
        3  (C1, H3)    109.30         s
        4  (C2, C3)    150.76         s

        """
        data = []
        # valence bond information are extracted from the molecules graph
        molecular_graph = self.mol_graph
        for atom1 in molecular_graph:
            for atom2 in molecular_graph[atom1]:
                # calculation of the distance between two points in 3-dim space is based on the Pythagorean theorem
                distance = round(
                    sqrt(((atom1.x - atom2.x)**2 + (atom1.y - atom2.y)**2 +
                          (atom1.z - atom2.z)**2)) * 100,
                    2,
                )
                # prevent repetitions and sort atoms1 vs atom2:
                if (
                    (atom1, atom2),
                        distance,
                        molecular_graph[atom1][atom2],
                ) not in data and (
                    (atom2, atom1),
                        distance,
                        molecular_graph[atom1][atom2],
                ) not in data:
                    data.append(((atom2, atom1), distance,
                                 molecular_graph[atom1][atom2]))
        # binary search algorithm is used to sort the atoms (possible due to Atom's comparison operators)
        result = []
        for e in data:
            insort(result, e)
        df = DataFrame(result, columns=["atoms", "distance", "bondorder"])
        return df

    def visualize_as_3d_plot(
        self,
    ):  # todo: display bond distances as annotation as well: display_bond_dist=true/false
        """ Visualizes the molecule as 3D-plot with Matplotlib. """

        # use dataframe for connected atoms and cut out atom-tuples to generate the bonds:
        df1 = self.bond_distances_as_dataframe()
        only_atoms = df1["atoms"]

        # use the atoms internal list of contained atoms to generate the atoms ina  scatter plot:
        atom_list = self.contained_atoms

        # plotting the graph:
        fig = plt.figure(figsize=plt.figaspect(1) * 1.5,
                         dpi=80)  # dpi for changing the atom balls sizes
        # ax = fig.add_subplot(111, projection='3d')  # what does the 111 mean??
        ax = fig.gca(projection="3d")
        ax.set_axis_off()

        # plot each bond/edge between two atoms:
        for atom1, atom2 in only_atoms:
            ax.plot3D(
                [atom1.x, atom2.x],
                [atom1.y, atom2.y],
                [atom1.z, atom2.z],
                "#525252",
                alpha=0.7,
            )

        # list of all atoms, scatter the atoms and map the element to certain colors
        element_colors = {
            "H": "#d9d9d9",
            "B": "#ffbf80",
            "C": "#000000",
            "N": "#0000ff",
            "O": "#ff3300",
            "F": "#00e600",
            "Si": "#4d4d4d",
            "P": "#ff9900",
            "S": "#e6e600",
            "Cl": "#00e600",
            "Br": "#990000",
            "I": "#6600cc",
        }
        for atom in atom_list:
            color = element_colors[atom.element]
            ax.scatter(
                atom.x,
                atom.y,
                atom.z,
                s=150,
                label="True Position",
                color=color,
                alpha=0.7,
                depthshade=False,
            )

        # add vertices/atoms annotation
        for atom in atom_list:
            coords = (atom.x, atom.y, atom.z)
            annotate3d(
                ax,
                s=str(atom),
                xyz=coords,
                fontsize=10,
                xytext=(-3, 3),
                textcoords="offset points",
                ha="right",
                va="bottom",
            )

        # scale 3d-plot to (solution from sebix last comment on)
        # https://stackoverflow.com/questions/8130823/set-matplotlib-3d-plot-aspect-ratio
        scaling = np.array(
            [getattr(ax, "get_{}lim".format(dim))() for dim in "xyz"])
        # factor = 0.9
        factor = 0.5
        ax.auto_scale_xyz(
            *[[factor * np.min(scaling), factor * np.max(scaling)]] * 3)

        plt.tight_layout()
        plt.show()


