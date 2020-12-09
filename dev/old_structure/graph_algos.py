
from src.molecule_classes_and_functions import Atom, Molecule
from pprint import pprint

graph = {
    'A': ['B', 'C'],
    'B': ['D', 'E'],
    'C': ['F'],
    'D': [],
    'E': ['F'],
    'F': [],
}

# pprint(graph, width=20)

# visited = []
# queue = []
#
#
# def breadth_first_search(visited, graph, node):
#     visited.append(node)
#     queue.append(node)
#
#     while queue:
#         print(f'visited: {visited}')
#         print('yet unvisited:', queue)
#         s = queue.pop(0)
#         print(f'\ntraversal on element {s!r}:')
#
#         for neighbor in graph[s]:
#             print(f'{neighbor } is neighbor of {s}')
#
#             if neighbor not in visited:
#                 print(f'{neighbor} was not yet visited!')
#                 visited.append(neighbor)
#                 queue.append(neighbor)


# breadth_first_search(visited, graph, 'A')

# # cool feature:
# ma_lst = [1, 2, 3]
# while ma_lst:
#     print(f'element {str(ma_lst.pop(0))!r} of the queue is removed')


bidir_graph = {
    'A': ['B', 'C'],
    'B': ['A', 'D', 'E'],
    'C': ['A', 'F'],
    'D': ['B'],
    'E': ['B', 'F'],
    'F': ['C', 'E'],
}

# breadth_first_search(visited, bidir_graph, 'A')

bidict_graph = {
    'A': {'B': 1, 'C': 1},
    'B': {'A': 1, 'D': 1, 'E': 1},
    'C': {'A': 1, 'F': 1},
    'D': {'B': 1},
    'E': {'B': 1, 'F': 1},
    'F': {'C': 1, 'E': 1},
}

# breadth_first_search(visited, bidict_graph, 'A')

# -----------------

# visited = []
# queue = []
#
#
# def breadth_first_search(visited, graph, node):
#     visited.append(node)
#     queue.append(node)
#
#     while queue:
#         print(f'visited: {visited}')
#         print('queue:', queue)
#
#         # traversal on an element from queue:
#         s = queue.pop(0)
#         print(f"\ntraversal on '{s!r}':")
#
#         for neighbor in graph[s]:
#             print(f"'{neighbor}' is neighbor of '{s}'")
#
#             if neighbor not in visited:
#                 print(f"'{neighbor}' was not yet visited! Though appended to visited/queue.")
#                 visited.append(neighbor)
#                 queue.append(neighbor)

# ethane = Molecule.from_smiles('CCOC#C')
# breadth_first_search(visited, ethane.mol_graph, Atom('H', index=1))

# visited atoms are a list of all connected molecules -> generate lists of stereochem/charge etc
# cleave a bond and apply BFS on fragment, copy graph, reindex (update properties etc)
# combine two molecule fragments to form a new one (diverse vs focused library)

# ----------------------

# methane = Molecule.from_smiles('CC')
# test_molgraph = methane.mol_graph
# pprint(test_molgraph)
# print()
#
# breadth_first_search(visited, test_molgraph, Atom('H', index=1))


# ----------
# BFS creating a new moelcule using a fragment of the molecule


# def breadth_first_search(graph, node):
#     """ using BFS to copy the whole molecule """
#     visited = []
#     queue = []
#     visited.append(node)
#     queue.append(node)
#
#     new_mol = Molecule()
#
#     while queue:
#         print(f'visited: {visited}')
#         print('queue:', queue)
#
#         # traversal on an element from queue:
#         s = queue.pop(0)
#         print(f"\ntraversal on '{s!r}':")
#
#         # add traversed atom to the outer dict:
#         if s not in new_mol.mol_graph:
#             print(f'{s} was added to the new molecules outer dict "mol_graph"')
#             new_mol.mol_graph[s] = dict()
#
#         for neighbor in graph[s]:
#             print(f"'{neighbor}' is neighbor of '{s}' and their bondorder is: {graph[s][neighbor]}")
#
#             # add traversed atom to inner dict (with it's bond order):
#             if neighbor not in new_mol.mol_graph[s]:
#                 new_mol.mol_graph[s][neighbor] = graph[s][neighbor]
#
#             if neighbor not in visited:
#                 print(f"'{neighbor}' was not yet visited! Though appended to visited/queue.")
#                 visited.append(neighbor)
#                 queue.append(neighbor)
#
#     return new_mol
#
# methane = Molecule.from_smiles('CC')
# test_molgraph = methane.mol_graph
# pprint(test_molgraph)
# print()
#
# mol = breadth_first_search(test_molgraph, Atom('C', index=2))  # selected_atom
# print(f'\n\n\n')
# print(mol.contained_atoms)
# pprint(mol.mol_graph, width=80)
# print(mol)

# -----
# use BFS to copy only a part of the molecule
import copy


def reindex(molecule):
    """ reindexes a molecule (total index and index) """
    total_index, temp_saved_atoms = 0, []
    molecule.element_counter_index = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    for atom in molecule.contained_atoms:
        total_index += 1
        atom.total_index = total_index
        temp_saved_atoms.append(atom)
    for atom in temp_saved_atoms:
        molecule.add_atom(atom)
    return molecule


def BFS_cleave_single_bond(molecule, atom1, atom2):
    """ using BFS to cleave a bond and create 2 fragments. """

    molecule = copy.deepcopy(molecule)

    # accessing the selected atoms:
    atom1 = molecule.contained_atoms[molecule.contained_atoms.index(atom1)]
    atom2 = molecule.contained_atoms[molecule.contained_atoms.index(atom2)]

    graph = molecule.mol_graph

    # specified bond is cleaved:
    del graph[atom1][atom2]
    del graph[atom2][atom1]

    pprint(graph)

    fragments = []
    for atom in [atom1, atom2]:
        queue, visited = [], []
        # visited.append(atom)
        visited.append(atom1)
        visited.append(atom2)

        queue.append(atom)
        new_mol = Molecule()
        while queue:
            s = queue.pop(0)
            if s not in new_mol.contained_atoms:
                new_mol.contained_atoms.append(s)
            if s in molecule.contained_atoms:
                molecule.contained_atoms.pop(molecule.contained_atoms.index(s))
            if s not in new_mol.mol_graph:
                new_mol.mol_graph[s] = dict()
            for neighbor in graph[s]:
                if neighbor not in new_mol.mol_graph[s]:
                    new_mol.mol_graph[s][neighbor] = graph[s][neighbor]
                if neighbor not in visited:
                    visited.append(neighbor)
                    queue.append(neighbor)
        new_mol = reindex(new_mol)
        fragments.append(new_mol)
    return fragments


ethane = Molecule.from_smiles('CCOC#C')

mol1, mol2 = BFS_cleave_single_bond(ethane, Atom('C', index=1), Atom('C', index=2))  # selected_atom
print('mol1:')
print(mol1.contained_atoms)
print([(atom.x, atom.y, atom.z) for atom in mol1.contained_atoms])
pprint(mol1.mol_graph, width=60)

print()
print('mol2:')
print(mol2.contained_atoms)
pprint([(atom.x, atom.y, atom.z, atom.stereochem, atom.radical, atom.element, atom.charge) for atom in mol2.contained_atoms])
# mol2.contained_atoms[0].index = 9
# mol2.contained_atoms[1].index = 10
# mol2.contained_atoms[4].index = 11
pprint(mol2.mol_graph, width=60)

# # still keyerror
# for atom in mol2.mol_graph:
#     print(atom)
#     for atom2 in mol2.mol_graph[atom]:
#         print(atom2)
#


import copy


def BFS_reindexing(graph):
    """ using BFS to reindex the molecule"""

    element_counter_index = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    total_index = 0
    atoms_lst = []

    # pick random atom:
    import random
    entry_list = list(graph)
    rnd_atom = random.choice(entry_list)

    visited = []
    queue = []
    visited.append(rnd_atom)
    queue.append(rnd_atom)
    new_mol = Molecule()

    mapping_atoms = {}
    for atom in graph:
        new = copy.deepcopy(atom)
        new_mol.add_atom(new)
        new_atom = new_mol.contained_atoms[-1]
        mapping_atoms[atom] = new_atom
    print('mapped atoms:', mapping_atoms)

    while queue:
        print(f'visited: {visited}')
        print('queue:', queue)
        s = queue.pop(0)
        print(f"\ntraversal on '{s!r}' | ele_index: {s.index} | tot_index: {s.total_index} | ele: {s.atoms_internal_number}")
        print(f'old atom: {s}, new, mapped atom: {mapping_atoms[s]}')

        if s not in new_mol.mol_graph:
            new_mol.mol_graph[mapping_atoms[s]] = dict()

        for neighbor in graph[s]:
            print(f"'{neighbor}' is neighbor of '{s}' and their bondorder is: {graph[s][neighbor]}")

            new_mol.mol_graph[mapping_atoms[s]][mapping_atoms[neighbor]] = graph[s][neighbor]

            if neighbor not in visited:
                print(f"'{neighbor}' was not yet visited! Though appended to visited/queue.")
                visited.append(neighbor)
                queue.append(neighbor)
    return new_mol



import random

def select_random_atom(self):
    """ selects an atom of the molecule randomly. """
    return random.choice(self.contained_atoms)

# delete coordinates and import new ones?!


from math import sqrt
from pandas import DataFrame
from bisect import insort


def bond_distances_as_dataframe(molecule):
    """ The distances of each valence bond of a molecule is calculated and returned as dataframe.
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
    molecular_graph = molecule.mol_graph
    pprint(molecular_graph)
    for atom1 in molecular_graph:
        for atom2 in molecular_graph[atom1]:
            # calculation of the distance between two points in 3-dim space is based on the Pythagorean theorem
            distance = round(sqrt(((atom1.x - atom2.x) ** 2 + (atom1.y - atom2.y) ** 2 + (atom1.z - atom2.z)
                                   ** 2)) * 100, 2)
            # prevent repetitions and sort atoms1 vs atom2:
            if ((atom1, atom2), distance, molecular_graph[atom1][atom2]) not in data and \
                    ((atom2, atom1), distance, molecular_graph[atom1][atom2]) not in data:
                data.append(((atom2, atom1), distance, molecular_graph[atom1][atom2]))
    # binary search algorithm is used to sort the atoms (possible due to Atom's comparison operators)
    result = []
    for e in data:
        insort(result, e)
    df = DataFrame(result, columns=['atoms', 'distance', 'bondorder'])
    return df


# def BFS_atom_distances(molecule, atom):
#     """ using BFS to cleave a bond and create 2 fragments. """
#     atom = molecule.contained_atoms[molecule.contained_atoms.index(atom)]
#     graph = molecule.mol_graph
#     queue, visited = [], []
#     visited.append(atom)
#     queue.append(atom)
#     while queue:
#         s = queue.pop(0)
#         print(s)
#         for neighbor in graph[s]:
#             if neighbor not in visited:
#                 visited.append(neighbor)
#                 queue.append(neighbor)
#

