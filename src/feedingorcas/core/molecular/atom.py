
from typing import Type


class Atom(object):
    """The purpose of the Atom-class' is to improve the data encapsulation of a molecule. Properties which are
    atom-related (element type and its derived properties like atomic number, weights, valence electrons etc.) are
    stored within the atom objects. The instantiated atoms are organized within a molecule object (See Molecule-class).
    The Atom-class lets you generate 12 different atom types (see atom_symbols). The atom-types are connected via their
    index to their respective derived properties.
    """

    element_symbols = ("H", "B", "C", "N", "O", "F", "Si", "P", "S", "Cl",
                       "Br", "I")
    atomic_numbers = (1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53)
    atom_weights = (
        1.008,
        10.811,
        12.0096,
        14.007,
        15.999,
        18.998,
        28.086,
        30.0738,
        32.06,
        35.45,
        79.904,
        126.904,
    )
    valence_electrons = (1, 3, 4, 5, 6, 7, 4, 5, 6, 7, 7, 7)
    sorted_element_symbols = (
        "C",
        "Si",
        "B",
        "N",
        "P",
        "O",
        "S",
        "I",
        "Br",
        "Cl",
        "F",
        "H",
    )

    def __init__(
        self,
        element_symbol: str,
        x=None,
        y=None,
        z=None,
        stereochem=None,
        charge=0,
        radical=False,
        index=None,
        total_index=None,
    ):
        """To initialize an atom only defining its element-type is mandatory. The element-types internal number (from
        0 to 11) is then generated which is used to access the derived atom-properties (stored as class attributes).
        Optional Arguments are the atoms cartesian coordinates (x, y, z), its stereo information (stereochem) and
        charge (charge).

        :argument

            - element_symbol ('str'):  Element of the atom                              (example: 'C')
            -            x ('float'):  x-component of coordinates, unit is Angström     (example: -0.873)
            -            y ('float'):  y-component                                      (see above)
            -            z ('float'):  z-component                                      (see above)
            -     stereochem ('str'):  chiral atoms, R- or S-configuration (CIP-Rule)   (example: 'R')
            -         charge ('int'):  N has 5 valence electrons, if charge is -1 -> 6  (example: -1)
            -          index ('int'):  atoms index (CH4: C1, H1, H2, H3, H4)            (example: 4)

        """

        if element_symbol not in Atom.element_symbols:
            raise ValueError(
                "Atom could not be initialized. You must choose "
                "between 'H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I'."
            )
        self.element = element_symbol
        self.atoms_internal_number = Atom.element_symbols.index(element_symbol)
        self.x = x
        self.y = y
        self.z = z
        self.stereochem = stereochem
        self.charge = charge
        self.radical = radical
        self.index = index
        self.total_index = total_index

    def __eq__(self, other):
        """The equal operator is defined for comparing if two atoms are equal/identical. Each instantiated
        atom can be discriminated by it's index and element type."""
        if (isinstance(other, Atom) and self.element == other.element
                and self.index == other.index):
            return True
        else:
            return False

    def __ne__(self, other):
        """ The not-equal operator is defined for comparing if two atoms are not equal/identical (see eq-operator). """
        if isinstance(other, Atom) and (self.element != other.element
                                        or self.index != other.index):
            return True
        else:
            return False

    def __lt__(self, other):
        """The lower-than operator is defined for sorting the atoms within the Molecule-class. Atoms are sorted
        based on their index (C1, C2, ... H1, H2, ...) and their bond valency (C, ... N, ... O, ... H)."""
        if isinstance(other, Atom) and self.element != other.element:
            return self.sorted_element_symbols.index(
                self.element) < self.sorted_element_symbols.index(
                    other.element)
        if isinstance(other, Atom) and self.element == other.element:
            return self.index < other.index

    def __gt__(self, other):
        """ The lower-than operator is defined for sorting the atoms within the Molecule-class (see lt-operator). """
        if isinstance(other, Atom) and self.element != other.element:
            return self.sorted_element_symbols.index(
                self.element) > self.sorted_element_symbols.index(
                    other.element)
        if isinstance(other, Atom) and self.element == other.element:
            return self.index > other.index

    def __str__(self):
        """Returns a string-representation of a molecule's atom to let you discriminate the atoms within a molecule
        from each other. Example: 'C3' (carbon-atom with index 3)"""
        return Atom.element_symbols[self.atoms_internal_number] + str(
            self.index)

    def __repr__(self):
        """ Gives the atom a representation within collections. example: 'C3' (carbon-atom with index 3) """
        return Atom.element_symbols[self.atoms_internal_number] + str(
            self.index)

    def __hash__(self):
        """Objects must be hashable for usage in a dict. Each atom's hash is generated from its
        index (int) and element (str)."""
        return hash(f"{self.element}{self.index}")

