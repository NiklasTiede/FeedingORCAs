"""
**********
Exceptions
**********

Base exceptions and errors for MoleculeDataHandler.
"""

# 'from exceptions import *'-statement is affected by __all__ -> _fct (with an underscore) are normally not imported
# but the __all_ list of strings enables this. These 2 comment lines can be deleted in the end
__all__ = [
    "AtomNotFoundError",
    "__MoleculeDataHandlerException",
]


class __MoleculeDataHandlerException(Exception):
    """Base class for exceptions in MoleculeDataHandler."""


class AtomNotFoundError(__MoleculeDataHandlerException):
    """Throws an Exception if no valid atom of a molecule is selected.

    Example: methane (CH4) contains the following atoms: [C1, H1, H2, H3, H4]
             selecting the atom 'C2' will throw the AtomNotFoundError
    """
    def __init__(self, element_symbol, index, contained_atoms):
        self.element = element_symbol
        self.index = index
        self.contained_atoms = contained_atoms

    def __str__(self):
        return str(
            f"The selected atom {repr(self.element+str(self.index))} could not be found "
            f"in the molecules list of contained atoms: \n{self.contained_atoms}"
        )
