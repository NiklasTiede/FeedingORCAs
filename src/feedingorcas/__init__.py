"""
feedingORCAs - ORCA interface
=============================
A small interface to the QM software ORCA for calculating molecular
properties, storing them in a MongoDB and visualization.

"""
# Copyright (c) 2020, Niklas Tiede.
# All rights reserved. Distributed under the MIT License.

from .core import Molecule, MoleculeList, MoleculeMongoDB
from .__main__ import __version__, __author__, __author_email__, __doc_url__, __src_url__

