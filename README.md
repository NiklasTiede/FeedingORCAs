
|     Author:    |            Niklas Tiede             |
|----------------|-------------------------------------|
| Documentation: | https://feedingorcas.readthedocs.io |


feedingORCAs
=============

![image](docs/molecule_animation.gif)


[comment]: <> ([![conda version]&#40;https://img.shields.io/&#41;]&#40;https://anaconda.org/&#41;    # https://shields.io/category/version)

[comment]: <> ([![platform]&#40;https://img.shields.io/&#41;]&#40;&#41;    # https://shields.io/category/platform-support)

![Read the Docs](https://img.shields.io/readthedocs/feedingorcas)

[comment]: <> ([![license]&#40;https://img.shields.io/conda/&#41;]&#40;&#41;    # https://shields.io/category/license)


[comment]: <> ([![Travis CI]&#40;https://img.shields.io/&#41;]&#40;https://travis-ci.com/github/numpy/numpy&#41;     # https://shields.io/category/build)

[comment]: <> ([![codecov]&#40;https://img.shields.io/&#41;]&#40;https://codecov.io/&#41;    # https://shields.io/category/coverage)

[comment]: <> ([![codacy]&#40;https://img.shields.io/&#41;]&#40;&#41;    # https://shields.io/category/analysis)

[comment]: <> (![total lines]&#40;https://img.shields.io/&#41;    # https://shields.io/category/size)

[comment]: <> (![repo size]&#40;https://img.shields.io/&#41;    # https://shields.io/category/size)



goal: an interface for ORCA and a db, rdkit adds also functionality to 
let you work with chemical data.

The goal of this small project is to populate a relational database with
molecular data using RDKit and the quantum chemistry software ORCA.
Molecules can be created using different chemical file formats and
viewed as 3D plot. Several quantum chemical methods can be used for the
calculation of molecule coordinates and other properties.

Requirements
------------

You need Anaconda, a linux environment and the quantum chemistry
software ORCA which is accessible for academic purposes for free. But
even though if you have no access to ORCA, the project lets you build
your own database for molecules by retrieving and calculating molecule
data.

Installation
------------

The dependencies are installed using the `environment_linux.yml` file.
After downloading the code, you can go into the projects folder and
create & activate the conda environment:

``` {.sourceCode .bash}
1. conda env create --file environment_linux.yml
2. conda activate molecule-data-handler
```

Currently the quantum chemistry software ORCA is used to calculate
quantum chemical properties of the data sets. The software can be
requested for free at <https://orcaforum.kofo.mpg.de/app.php/portal> if
used for academic purposes only.

The jupyter notebook `Tutorial.ipynb` gives you an introduction to some
basic functionality of the project.
