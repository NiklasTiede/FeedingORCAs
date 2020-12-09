
|     Author:    |            Niklas Tiede             |
|----------------|-------------------------------------|
| Documentation: | https://feedingorcas.readthedocs.io |


feedingORCAs
=============

![image](docs/molecule_animation.gif)

<img src="https://github.com/favicon.ico" width="48">

![CI](https://github.com/fastai/fastai/workflows/CI/badge.svg)
[![PyPI](https://img.shields.io/pypi/v/fastai?color=blue&label=pypi%20version)](https://pypi.org/project/fastai/#description) 
[![Conda (channel only)](https://img.shields.io/conda/vn/fastai/fastai?color=seagreen&label=conda%20version)](https://anaconda.org/fastai/fastai)
[![Build fastai images](https://github.com/fastai/docker-containers/workflows/Build%20fastai%20images/badge.svg)](https://github.com/fastai/docker-containers) 
![docs](https://github.com/fastai/fastai/workflows/docs/badge.svg)

![Build fastai images](https://github.com/fastai/docker-containers/workflows/Build%20fastai%20images/badge.svg) 
![Build nbdev images](https://github.com/fastai/docker-containers/workflows/Build%20nbdev%20images/badge.svg) 
![Build nbdev docs](https://github.com/fastai/docker-containers/workflows/Build%20nbdev%20docs/badge.svg) 
![Build CI Containers](https://github.com/fastai/docker-containers/workflows/Build%20CI%20Containers/badge.svg)

[![Travis CI](https://img.shields.io/travis/com/numpy/numpy/master?label=Travis%20CI)](
    https://travis-ci.com/github/numpy/numpy)
[![Azure](https://dev.azure.com/numpy/numpy/_apis/build/status/azure-pipeline%20numpy.numpy)](
    https://dev.azure.com/numpy/numpy/_build/latest?definitionId=5)
[![codecov](https://codecov.io/gh/numpy/numpy/branch/master/graph/badge.svg)](
    https://codecov.io/gh/numpy/numpy)

[![Python](https://img.shields.io/pypi/pyversions/tensorflow.svg?style=plastic)](https://badge.fury.io/py/tensorflow)
[![PyPI](https://badge.fury.io/py/tensorflow.svg)](https://badge.fury.io/py/tensorflow)





[![conda version]()]()    # https://shields.io/category/version
[![platform]()]()    # https://shields.io/category/platform-support
[![docs]()]()    # https://shields.io/category/build
[![license](https://img.shields.io/conda/l/pandas/pandas)]()    # https://shields.io/category/license


[![Travis CI]()]()     # https://shields.io/category/build
[![codecov]()]()    # https://shields.io/category/coverage
[![codacy]()]()    # https://shields.io/category/analysis

![total lines]()    # https://shields.io/category/size
![repo size]()    # https://shields.io/category/size


[comment]: <> (test )


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
