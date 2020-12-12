import pathlib

import setuptools

# from src.feedingorcas import (__author__, __author_email__, __doc_url__,
#                               __src_url__, __version__)

__version__ = "0.1.0"
__author__ = 'Niklas Tiede'
__author_email__ = 'niklastiede2@gmail.com'
__doc_url__ = 'https://feedingorcas.readthedocs.io'
__src_url__ = 'https://github.com/NiklasTiede/feedingORCAs'

setuptools.setup(
    name="feedingorcas",
    version=__version__,
    author=__author__,
    author_email=__author_email__,
    description="A small package for playing around with molecules",
    long_description=pathlib.Path("README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    url=__src_url__,
    project_urls={
        'Documentation': __doc_url__,
        'Source Code': __src_url__,
    },
    license="MIT",
    package_dir={"": "src"},
    # packages=setuptools.find_packages(where="src"),
    packages=setuptools.find_packages("./src"),
    install_requires=[  # production install: python setup.py install or pip install .
        "rdkit",
        "matplotlib",
        "pandas",
        "numpy",
        "pymongo",
    ],
    extras_require={   # pip install -e .[dev]    to install requirements AND dev-packages
        'dev': [
            'pytest',
            # 'pytest-pep8',
            # 'pytest-cov',
        ],
    },
    # extras_require=[
    #         'pytest',
    #         # 'pytest-pep8',
    #         # 'pytest-cov',
    #         # 'yapf',
    #         # 'isort',
    #         # 'sphinx',
    # ],
    platforms="linux",
    python_requires=">=3.5",
    entry_points={"console_scripts": ["feedingorcas = feedingorcas:run_main"]},
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
# entry_points = {"console_scripts": ["feedingorcas = feedingorcas.__main__:run_main"]},
