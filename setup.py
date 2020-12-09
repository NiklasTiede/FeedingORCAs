import pathlib
import re
from os.path import join

import setuptools


def get_version():
    with open(join("src", "feedingorcas", "__init__.py"), "r") as f:
        content = f.read()
    p = re.compile(r'^__version__ = [\'"]([^\'\"]*)[\'"]', re.M)
    return p.search(content).group(1)

from src.feedingorcas import __author__, __author_email__, __url__, __version__

setuptools.setup(
    name="feedingORCAs",
    version=get_version(),
    author=__author__,
    author_email=__author_email__,
    description="A small package for playing around with molecules",
    long_description=pathlib.Path("README.md").read_text(encoding="utf-8"),
    long_description_content_type="text/markdown",
    platforms="linux",
    url="https://github.com/NiklasTiede/feedingORCAs",
    # project_urls={
    #     'Documentation': 'https://feedingorcas.readthedocs.io',
    #     'Source Code': 'https://github.com/NiklasTiede/feedingORCAs',
    # },
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    install_requires=(
        "matplotlib",
        "pandas",
        "numpy",
        "pymongo",
    ),
    license="MIT license",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["feedingorcas = feedingorcas:run_main"]},
)
# entry_points = {"console_scripts": ["feedingorcas = feedingorcas.__main__:run_main"]},
