import os
import sys

from setuptools import find_packages, setup

# needed for setuptools.build_meta to pickup vendored versioneer.py
sys.path.insert(0, os.path.dirname(__file__))
import versioneer  # noqa: E402

# NOTE: This file must remain Python 2 compatible for the foreseeable future,
# to ensure that we error out properly for people with outdated setuptools
# and/or pip.
min_version = (3, 10)
if sys.version_info < min_version:
    error = """
edrixs does not support Python {0}.{1}.
Python {2}.{3} and above is required. Check your Python version like so:

python3 --version

This may be due to an out-of-date pip. Make sure you have pip >= 9.0.1.
Upgrade pip like so:

pip install --upgrade pip
""".format(
        *sys.version_info[:2], *min_version
    )
    sys.exit(error)

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, "README.rst"), encoding="utf-8") as readme_file:
    readme = readme_file.read()

with open(os.path.join(here, "requirements.txt")) as requirements_file:
    # Parse requirements.txt, ignoring any commented-out lines.
    requirements = [
        line
        for line in requirements_file.read().splitlines()
        if not line.startswith("#")
    ]


setup(
    name="edrixs",
    version=versioneer.get_version(),
    description="An open source toolkit for simulating RIXS spectra based on ED",
    long_description=readme,
    author="Brookhaven National Lab",
    author_email="yilinwang@bnl.gov",
    url="https://github.com/NSLS-II/edrixs",
    packages=find_packages(exclude=["docs", "tests", "bin", "examples", "src"]),
    include_package_data=True,
    package_data={
        "edrixs": [
            # When adding files here, remember to update MANIFEST.in as well,
            # or else they will not be included in the distribution on PyPI!
            # 'path/to/data_file',
            "edrixs/*",
            "edrixs/atom_data/*",
        ]
    },
    install_requires=requirements,
    python_requires=">=3.10",
    license="GPL-3.0-or-later",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Development Status :: 2 - Pre-Alpha",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
    ],
)
