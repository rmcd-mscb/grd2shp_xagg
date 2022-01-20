"""Setuptools configuration for Grd2Shp_Xagg."""
from os import path

from setuptools import find_packages
from setuptools import setup

with open("README.md") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

here = path.abspath(path.dirname(__file__))
with open(path.join(here, "requirements.txt")) as requirements_file:
    # Parse requirements.txt, ignoring any commented-out lines.
    requirements = [
        line
        for line in requirements_file.read().splitlines()
        if not line.startswith("#")
    ]

setup_requirements = [
    "pytest-runner",
]

test_requirements = [
    "pytest>=3",
]

setup(
    author="rmcd@usgs.gov",
    author_email="rmcd@usgs.gov",
    python_requires=">=3.8",
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    description="Interpolate gridded input to geometry polygons",
    entry_points={
        "console_scripts": [
            "grd2shp_xagg=grd2shp_xagg.__main__:main",
        ],
    },
    install_requires=requirements,
    license="MIT",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="grd2shp_xagg,grd2shp_xagg",
    name="grd2shp_xagg",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/rmcd-mscb/grd2shp_xagg",
    version="0.0.1-dev0",
    zip_safe=False,
)
