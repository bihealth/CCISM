#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import os.path

from setuptools import setup, find_packages


def parse_requirements(path):
    """Parse ``requirements.txt`` at ``path``."""
    requirements = []
    with open(path, "rt") as reqs_f:
        for line in reqs_f:
            line = line.strip()
            if line.startswith("-r"):
                fname = line.split()[1]
                inner_path = os.path.join(os.path.dirname(path), fname)
                requirements += parse_requirements(inner_path)
            elif line != "" and not line.startswith("#"):
                requirements.append(line)
    return requirements


with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

test_requirements = parse_requirements("requirements/test.txt")
install_requirements = parse_requirements("requirements/base.txt")

setup(
    author="Benedikt Obermayer",
    author_email="benedikt.obermayer@bih-charite.de",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    entry_points={"console_scripts": ("scitcem = scitcem.cli:main",)},
    description="Single-cell Identification of Tumor Cells by EM",
    install_requires=install_requirements,
    license="MIT license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="single-cell, somatic variants, bioinformatics",
    name="scitcem",
    packages=find_packages(include=["scitcem"]),
    test_suite="tests",
    tests_require=test_requirements,
    url="https://github.com/bihealth/scitcem",
    zip_safe=False,
)
