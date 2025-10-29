#!/usr/bin/env python3
"""
Setup script for the EC Proteomics Analysis package.
"""

from setuptools import setup, find_packages
import pathlib

# Get the long description from the README file
HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

setup(
    name="ec-proteomics",
    version="1.0.0",
    description="A comprehensive data science pipeline for analyzing endothelial cell proteomics data",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/sdurot/EC_proteomics",
    author="sdurot",
    author_email="",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "scikit-learn>=1.0.0",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "plotly>=5.0.0",
        "statsmodels>=0.12.0",
        "multipy>=0.15",
        "openpyxl>=3.0.0",
        "xlrd>=2.0.0",
        "tqdm>=4.60.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "black>=21.0.0",
            "flake8>=3.9.0",
            "jupyter>=1.0.0",
            "ipykernel>=6.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "ec-proteomics=proteomics_analysis:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)