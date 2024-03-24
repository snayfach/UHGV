from setuptools import find_packages, setup

setup(
    name="uhgv",
    version="0.0.1",
    packages=find_packages(),
    license="Modified BSD",
    description="Toolkit associated with the Unified Human Gut Virome Catalog",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    install_requires=[
        "biopython",
        "importlib-metadata>=0.12; python_version<'3.8'",
        "psutil"
    ],
    python_requires=">=3.6",
    entry_points={"console_scripts": ["uhgv-tools=uhgv.cli:cli"]},
    url="https://github.com/snayfach/UHGV-tools",
    keywords=["bioinformatics", "genomics", "metagenomics", "viromics"],
    author="Stephen Nayfach",
    classifiers=[
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Topic :: Software Development :: Libraries",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
    ],
)
