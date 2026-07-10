molecbio
========

This Python package was written to handle various bioinformatics tasks. It defines the following modules:
- align - Quantify and manipulate DNA/RNA and protein sequence alignments.
- blosum - Data objects representing the BLOSUM30, 45, 50, 62, and 80 matrices.
- Cluster - Functions to align protein structures in 3D space, and measure the differences between them.
- parseblastxml - Parse the results of a BLAST results XML file.
- PDB - Functions to parse and manipulate protein structures in PDB files.
- phylo - Parse, quantify, and manipulate phylogenetic tree files in several formats.
- sequ - Functions and classes to work with DNA/RNA or protein sequence files.


INSTALLATION
------------

Though this package is not available on PyPI, the best installation method is to use the tool pip:

    pip install git+https://github.com/dave-the-scientist/molecbio.git

If you have an earlier version installed, you first need to uninstall it with:

    pip uninstall -y molecbio

If you have a very early version, you may need to uninstall it by manually deleting the molecbio folder. You can find where it is installed by importing molecbio in a Python session, and examining the variable `molecbio.__file__`.
