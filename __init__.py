"""
Several modules written to handle molecular biology tasks.

-- Currently defines the following modules or packages:
   align -- Functions for dealing with sequence alignments.
   blosum -- the BLOSUM matrices as dictionaries.
   pdb -- Functions for dealing with protein structure files.
   sequ -- Functions for DNA/protein sequences, alignments, etc.
   rosetta -- A package to run and parse various rosetta programs.
     Note that this package is not implicitly loaded with molecbio.
"""
__version__ = '0.3'
__author__ = 'Dave Curran (curran.dave.m@gmail.com)'
__all__ = ['PDB', 'Cluster', 'sequ', 'blosum', 'align', 'rosetta', 'aligners']

try: import Cluster
except ImportError:
    print "Could not import the Cluster module, so its functionality has been disabled. This is usually because the 'numpy' module could not be found."
    Cluster = None
import PDB, sequ, blosum, align, aligners
