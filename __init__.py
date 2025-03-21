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
__version__ = '0.5'
__author__ = 'Dave Curran (daves.binf.tools@gmail.com)'
__all__ = ['PDB', 'Cluster', 'sequ', 'blosum', 'align', 'rosetta', 'aligners', 'phylo']

# when I rewrite this, include the alignment functionality from /home/dave/Desktop/work/p_multocida_slp_diversity/strain_prevalence/calculate_alignment_identity.py and the visualization functionality from /home/dave/Desktop/work/p_multocida_slp_diversity/strain_prevalence/analyze_hits.py and the blast sequence extraction from /home/dave/Desktop/work/p_multocida_slp_diversity/strain_prevalence/extract_sequences.py

try:
    from molecbio import Cluster
except ImportError:
    print("Could not import the Cluster module, so its functionality has been disabled. This is usually because the 'numpy' module could not be found.")
    Cluster = None
from molecbio import PDB, sequ, blosum, align, aligners
