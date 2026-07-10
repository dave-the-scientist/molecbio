"""
Several modules written to handle molecular biology tasks.

-- Currently defines the following modules or packages:
   align -- Functions for dealing with sequence alignments.
   blosum -- the BLOSUM matrices as dictionaries.
   Cluster -- Functions to align and measure PDB files.
   parseblastxml -- Parse BLAST files.
   PDB -- Functions for dealing with protein structure files.
   phylo -- Functions to work with phylogenetic tree files.
   sequ -- Functions for DNA/protein sequences, alignments, etc.

"""
__version__ = '0.7.0'
__all__ = ['align', 'blosum', 'Cluster', 'parseblastxml', 'PDB', 'phylo', 'sequ']

from molecbio import align, blosum, Cluster, parseblastxml, PDB, phylo, sequ

# when I rewrite this, include the visualization functionality from /home/dave/Desktop/work/p_multocida_slp_diversity/strain_prevalence/analyze_hits.py and the blast sequence extraction from /home/dave/Desktop/work/p_multocida_slp_diversity/strain_prevalence/extract_sequences.py
# Move the Cluster functionality into PDB