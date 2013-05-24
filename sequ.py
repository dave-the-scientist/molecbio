#!/usr/bin/python
# Author: Dave Curran
# Date:   April 2012
# # # # # # # # # # #
"""
Various functions and variables to deal with DNA or protein sequences.

Class:
    -- Sequence(self, name='Unnamed sequence', sequence='',
                allowedChars='-_*?', onlyUpper=True) -- An objected created
       to manipulate and store sequence data. Several optional arguments can
       be given on instantiation.
       
       Attributes:
       -- name -- The name as a string.
       -- sequence / seq -- The sequence as a string.
       
       Methods:
       -- translate() -- Returns DNA-amino acid translation as a string.
       -- invcomplement() -- Returns inverse complement of the DNA sequence.
       -- fasta() -- Returns the name and sequence in a fasta-formatted string.
       -- append(sequence) -- Appends the sequence from another Sequence object
          to its own sequence.
       -- isNucleotide() -- Tests if its sequence contains only A, C, T, G or U,
          otherwise returns False. This means is likely a protein sequence.

Variables:
    -- complement -- Dictionary to find the complement to some DNA/RNA base.
    -- full_complement -- Dictionary to find the complement to some DNA/RNA base.
       Also handles non-specific base codes like R, Y, W, S, etc.
    -- codontable -- Dictionary to translate codons to amino acids.
    -- peptide3to1 -- Dictionary to translate 3-letter residue codes to 1-letter.

Functions:
    Fasta functions:
    -- fasta(sequence, line=60, spaces=True, numbers=True) -- Converts some
       sequence into a fasta format, returns as a string.
    -- parsefasta(textobj) -- Takes a file object or list of text lines and
       returns a list of Sequence objects.
    -- loadfasta(filepath) -- Takes a file path string, returns list of Sequences.
    -- savefasta(seqList, filepath) -- Saves list of Sequences to filepath.
    -- cleanfasta(filepath) -- Overwrites sequences at filepath, formatting them.
    
    Data functions:
    -- translate(sequence, unknownChar='?', stopcodonChar='_') -- Takes a string
       or list of strings, returns amino acid translation.
    -- invcomplement(sequence) -- Returns inverse complement of the given string.
    -- findsub(seqList, sub) -- Takes a sequence string 'sub', and searches the
       Sequence objects in seqList for a match. Returns a list of names.
    -- translatepeptide3to1(sequence, unknownChar='?') -- Translates a list of
       3 letter residue codes to 1 letter, returns as a string.
    -- findORFs(sequence, minLength=60, negStrand=False) -- Takes a string of
       DNA sequence, looks for open reading frames greater than 'minLength', and
       returns a list of Sequence objects named for their location. If negStrand
       is True, searches the inverse complement instead.
    -- calcIdentity(sequence1, sequence2) -- Does not perform an alignment, but
       reads through both sequences, counting each match. If both sequences have
       a gap at the same place, it is ignored, not counting for a match or for
       the total count. Takes two sequence objects, or any two iterables, returns
       a float and a string, the percentage identity and the matches / total.
"""
from __future__ import with_statement # Needed for python 2.5
import itertools

# # # # # # # # # #  Variables  # # # # # # # # # #
"""Dictionary of complementary bases for DNA or RNA."""
complement = {'A':'T','T':'A','C':'G','G':'C','U':'A'}
full_complement = {
'A':'T','T':'A','C':'G','G':'C','U':'A',
'R':'Y','Y':'R','W':'S','S':'W','M':'K','K':'M',
'B':'V','V':'B','H':'D','D':'H' }
codontable = {  
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',  
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',  
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',  
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',  
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',  
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',  
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',  
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',  
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',  
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',  
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',  
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',  
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',  
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',  
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',  
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',  
    }
"""Dictionary to translate DNA into amino acids."""
peptide3to1 = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
    'CYS':'C', 'GLU':'E', 'GLN':'Q', 'GLY':'G',
    'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
    'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
    'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'
    }
"""Dictionary to change 3-letter amino acid codes to 1-letter."""

# # # # # # # # # #  I/O Functions  # # # # # # # # # #
def fasta(sequence, line=60, spaces=True, numbers=True):
    """FASTA-formats some sequence string."""
    if not spaces: seq = '\n'.join(__chunksequence(sequence, line))
    else:
        l = []
        if not numbers:
            for i, s in enumerate(__chunksequence(sequence, line)):
                l.append('          %s' % (' '.join(__chunksequence(s,10)) ))
        else:
            for i, s in enumerate(__chunksequence(sequence, line)):
                l.append('%9d %s' % (i*line+1, ' '.join(__chunksequence(s,10)) ))
        seq = '\n'.join(l)
    return seq

def parsefasta(textobj, onlyThese=None):
    """Returns a list of Sequence objects from the lines or file object.
    If onlyThese is a list of strings, only those sequences that start
    with a string in that list will be collected."""
    if onlyThese:
        onlyThese = tuple(onlyThese)
    seqs, buff, curName, curDescript = [], [], None, ''
    for line in textobj:
        if line.startswith('>'):
            if buff:
                seqs.append(Sequence(name=curName, description=curDescript,
                                     sequence=''.join(buff)))
                buff = []
            if not onlyThese or line[1:].startswith(onlyThese):
                curName, _, curDescript = line[1:].strip().partition(' ')
            else:
                curName, curDescript = None, ''
        else:
            if not curName: continue
            line = line.strip()
            if line: buff.append(line)
    if buff and curName:
        seqs.append(Sequence(name=curName, description=curDescript,
                             sequence=''.join(buff)))
    return seqs
def loadfasta(filepath, onlyThese=None):
    """Returns a list of Sequence objects from the filepath."""
    with open(filepath, 'rb') as f:
        seqs = parsefasta(f, onlyThese)
    return seqs
def savefasta(seqList, filepath, line=60, spaces=True, numbers=True):
    """Saves the given list of Sequences to filepath in fasta format."""
    seqStr = '\n'.join(seq.fasta(line, spaces, numbers) for seq in seqList)
    with open(filepath, 'wb') as f:
        f.write(seqStr)
def cleanfasta(filepath):
    seqs = loadfasta(filepath)
    buff = [seq.__str__() for seq in seqs]
    #buff = map(str, seqs)
    f = open(filepath, 'wb')
    f.write('\n'.join(buff)); f.close()

# # # # # # # # # #  Basic Functions  # # # # # # # # # #
def invcomplement(sequence):
    """Returns the inverse complement of some DNA sequence."""
    return ''.join([complement.get(c,'N') for c in sequence[:].upper()][::-1])

def translate(sequence, unknownChar='?', stopcodonChar='_'):
    """Translates the sequence from DNA to amino acids."""
    seq = ''.join(codontable.get(codon.upper(), unknownChar) for codon
                  in __chunksequence(sequence,3,True))
    if stopcodonChar != '_': return seq.replace('_', stopcodonChar)
    else: return seq

def translatepeptide3to1(sequence, unknownChar='?'):
    """Translates three-letter residue codes to the one-letter code."""
    return ''.join(peptide3to1.get(res.upper(), unknownChar) for res in sequence)

def findsub(seqList, sub):
    """Checks each Sequence or string in seqList for the string sub."""
    sub = sub.upper()
    matches = []
    for seq in seqList:
        if sub in seq: matches.append(seq.header)
    return matches

def findORFs(sequence, minLength=60, negStrand=False):
    def tripScanner(seq):
        subs = itertools.tee(seq, 3)
        for i, sub in enumerate(subs): next(itertools.islice(sub, i, i), None)
        return itertools.izip(*subs)
    stopCodons = ('TGA', 'TAA', 'TAG')
    orfs = []
    sequence = sequence.upper()
    if negStrand: sequence = invcomplement(sequence)
    seqLen = len(sequence)
    for i, trip in enumerate(tripScanner(sequence)):
        if trip == ('A', 'T', 'G'):
            seq = ''.join(itertools.takewhile(lambda codon: codon not in stopCodons,
                                              __chunksequence(sequence[i:], 3, True)))
            if len(seq) < minLength: continue
            if negStrand:
                name = '%s_to_%s' % (seqLen-i, seqLen-i-len(seq)+1)
                desc = '(negative strand)'
            else:
                name = '%i_to_%i' % (i+1, i+len(seq))
                desc = ''
            orfs.append(Sequence(name=name, description=desc, sequence=seq))
    return orfs

def calcIdentity(sequence1, sequence2):
    matches, total = 0, 0
    for c1, c2 in zip(sequence1, sequence2):
        if c1 == c2:
            if c1 == '-': continue
            matches += 1
        total += 1
    numStr = '%i / %i' % (matches, total)
    percent = float(matches) / total * 100
    return percent, numStr
            

# # # # # # # # # #  Private Functions  # # # # # # # # # #
def __chunksequence(sequence, chunksize, only_complete=False):
    if only_complete: length = len(sequence)-(chunksize-1)
    else: length = len(sequence)
    for i in xrange(0, length, chunksize):
        yield sequence[i:i+chunksize]

def parsefasta_OLD(textobj):
    """Returns a list of Sequence objects from the lines or file object.
    If onlyThese is a list of strings, only those sequences that start
    with a string in that list will be collected."""
    seqs, buff = [], []
    for line in textobj:
        if line.startswith('>'):
            if buff:
                if not seqs:
                    seqs.append(Sequence())
                seqs[-1].seq = ''.join(buff)
                buff = []
            nm, _, desc = line[1:].strip().partition(' ')
            seqs.append(Sequence( name=nm, description=desc ))
        else:
            line = line.strip()
            if line == '': continue
            buff.append(line)
    if buff:
        if seqs: seqs[-1].seq = ''.join(buff)
        else: seqs.append(Sequence(name='', sequence=''.join(buff)))
    return seqs


class Sequence(object):
    """DocString

    seq is an alias for the sequence attribute.
    Setting seq or sequence involves several filtering steps. If speed is an issue,
    the _sequence attribute can be used which bypasses these.
    """
    def __init__(self, name='Unnamed sequence', description='', sequence='', allowedChars='-_*?', onlyUpper=True):
        self.allowedChars = allowedChars; self.onlyUpper = onlyUpper
        self.name = name
        self.description = description
        self.sequence = sequence
        
    # # # # #  Public Methods  # # # # #
    def translate(self):
        return translate(self)
    def invcomplement(self):
        return invcomplement(self)
    def fasta(self, line=60, spaces=True, numbers=True):
        return '>%s\n%s\n' % (self.header, fasta(self, line, spaces, numbers))
    def append(self, sequence):
        self.seq += sequence
    def isNucleotide(self):
        if set(filter(str.isalpha, self.seq)) <= set(['A', 'C', 'G', 'T', 'U']): return True
        return False

    # # # # #  Public String Methods  # # # # #
    def count(self, sub, *args):
        if self.onlyUpper: sub = sub.upper()
        return self.sequence.count(sub, *args)
    def startswith(prefix, *args):
        if self.onlyUpper: prefix = prefix.upper()
        return self.sequence.startswith(prefix, *args)
    def endswith(suffix, *args):
        if self.onlyUpper: suffix = suffix.upper()
        return self.sequence.endswith(suffix, *args)
    def find(sub, *args):
        if self.onlyUpper: sub = sub.upper()
        return self.sequence.find(sub, *args)
    def rfind(sub, *args):
        if self.onlyUpper: sub = sub.upper()
        return self.sequence.rfind(sub, *args)
    def index(sub, *args):
        if self.onlyUpper: sub = sub.upper()
        return self.sequence.index(sub, *args)
    def rindex(sub, *args):
        if self.onlyUpper: sub = sub.upper()
        return self.sequence.rindex(sub, *args)

    # # # # #  Private Methods  # # # # #
    
    # # # # #  Under-the-hood Methods  # # # # #
    def __getName(self): return self.__name
    def __setName(self, name):
        if not name: name = ''
        else: name = str(name)
        name = name.strip()
        if name.startswith('>'): name = name[1:]
        self.__name = name
    def __getDescription(self): return self.__description
    def __setDescription(self, desc):
        if not desc: desc = ''
        else: desc = str(desc)
        desc = desc.strip()
        self.__description = desc
    def __getHeader(self):
        header = self.name
        if self.description: header += ' %s' % self.description
        return header
    def __getSequence(self): return self.__sequence
    def __setSequence(self, sequence):
        if not sequence: sequence = ''
        sequence = filter(lambda c: c in self.allowedChars or c.isalpha(), sequence)
        if self.onlyUpper: sequence = sequence.upper()
        self.__sequence = ''.join(sequence)
    def __getFastseq(self): return self.__sequence
    def __setFastseq(self, sequence): self.__sequence = sequence
    name = property(__getName, __setName)
    description = property(__getDescription, __setDescription)
    header = property(__getHeader)
    sequence = property(__getSequence, __setSequence)
    seq = sequence
    _sequence = property(__getFastseq, __setFastseq)
    def __len__(self): return len(self.__sequence)
    def __getitem__(self, key): return self.sequence[key]
    def __setitem__(self, key, value):
        s = list(self.sequence)
        s[key] = value
        self.sequence = ''.join(s)
    def __delitem__(self, key):
        s = list(self.sequence)
        del s[key]
        self.sequence = ''.join(s)
    def __iter__(self): return iter(self.sequence)
    def __reversed__(self): return reversed(self.sequence)
    def __contains__(self, item): return ''.join(item) in self.sequence
    def __repr__(self): return '%s: %s' % (self.header, self.sequence)
    def __str__(self): return '>%s\n%s\n' % (self.header, self.sequence)
    def __eq__(self, other):
        if type(other) == str: s = other
        else: s = other.sequence
        return self.sequence == s
    def __ne__(self, other):
        return not self == other
