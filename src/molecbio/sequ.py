#!/usr/bin/python
# Author: Dave Curran
# Date:   April 2012
# # # # # # # # # # #

# TODO:
# Replace this with work/antigen_analysis/sequ.py once it's done. Probably a good idea to keep this file around, maybe as something like sequ_legacy.py so it can still be accessed.
# When rewritten, would be great to also have a command line tool to do many common tasks.
# - Like to remove empty sequence. Or remove non-unique sequences. Extract a portion of an alignment, etc.

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
    -- fasta(sequence, line=60, spaces=False, numbers=False) -- Converts some
       sequence into a fasta format, returns as a string.
    -- parsefasta(textobj) -- Takes a file object or list of text lines and
       returns a list of Sequence objects.
    -- loadfasta(filepath, onlyThese=[]) -- Takes a file path string, returns list of Sequences.
    -- savefasta(seqList, filepath) -- Saves list of Sequences to filepath.
    -- cleanfasta(filepath) -- Overwrites sequences at filepath, formatting them.
    -- filter_unique(seqs, to_keep=None, duplicate_names=False, return_replaced=False, compare_gaps=False, compare_terminal_stop=False) -- Takes a list of Sequence objects, filtering out repeated sequences.

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
def fasta(sequence, line=60, spaces=False, numbers=False):
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
def loadfasta(filepath, onlyThese=[]):
    """Returns a list of Sequence objects from the filepath."""
    with open(filepath, 'r') as f:
        seqs = parsefasta(f, onlyThese)
    return seqs
def savefasta(seqList, filepath, line=60, spaces=False, numbers=False):
    """Saves the given list of Sequences to filepath in fasta format."""
    seqStr = '\n'.join(seq.fasta(line, spaces, numbers) for seq in seqList)
    with open(filepath, 'w') as f:
        f.write(seqStr)
def cleanfasta(filepath):
    seqs = loadfasta(filepath)
    buff = [seq.__str__() for seq in seqs]
    #buff = map(str, seqs)
    with open(filepath, 'w') as f:
        f.write('\n'.join(buff))
        f.close()

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

def filter_unique(aln, to_keep=None, duplicate_names=False, return_replaced=False, compare_gaps=False, compare_terminal_stop=False, **kwargs):
    """Sequentially parses aln (a list of Sequence objects, aligned or not), filtering out any with a sequence that has been seen before. Returns a list of filtered Sequence objects in the same order as in aln (except that sequences in `to_keep` will be moved to the location of the first identical non-to_keep sequence.
    `to_keep` can be None or a sequence of strings representing sequence names to keep. If given, any such sequences will never be filtered whether their sequences are unique or not. These sequences will also be preferentially used to replace other identical sequences, over sequences not in `to_keep`. If none of a set of identical sequences are in `to_keep`, the one found earliest in the list will be kept.
    If `duplicate_names` is False, a sequence object will be completely ignored if another sequence has already been encountered with the same name, no matter the sequence itself. This will also apply to any names in `to_keep`. If True, this second sequence will only be kept if its sequence is different from other sequences (unless it is in `to_keep`); if multiple sequences with the same name are to be kept, a "_DUPLICATE_X" counter tag will be added to the sequence name in the returned `replaced` dict (but not in the returned list of Sequence objects, so they may not match).
    If `return_replaced` is True, this function will return a list of filtered Sequence objects and a dict. This dict describes the filtered sequences: {'kept_name':['filtered_name1', 'filtered_name2', ...], ...}. Every filtered sequence name will be present once in the values of this dict (except for sequences with identical names if `duplicate_names` is False), but a kept sequence will only be in the keys if it was identical to at least 1 filtered sequence.
    If `compare_gaps` is False, gap characters will be removed for the sake of the comparison, but will still be present in the returned sequence objects. If True, the sequences will be compared exactly as they are.
    If `compare_terminal_stop` is False, a single terminal "_" or "*" character will be removed for the sake of comparison, but will still be present in the returned sequence objects. If True, the sequences will be compared exactly as they are.
    This function may safely be called with a dict of **kwargs containing unused arguments."""
    # Formatting functions
    def format_duplicate_name(name, cntr):
        return f'{name}_DUPLICATE_{cntr}'
    def format_sequence(seq_str):
        if compare_gaps == False:
            seq_str = seq_str.replace('-', '') # Remove gaps in case of ambiguous placement
        if compare_terminal_stop == False and seq_str[-1] in ('_', '*'):
            seq_str = seq_str[:-1]
        return seq_str
    
    # Format argument
    if not to_keep:
        to_keep = set()
    else:
        to_keep = set(to_keep)
    
    # Get the important sequences first
    names_cntr = {} # To track unique names
    seqdict = {} # {"SEQUENCE":(kept_name, [filtered_name1, ...]), ...}
    if to_keep:
        for seq in aln:
            seqname = seq.name
            if seqname not in to_keep:
                continue
            # Unique names
            dup_name = seqname in names_cntr
            if dup_name:
                if duplicate_names == False:
                    continue
            else:
                names_cntr[seqname] = 1
            # Get comparison sequence
            compseq = format_sequence(seq.seq)
            # Check for comparison sequence uniqueness
            if compseq not in seqdict: # First occurance
                if dup_name:
                    dup_name = format_duplicate_name(seqname, names_cntr[seqname])
                    seqdict[compseq] = (dup_name, []) # List tracks filtered seq names with same compseq
                    names_cntr[seqname] += 1
                else:
                    seqdict[compseq] = (seqname, [])

    # Identify unique & redundant sequences, filling out new_aln list
    names_cntr = {} # Reset name counter
    new_aln = [] # Final list of filtered sequence objects to return
    for seq in aln:
        seqname = seq.name
        # Unique names
        dup_name = seqname in names_cntr
        if dup_name:
            if duplicate_names == False:
                continue
        else:
            names_cntr[seqname] = 1
        
        # Important sequence
        if seqname in to_keep:
            new_aln.append(seq) # Already processed above; don't need to count name
            continue

        # Get comparison sequence
        compseq = format_sequence(seq.seq)
        # Check for comparison sequence uniqueness
        if compseq not in seqdict: # First occurance, keep it
            if dup_name:
                dup_name = format_duplicate_name(seqname, names_cntr[seqname])
                seqdict[compseq] = (dup_name, []) # List tracks filtered seq names with same compseq
                names_cntr[seqname] += 1
            else:
                seqdict[compseq] = (seqname, [])
            new_aln.append(seq)
        else: # Not first occurance, filter it
            if dup_name:
                dup_name = format_duplicate_name(seqname, names_cntr[seqname])
                seqdict[compseq][1].append(dup_name)
                names_cntr[seqname] += 1
            else:
                seqdict[compseq][1].append(seqname)
    
    # Return final objects
    if return_replaced == True:
        # Fill out replaced dict
        replaced = {} # Record which sequences were filtered
        for (kept_name, filtered_names) in seqdict.values():
            if len(filtered_names) == 0: # Unique sequence, nothing filtered
                continue
            replaced[kept_name] = filtered_names
        return new_aln, replaced
    else:
        return new_aln


# # # # # # # # # #  Private Functions  # # # # # # # # # #
def __chunksequence(sequence, chunksize, only_complete=False):
    if only_complete: length = len(sequence)-(chunksize-1)
    else: length = len(sequence)
    for i in range(0, length, chunksize):
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
    def fasta(self, line=60, spaces=False, numbers=False):
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
        sequence = ''.join(filter(lambda c: c in self.allowedChars or c.isalpha(), sequence))
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
