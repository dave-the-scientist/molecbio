#!/usr/bin/python
# Author: Dave Curran
# Date:   April 2012
# # # # # # # # # # #
"""
A module made to manipulate and perform several functions on sequence alignments.

Defines the following:
    I/O Functions:
    -- open_clustal(filepath) -- Reads in a clustal-formatted alignment
    file, returning a list of sequ.Sequence objects.
    -- open_fasta(filepath) -- Reads in a fasta-formatted alignment file,
    returning a list of sequ.Sequence objects.
    -- saveas_clustal(seqList, filepath, nameWith=18, seqWidth=60) -- Saves
    the given list of sequ.Sequence objects as a clustal-formatted file.
    -- saveas_fasta(seqList, filepath) -- Saves the given list of sequ.Sequence
    objects as a clustal-formatted file.
    -- saveas_binned(filepath, scores, thresholds=(3.33,6.66)) -- Divides the
    given list of scores into bins specified by thresholds, and saves them as
    a tab-separated null-padded file.

    Data Manipulation Functions:
    -- fasta_to_clustal(seqList, nameWidth=18, seqWidth=60) -- Takes a list
    of Sequence objects, returning a clustal-formatted string.
    -- sliding_average(values, size=3) -- Returns a list of average values
    calculated for a sliding window of the given size.
    -- bin_values(values, thresholds=(3.33, 6.66)) -- Splits the list of values
    into bins specified by thresholds.

    Alignment Functions:
    -- quality(seqList, qualityCalc='sd', normCalc='new', maxVal=10.0,
               blosumD=blosum.blosum62) -- This returns a list of quality scores
    for the given sequence objects.
    -- protein_align_to_dna(proteinAlnFile, dnaFastaFile) -- Takes filepaths for
    a protein alignment and DNA sequences, returns a list of Sequence objects
    containing the aligned DNA codons.
    -- seqIdentity(seq1, seq2) -- Does not perform an alignment, but calculates the
    percent identity between 2 sequences. Returns a float and a string, the identity
    and the matches / total.
    -- pairwiseIdentity(seqList) -- Calculates percent identity between all pairs of
    sequences in seqList, returning the results as a string.
"""
from __future__ import with_statement # Needed for python 2.5
import itertools, math, operator
import sequ, blosum

# # # # #  I/O Functions  # # # # #
def open_clustal(filepath):
    """Parses a clustal alignment file, returning a list of Sequences."""
    names, seqs = __parseClustal(filepath)
    if not names or not seqs: return []
    return [sequ.Sequence(name, sequence=''.join(seqs[name])) for name in names]
def open_fasta(filepath):
    """Returns a list of Sequence objects from the filepath."""
    seqs = []
    with open(filepath, 'rb') as f:
        seqs = sequ.parsefasta(f)
    return seqs

def saveas_clustal(seqList, filepath, nameWidth=18, seqWidth=60):
    """Saves the given list of Sequences to filepath, in clustal format.

    The nameWidth argument controls how many characters to use from each name,
    while seqWith controls how many sequence characters are written per line.
    Both have default values."""
    seqStr = fasta_to_clustal(seqList, nameWidth, seqWidth)
    with open(filepath, 'wb') as f:
        f.write(seqStr)
def saveas_fasta(seqList, filepath):
    """Saves the given list of Sequences to filepath in fasta format."""
    seqStr = '\n'.join(seq.fasta() for seq in seqList)
    with open(filepath, 'wb') as f:
        f.write(seqStr)
def saveas_phylip(seqList, filepath):
    """Formats for use with PhyML, which is a little different than
    standard PHYLIP format when it comes to name conventions."""
    perLine = 80
    seqLen = len(seqList[0].seq)
    nameLen = 0
    names = []
    sequences = []
    for seq in seqList:
        name = seq.header
        name = '_'.join(name.strip().split())
        name = name.replace('(','').replace(')','').replace(',','').replace(':','').replace('.','')
        names.append(name)
        if len(name) > nameLen: nameLen = len(name)
        sequences.append([seq.seq[i:i+perLine] for i in range(0,seqLen,perLine)])
    nameLen = min(nameLen, 100)
    fmtstr = '%%-%is %%s' % nameLen
    for i, name in enumerate(names):
        sequences[i][0] = fmtstr % (name[:nameLen], sequences[i][0])
    sequences = ['\n'.join(block)+'\n' for block in zip(*sequences)]
    buff = ['%i %i' % (len(names), seqLen)]
    buff.extend(sequences)
    with open(filepath, 'wb') as f:
        f.write('\n'.join(buff))
    
def saveas_binned(filepath, scores, thresholds=(3.33333, 6.66666)):
    """Saves list of scores as spreadsheet-readable bins.

    The scores will be divided up into the bins specified by the thresholds
    tuple. Each bin will be null-padded and tab-separated, and saved on its
    own line the the given filepath."""
    nan = "#N/A" # Null value for MS Excel.
    scores = bin_values(scores, thresholds)
    f = open(filepath, 'wb')
    f.write('\n'.join('\t'.join('%.3f'%num if num else nan for num in binn) for binn in scores))
    f.close()

# # # # #  Data Manipulation Functions  # # # # #
def fasta_to_clustal(seqList, nameWidth=18, seqWidth=60):
    segs = __fasta_to_segs(seqList, nameWidth, seqWidth)
    return '\n\n'.join(segs)
def sliding_average(values, size=3):
    """Calculates sliding average over the given list."""
    def window(seq, size):
        """Posted by Daniel DiPaolo on StackOverflow.com"""
        it = iter(seq)
        result = tuple(itertools.islice(it, size))
        if len(result) == size: yield result    
        for elem in it:
            result = result[1:] + (elem,)
            yield result
    def avg(nums):
        return sum(nums) / size
    size = float(size)
    return map(avg, window(values, size))

def bin_values(values, thresholds=(3.33333, 6.66666)):
    """Divides the list of values into the bins specified by thresholds.
    
    Each bin will be zero-padded to the same length as the original, and
    excludes the upper boundary value. So the default bins are [-inf, 3.333),
    [3.333, 6.666), [6.666, +inf]."""
    def binNum(num):
        for i, thresh in enumerate(thresholds):
            if num < thresh: break
        else: i += 1
        return nans[:i] + (num,) + nans[i+1:]
    nans = (0.0,) * (len(thresholds)+1)
    l = (binNum(val) for val in values)
    return zip(*l)

# # # # #  Alignment Functions  # # # # #
def quality(seqList, qualityCalc='sd', normCalc='new', maxVal=10.0,
            blosumD=blosum.blosum62):
    """Calculates the quality score for the given alignment.

    Returns a list of quality scores for the list of sequ.Sequence objects given
    as seqList. The qualityCalc argument chooses the quality algorithm, normCalc
    specifies the normalization calculation, maxVal allows different alignments
    to be compared against eachother by putting the scores in the same range, and
    blosumD allows a different BLOSUM matrix to be used.

    A quality score is a measure of the differences of the amino acids at one
    position in a protein alignment. The traditional way considers each residue
    as a 20-dimensional vector, equal to the BLOSUM score for that residue against
    all others. At each column the quality score is actually the absolute mean
    deviation (MD) of those vectors. Passing 'md' as the metric argument will use
    this original algorithm, while passing 'sd' or nothing will use my new method
    that calculates the standard deviation instead.

    The normCalc allows different algorithms to be used when normalizing the data.
    The 'old' way is to find the difference between the score and the maximum score,
    and multiply by the fraction of residues that are not gaps. The new way is almost
    the same, but it uses the square of the fraction of not-gaps for a position. This
    is to lessen the penalty of having a gap. Finally, passing anything other than
    'old' or 'new' will use no normalization, and return the raw score values."""
    columns = [column for column in itertools.izip(*seqList) if len(filter(str.isalpha, column)) > 0]
    calcQuality = __chooseQualityCalc(qualityCalc, blosumD)
    normalize = __chooseNormCalc(normCalc, maxVal, columns)
    scores = map(calcQuality, columns)
    return normalize(scores)

def protein_align_to_dna(proteinAlnFile, dnaFastaFile):
    """Aligns DNA codons to a protein sequence alignment.

    Pass in the filepaths to a protein alignment (either in clustal or fasta format)
    and a DNA sequence file (in fasta format). Every protein sequence must have its
    corresponding DNA sequence, or this will return False. A list of sequ.Sequence
    objects is returned with the same names and in the same order as the protein
    sequences from the file. These will have a DNA alignment as their sequence, that
    is organized into codons."""
    protSeqs = open_clustal(proteinAlnFile)
    if not protSeqs: protSeqs = open_fasta(proteinAlnFile)
    dnaSeqs = open_fasta(dnaFastaFile)
    if not protSeqs or not dnaSeqs or len(protSeqs) > len(dnaSeqs): return False
    matchedSeqs = __matchProtDna(protSeqs, dnaSeqs)
    alnDnaSeqs = __generateDnaAlign(matchedSeqs)
    return alnDnaSeqs

def seqIdentity(seq1, seq2):
    """ """
    matches, total = 0, 0
    for c1, c2 in itertools.izip(seq1, seq2):
        if c1 == c2:
            if c1 == '-': continue
            matches += 1
        total += 1
    numStr = '%i / %i' % (matches, total)
    percent = float(matches) / total * 100
    return percent, numStr

def pairwiseIdentity(seqList):
    """ """
    names = []
    buff = ['      ' + ' '.join('%5i' % i for i in reversed(xrange(1, len(seqList))))]
    for i, seq1 in enumerate(seqList[:-1]):
        buff.append('%5i %s\n' % (i, ' '.join('%5.1f' % seqIdentity(seq1,seq2)[0] for seq2 in seqList[:i:-1])))
        names.append('%i: %s' % (i, seq1.name))
    names.append('%i: %s' % (i+1, seqList[-1].name))
    buff.append('\n\nID numbers:')
    buff.extend(names)
    return '\n'.join(buff)

# # # # #  Private Functions  # # # # #
def __parseClustal(filename):
    f = open(filename, 'rb')
    line = f.readline()
    if not line.startswith('CLUSTAL'): return (False, False)
    for line in f:
        if line != '\n' and not line.startswith('CLUSTAL'):
            break
    line = line.split(); name = line[0]; seq = line[1]
    names = [name]
    d = {name:[seq]}
    for line in f:
        if line[0].isspace(): continue
        line = line.split(); name = line[0]; seq = line[1]
        if name not in d:
            names.append(name)
            d[name] = [seq]
        else:
            d[name].append(seq)
    f.close()
    return names, d
def __fasta_to_segs(seqList, nameWidth, seqWidth):
    strong = (set('STA'), set('NEQK'), set('NHQK'), set('NDEQ'), set('QHRK'),
              set('MILV'), set('MILF'), set('HY'), set('FYW'))
    weak = (set('CSA'), set('ATV'), set('SAG'), set('STNK'), set('STPA'),
            set('SGND'), set('SNDEQK'), set('NDEQHK'), set('NEQHRK'),
            set('FVLIM'), set('HFY'))
    def lineGen(seq):
        s = '%%-%is' % nameWidth
        name = s % seq.name[:nameWidth]
        total = 0
        for i in xrange(0, len(seq), seqWidth):
            s = seq[i:i+seqWidth]
            total += (len(s) - s.count('-'))
            yield '%s %s %i' % (name, s, total)
    def consGen():
        def nucConsSymbol(column):
            if '-' not in column and len(set(column)) == 1:
                return '*'
            else: return ' '
        def residueConsSymbol(column):
            if '-' in column: return ' '
            c = set(column)
            if len(c) == 1: return '*'
            for s in strong:
                if s >= c: return ':'
            for w in weak:
                if w >= c: return '.'
            return ' '
        s = ' ' * nameWidth
        buff = []
        if allNucleotides: it = map(nucConsSymbol, itertools.izip(*seqList))
        else: it = map(residueConsSymbol, itertools.izip(*seqList))
        for i in xrange(0, len(seqList[0]), seqWidth):
            yield '%s %s' % (s, ''.join(it[i:i+seqWidth]))
    allNucleotides = all(map(sequ.Sequence.isNucleotide, seqList))
    buff = ['CLUSTAL W multiple sequence alignment']
    gens = map(lineGen, seqList)
    gens.append(consGen())
    for seg in itertools.izip(*gens):
        buff.append('\n'.join(seg))
    buff.append('\n')
    return buff
def __chooseQualityCalc(qualityCalc, blosumD):
    """Both calcMD and calcSD are optimized for speed, so it may be difficult to see
    exactly what is being calculated. Especially with calcSD, where the SD
    formula has been algebreaically manipulated to be more efficient.
    Original MD calculations were adapted from
    http://bips.u-strasbg.fr/fr/Documentation/ClustalX/#Q"""
    alpha = blosumD.alphabet
    srs = dict((res, [blosumD[res, r] for r in alpha]) for res in alpha)
    def calcMD(column):
        """Calculates the absolute mean deviation of the S-scores for a column.

        The average distance each residue lies from the mean point, in R-
        dimensional space."""
        def dist(sr):
            return math.sqrt(sum(dr*dr for dr in itertools.imap(operator.sub, sr, X)))
        column = filter(str.isalpha, column)
        n = float(len(column))
        ss = [srs[res] for res in column]
        X = [sum(dimension)/n for dimension in itertools.izip(*ss)]
        return sum(dist(sr) for sr in ss) / n
    def calcSD(column):
        """Calculates the standard deviation of the S-scores for a column.

        An indication of the spread of the residues in R-dimensional space. The
        formula for standard deviation has been algebraically manipulated to
        optimize it for speed."""
        column = filter(str.isalpha, column)
        n = float(len(column))
        ss = [srs[res] for res in column]
        tot, tot2 = 0, 0
        for dimension in itertools.izip(*ss):
            tot += sum(dimension)**2
            tot2 += sum(s*s for s in dimension)
        var = (tot2 - tot / n) / n
        return math.sqrt(var)
    if qualityCalc == 'sd': return calcSD
    elif qualityCalc == 'md': return calcMD
    else: raise ValueError("Invalid argument for quality function. The 'metric' argument must be set to 'sd' or 'md'.")
def __chooseNormCalc(normCalc, maxVal, columns):
    """Normalizes to number of gaps and to max score."""
    def newNorm(scores):
        def norm(d, column):
            g = column.count('-')
            return (maxScore - d) * ( 1 - (g/num) * (g/num) ) * x
        maxScore = max(scores)
        x = maxVal / maxScore if maxScore else 1.0
        num = float(len(columns[0]))
        return [norm(d, column) for d, column in
                itertools.izip(scores, columns)]
    def oldNorm(scores):
        def norm(d, column):
            g = column.count('-')
            return (maxScore - d) * ( 1 - (g/num) ) * x
        maxScore = max(scores)
        x = maxVal / maxScore if maxScore else 1.0
        num = float(len(columns[0]))
        return [norm(d, column) for d, column in
                itertools.izip(scores, columns)]
    def noNorm(scores): return scores
    if normCalc == 'new': return newNorm
    elif normCalc == 'old': return oldNorm
    else: return noNorm

def __matchProtDna(protSeqs, dnaSeqs):
    """Goes through the protSeqs, finding the corresponding DNA sequence in
    dnaSeqs. Returns the list of protSeqs, where each object now has the
    attribute 'Sequence.dnaSequence'."""
    dnaSeqs = dnaSeqs[:]
    for protSeq in protSeqs:
        matchedSeq = None
        cleanSeq = filter(str.isalpha, protSeq.seq)
        for dnaSeq in dnaSeqs:
            if cleanSeq in dnaSeq.translate():
                if dnaSeq.name == protSeq.name:
                    protSeq.dnaSequence = dnaSeq.seq
                    dnaSeqs.remove(dnaSeq)
                    break
                if not matchedSeq: matchedSeq = dnaSeq
        if matchedSeq:
            protSeq.dnaSequence = matchedSeq.seq
            dnaSeqs.remove(matchedSeq)
        else: return False
    return protSeqs
def __generateDnaAlign(matchedSeqs):
    """Takes a list of Sequence objects for a protein alignment, that all
    have an attribute dnaSequence. A list of new sequence objects are
    returned, with the DNA sequence aligned to the protein alignment, on
    a codon-codon basis."""
    seqs = []
    for seq in matchedSeqs:
        cleanPep = filter(str.isalpha, seq.seq)
        trans = sequ.translate(seq.dnaSequence)
        i = trans.find(cleanPep)
        dnaSeq = seq.dnaSequence[i*3:]
        subIters = [iter(dnaSeq)] * 3
        dnaCodons = itertools.izip(*subIters)
        dnaAln = []
        for res in seq.seq:
            if res.isalpha(): dnaAln.extend(dnaCodons.next())
            elif res == '-': dnaAln.append('---')
        seqs.append(sequ.Sequence(name=seq.name, sequence=''.join(dnaAln)))
    return seqs
