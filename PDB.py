#!/usr/bin/python
# Author: Dave Curran
# Date:   April 2012
# # # # # # # # # # #
"""Several functions for dealing with and modifying .pdb files.

     Functions come in two kinds, those that automatically resave the data after
some manipulation, and those that return the data as a list (which have the same
name, but ending in 'data'). These data functions allow you to chain together
several functions, starting with parsepdb(), before saving the eventual output
using savepdb().

     This script can also be run from the command line. Running with no options
will display the help message.

Variables:
    -- trans3to1 -- Dict for changing 3-letter residue codes to 1-letter.

Functions:
    -- parsepdb(filename) -- Opens the filename, returning list of lines.
    -- savepdb(pdbdata, filename) -- Saves the data to the filename.
    -- cleanpdb(original, new, noH=False, prefixes=['ATOM','TER']) -- Filters extraneous lines.
    -- cleanpdbdata(pdbdata, noH=False, prefixes=['ATOM','TER']) -- Filters extraneous lines.
    -- changechainID(original, new, oldID, newID) -- Change the chainID of some chain.
    -- changechainIDdata(pdbdata, oldID, newID) -- Change the chainID of some chain.
    -- pdbtofasta(pdbFile, fastaFile, fillerChar='-') -- Extract protein sequence from pdb.
    -- resfromdata(pdbdata, fillerChar='-') -- Extract protein sequence from pdb.
    -- tripresfromdata(pdbdata, fillerChar='-') -- Extract protein sequence from pdb.
    -- renumber(original, new, chainID, difference) -- Renumber the residues.
    -- renumberdata(pdbdata, chainID, difference) -- Renumber the residues.
    -- mapconservation(pdbfile, alignmentfasta, newpdbfile) -- Replace the b-factors with alignment quality.
    -- pdbqualityscores((pdbdata, alignmentfasta)) -- Replace the b-factors with alignment quality.
"""
import os, sys
# # # # # # # # # #  Variables  # # # # # # # # # #
"""Dictionary to change 3-letter amino acid codes to 1-letter."""
trans3to1 = {
    'ALA':'A', 'ARG':'R', 'ASN':'N', 'ASP':'D',
    'CYS':'C', 'GLU':'E', 'GLN':'Q', 'GLY':'G',
    'HIS':'H', 'ILE':'I', 'LEU':'L', 'LYS':'K',
    'MET':'M', 'PHE':'F', 'PRO':'P', 'SER':'S',
    'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'
    }

# # # # # # # # # #  Functions  # # # # # # # # # #
def parsepdb(filename):
    """Reads the file with readlines(), returns the list."""
    f = open(filename, 'rb')
    lines = f.readlines()
    f.close()
    return lines

def savepdb(pdbdata, filename):
    """Writes the data to the filename, ensuring it has an END."""
    if not pdbdata: return False
    if not pdbdata[-1].startswith('END'): pdbdata.append('END\n')
    f = open(filename, 'wb'); f.write(''.join(pdbdata)); f.close()
    return True

def cleanpdb(original, new, noH=False, prefixes=['ATOM','TER']):
    """Filters the original .pdb file and saves it as 'new'. Returns boolean runcode.

    noH -- Remove all hydrogen atoms (default False).
    prefixes -- A list of strings indicating lines to keep (default ['ATOM','TER']).
    """
    if not os.path.isfile(original): return False
    f = open(original, 'rb')
    lines = cleanpdbdata(f, noH, prefixes)
    f.close()
    return savepdb(lines, new)

def cleanpdbdata(pdbdata, noH=False, prefixes=['ATOM','TER']):
    """Filters an iterable of lines from a pdb file. Returns a list of strings."""
    if noH:
        def renumberAtom(line, num):
            cur = line[6:11]
            if not cur.strip().isdigit(): return line
            cur = int(cur)
            if cur == num: return line
            return line[:6] + '%5i'%num + line[11:]
        lines, i = [], 1
        for line in pdbdata:
            if not any(line.startswith(pref) for pref in prefixes):
                continue
            if line[13] != 'H':
                lines.append(renumberAtom(line, i))
                i+=1
    else:
        lines = [line for line in pdbdata if any(line.startswith(pref) for pref in prefixes)]
    return lines

def changechainID(original, new, oldID, newID):
    """Changes the chain ID from 'oldID' to 'newID', and saves the pdb file as 'new'."""
    if not os.path.isfile(original): return False
    f = open(original, 'rb')
    lines = changechainIDdata(f, oldID, newID)
    f.close()
    return savepdb(lines, new)

def changechainIDdata(pdbdata, oldID, newID):
    """Modifies the chain ID from 'oldID' to 'newID' in the iterable pdbdata, returning a list."""
    if oldID == '': oldID = ' '
    if newID == '': newID = ' '
    if len(oldID) != 1 or len(newID) != 1:
        return False
    prefixes = ['ATOM  ', 'ANISOU', 'HETATM', 'TER   ']
    lines = []
    for line in pdbdata:
        if line[:6] in prefixes:
            if line[21:22] == oldID:
                line = line[:21] + newID + line[22:]
        lines.append(line)
    return lines

def pdbtofasta(pdbFile, fastaFile, fillerChar='-'):
    """Extracts the amino acid sequence from the pdb and saves it as fastaFile.

    Each chain in the pdb file will be its own entry in the fasta formatted file.
    Unknown or missing residues will be represented by fillerChar."""
    f = open(pdbFile, 'rb')
    seqs = resfromdata(f, fillerChar)
    f.close()
    seqName = os.path.basename(fastaFile)
    if '.' in seqName: seqName = seqName.rpartition('.')[0]
    f = open(fastaFile, 'wb')
    for chain in seqs:
        if len(seqs) == 1: name = '>' + seqName + '\n'
        else: name = '>' + seqName + '_%s\n' % chain[0]
        seq = ''.join(chain[1:]) + '\n\n'
        f.write(name + seq)
    f.close()

def resfromdata(pdbdata, fillerChar='-'):
    """Parses the iterable 'pdbdata' and returns the sequence as a list of lists.

    Each sublist represents one protein chain, with the identifier as the first
    entry, and each residue as a one letter code. Unknown or missing residues
    are represented by the fillerChar."""
    trips = tripresfromdata(pdbdata, fillerChar)
    seq = []
    for chain in trips:
        l = [chain[0]]
        for res in chain[1:]:
            l.append(trans3to1.get(res, fillerChar))
        seq.append(l)
    return seq

def tripresfromdata(pdbdata, fillerChar='-'):
    """Parses the iterable 'pdbdata' and returns the sequence as a list of lists.

    Each sublist represents one protein chain, with the identifier as the first
    entry, and each residue as a three letter code. Unknown or missing residues
    are represented by the fillerChar."""
    curChain, curNum, seq = '', 0, []
    for line in pdbdata:
        if line.startswith('ATOM'):
            chain = line[21]
            num = int(line[22:26])
            if chain != curChain:
                seq.append([chain])
                curChain, curNum = chain, 0
            if curNum == num:
                continue
            if curNum < num - 1:
                diff = num - 1 - curNum
                seq[-1].extend(fillerChar*diff)
                curNum += diff
            elif curNum == num - 1:
                res = line[17:20].strip().upper()
                seq[-1].append(res)
                curNum = num
            else:
                return False
    return seq

def renumber(original, new, chainID, difference):
    """Changes the residue numbers by adding 'difference' to each, saving as 'new'.

    If the initial residue is number 25, pass -24 as the difference to change it to
    1. You must pass the single letter chainID code to specify which pdb chain to
    modify."""
    f = open(original, 'rb')
    lines = renumberdata(f, chainID, difference)
    f.close()
    savepdb(lines, new)

def renumberdata(pdbdata, chainID, difference):
    """Changes the residue numbers by adding 'difference' to each. Returns a list.

    If the initial residue is number 25, pass -24 as the difference to change it to
    1. You must pass the single letter chainID code to specify which pdb chain to
    modify."""
    prefixes = ['ATOM  ', 'ANISOU', 'HETATM', 'TER   ']
    lines = []
    for line in pdbdata:
        if line[:6] in prefixes:
            chain = line[21:22]
            if chain == chainID:
                num = int(line[22:26])
                newNum = '%4i' % (num + difference)
                line = line[:22] + newNum + line[26:]
        lines.append(line)
    return lines

def mapconservation(pdbfile, alignmentfasta, newpdbfile):
    """Replaces the b-factors of the pdbfile with the quality scores from the
    given alignment, saving the structure as newpdbfile."""
    pdbdata = parsepdb(pdbfile)
    chain, scores = pdbqualityscores(pdbdata, alignmentfasta)
    buff = []
    scoresIter = iter(scores)
    maxScore, prevRes, curQual = max(scores), -1, 10.0
    for line in pdbdata:
        if not line.startswith('ATOM'):
            buff.append(line)
            continue
        curChain = line[21]
        if curChain != chain:
            buff.append(line)
            continue
        resNum = int(line[22:26])
        if resNum > prevRes:
            prevRes = resNum
            try:
                curQual = maxScore - scoresIter.next()
            except StopIteration:
                print('Error: the pdb sequence is longer than the calculated alignment scores. This usually means the pdb sequence is not present in the alignment file.')
                return
        buff.append('%s%6.2f%s' % (line[:60], curQual, line[66:]))
    savepdb(buff, newpdbfile)

def pdbqualityscores(pdbdata, alignmentfasta):
    """Returns a list of floats, where each is the calculated alignment quality
    of one of the residues of the pdb structure. Only the first chain in the
    structure is modified."""
    import align, aligners
    pdbseq = resfromdata(pdbdata, '-')[0]
    pdbchain = pdbseq[0]
    pdbseq = ''.join(pdbseq[1:])
    seqs = align.open_fasta(alignmentfasta)
    quality = align.quality(seqs)
    aligner = aligners.Needleman()
    identity, sequence = 0, None
    identities = aligner.one_vs_all_identity(pdbseq, [seq.seq for seq in seqs])
    for iden, seq in zip(identities, seqs):  # Find closest sequence in alignment
        if iden > identity:
            identity = iden
            sequence = seq
    aln = aligner.align(pdbseq, sequence.seq)
    pdbseqIter, seqIter = iter(aln[0].replace('?','-')), iter(aln[1].replace('?','-'))
    pdbscores = []
    for qscore, alnres in zip(quality, sequence.seq):
        if alnres == '-': continue  # Closest sequence not in that alignment column
        pdbres, seqres = pdbseqIter.next(), seqIter.next()
        if seqres == '-': pdbscores.append(0.0)  # Closest sequence has gap compared to pdb sequence
        elif pdbres != '-': pdbscores.append(qscore)  # Pdb sequence has missing residue
    return pdbchain, pdbscores

# # # # #  Command Line  # # # # #
__help__ = """
Usage: python PDB.py <command> [args]

Available commands:
    cleanpdb [--noH] old new         : Remove all lines from 'old' that do not begin
                                       with ATOM or TER, save as 'new'. If --noH is
                                       given, hydrogen atoms will be removed too.
    changechain oldID newID old new  : Opens 'old', changes the chain letter from
                                       'oldID' to 'newID', then saves to 'new'.
    renumber chainID diff pdb new    : Opens 'pdb', modifying the residue numbers for
                                       'chainID' by adding the value of 'diff', and
                                       then saving to 'new'.
    pdbtofasta pdb new               : Extracts the amino acid sequences from 'pdb',
                                       saving them in fasta format as 'new'.
"""

if __name__ == '__main__':
    args = sys.argv[1:]
    if not args or '--help' in args:
        print __help__
        exit()
    command = args[0]
    if command == 'cleanpdb' and 3 <= len(args) <= 4:
        if args[1] == '--noH': noH = True
        else: noH = False
        oldFile, newFile = args[-2:]
        cleanpdb(oldFile, newFile, noH)
    elif command == 'changechain' and len(args) == 5:
        oldFile, newFile = args[3:5]
        oldID, newID = args[1:3]
        changechainID(oldFile, newFile, oldID, newID)
    elif command == 'renumber' and len(args) == 5:
        oldFile, newFile = args[3:5]
        chainID, diff = args[1:3]
        renumber(oldFile, newFile, chainID, int(diff))
    elif command == 'pdbtofasta' and len(args) == 3:
        oldFile, newFile = args[1:3]
        pdbtofasta(oldFile, newFile)
    else:
        print __help__
        exit()
