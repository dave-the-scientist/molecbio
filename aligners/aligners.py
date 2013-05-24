"""
"""
import itertools, operator
from molecbio import sequ, blosum
import nwmodule

class Needleman(object):
    """Performs a Needleman-Wunsch global alignment on 2 sequences. Implemented
    in C, so about 100x faster than the python version.
    The align() function ensures the sequences are properly formatted, but this
    can be ~10% slower if they're already formatted, and takes twice as long if
    they're Sequence objects. Faster to format the sequences beforehand, and then
    use the align_raw() function."""
    def __init__(self, match=2, mismatch=-1, gap=-1,
                 score_matrix=blosum.blosum62, blosum_gap_open=-10, blosum_gap=-1):
        self.matchscore = match
        self.mismatchscore = mismatch
        self.gapscore = gap
        self.score_matrix = score_matrix
        self.blosum_gap_open = blosum_gap_open
        self.blosum_gap = blosum_gap

    def align(self, seq1, seq2):
        seq1 = self._filterSeq(seq1)
        seq2 = self._filterSeq(seq2)
        return self._align(seq1, seq2)
    
    def align_raw(self, seq1, seq2):
        return self._align(seq1, seq2)

    def pairwiseAlign(self, seqList):
        """Aligns each sequence to each other in the list, returning a list of
        tuples containing the alignments."""
        return list(self._pairwiseAlignGen(seqList))

    def calcPercentIdentity(self, seq1, seq2):
        """Runs an alignment, returning the percent identity between the sequences."""
        return self._percentIdentity(self.algin(seq1, seq2))

    def calcPairwisePercentIdentity(self, seqList):
        return [self._percentIdentity(seqs) for seqs in self._pairwiseAlignGen(seqList)]

    def oneVsAllIdentity(self, seq1, seqList):
        seq1 = self._filterSeq(seq1)
        seqList = self._filterSeqList(seqList)
        return [self._percentIdentity(self._align(seq1,seq2)) for seq2 in seqList]

    # # # # #  Overwritable Methods  # # # # #
    def _align(self, seq1, seq2):
        """Returns the aligned sequences as a tuple of strings."""
        return nwmodule.align(seq1, seq2, self.matchscore,
                              self.mismatchscore, self.gapscore)

    # # # # #  Class Private Methods  # # # # #
    def _filterSeq(self, seq):
        if type(seq) == sequ.Sequence: seq = seq.seq
        else: seq = str(seq)
        return filter(str.isalpha, seq).upper()
    def _filterSeqList(self, seqList):
        try: seqList = [filter(str.isalpha, s.seq).upper() for s in seqList]
        except: seqList = [filter(str.isalpha, s).upper() for s in seqList]
        return seqList
    def _pairwiseAlignGen(self, seqList):
        seqList = self._filterSeqList(seqList)
        for seq1, seq2 in itertools.combinations(seqList, 2):
            yield self._align(seq1, seq2)
    def _percentIdentity(self, align1, align2):
        matches = map(operator.eq, align1, align2)
        return sum(matches) * 100.0 / len(matches)
    

class PyNeedleman(Needleman):
    """Performs a Needleman-Wunsch global alignment on 2 sequences."""
    def __init__(self, match=2, mismatch=-1, gap=-1,
                 score_matrix=blosum.blosum62, blosum_gap_open=-10, blosum_gap=-1):
        Needleman.__init__(self, match, mismatch, gap, score_matrix,
                           blosum_gap_open, blosum_gap)
        # # #  Private variables
        self._scores = []
        self._paths = []

    # # # # #  Overwritable Methods  # # # # #
    def _align(self, seq1, seq2):
        self._scores = []
        self._paths = []
        self.__generateMatrices(seq1, seq2)
        align1, align2 = self.__backtrack(seq1, seq2)
        
        return align1, align2

    # # # # #  Private Methods  # # # # #
    def __generateMatrices(self, seq1, seq2):
        """For the paths matrix, 1 is diagonal, 2 is left, 3 is up, 0 is end.
        _scores is a n+1 x m+1 matrix, while _paths is n x m."""
        match = self.matchscore; mismatch = self.mismatchscore
        gap = self.gapscore
        prevScores = [gap * i for i in range(1, len(seq1)+1)]
        iter1 = itertools.cycle(seq1)
        
        for j, n in enumerate(seq2):
            diag = gap * j
            left = diag + gap
            scores = []
            paths = []
            for up, m in itertools.izip(prevScores, iter1):
                if m == n: score = match
                else: score = mismatch
                diagscore = diag + score
                leftscore = left + gap
                upscore = up + gap
                score = max(diagscore, upscore, leftscore)
                if score == diagscore: paths.append(1)
                elif score == leftscore: paths.append(2)
                elif score == upscore: paths.append(3)
                scores.append(score)
                diag = up
                left = score
            prevScores = scores
            self._scores.append(scores)
            self._paths.append(paths)

    def __backtrack(self, seq1, seq2):
        paths = self._paths
        align1, align2 = [], []
        iter1 = reversed(seq1); iter2 = reversed(seq2)
        i, j = len(seq1) - 1, len(seq2) - 1
        while i != -1 and j != -1:
            path = paths[j][i]
            if path == 1:
                align1.append(iter1.next())
                align2.append(iter2.next())
                i -= 1; j -= 1
            elif path == 2:
                align1.append(iter1.next())
                align2.append('-')
                i -= 1
            else:
                align1.append('-')
                align2.append(iter2.next())
                j -= 1
        for c1, c2 in itertools.izip_longest(iter1, iter2, fillvalue='-'):
            align1.append(c1)
            align2.append(c2)
        return ''.join(reversed(align1)), ''.join(reversed(align2))
        


