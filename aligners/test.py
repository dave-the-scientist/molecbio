import os, sys, time, itertools
os.system('python setup.py build_ext --inplace')
import nwmodule
import aligners
from molecbio import align

seqs = align.open_fasta("/Users/davecurran/Desktop/Master's/Sequence_Files/Ap_Our_TbpB.faa")
seqs = align.open_fasta("/Users/davecurran/Desktop/phd_wip/CYPs/Phylogenetics/Hc13_Ce88_Pp3_aln.fasta")[:10]
seqs = [filter(str.isalpha, s.seq).upper() for s in seqs]
#seq1 = str(seqs[0].seq)
#seq2 = str(seqs[4].seq)
match = 2; mismatch = -1; gap = -2
reportEvery = 5
align1, align2 = [], []

print '\nC version:'
t1 = time.time()
nw_c = aligners.Needleman(match, mismatch, gap)
for i, (seq1, seq2) in enumerate(itertools.combinations(seqs, 2)):
    align2.append(nw_c.align(seq1, seq2))
    if not i % reportEvery:
        sys.stdout.write("\r%i comparisons finished." % i)
        sys.stdout.flush()
print "\n%.4f seconds." % (time.time()-t1)

#exit()

print '\nPython version:'
t = time.time()
nw = aligners.PyNeedleman(match, mismatch, gap)
for i, (seq1, seq2) in enumerate(itertools.combinations(seqs, 2)):
    align1.append(nw.align(seq1, seq2))
    if not i % reportEvery:
        sys.stdout.write("\r%i comparisons finished." % i)
        sys.stdout.flush()
print "\n%.4f seconds." % (time.time()-t)

print "\nAlignments the same: %s." % (align1 == align2)
