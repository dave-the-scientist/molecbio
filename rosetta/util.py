# Author: Dave Curran
# Date:   April 2012
"""Utility classes used in the rosetta package.

Contains:
-- class OrderedDict
-- class ScoresDict

For the ScoresDict object, the different types of score files are defined by
the score file suffix as well as the measurements it contains.
-- docking: 'total_score' and 'I_sc' in measurements.
-- homology: 'total_score' and 'aln_len' in measurements.
-- floppytail:  'total_score' in measurements, suffix is '.sc'.
-- loopmodel: 'total_energy' in measurements.
-- abinitio: 'score' in measurements, suffix is '.fsc'.
"""
import os, sys, subprocess

def startFiles(*filepaths):
    """Takes a sequence of filepath strings."""
    ### Test with linux.
    if not filepaths: return
    subprocess.call(('open',) + filepaths)
def startWithArgs(filepath, args=[]):
    # test with linux.
    def isExecutable(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    if not isExecutable(filepath) and filepath.endswith('.app'):
        progName = os.path.basename(filepath)[:-4]
        exePath = os.path.join(filepath, 'Contents', 'MacOS', progName)
        if isExecutable(exePath): filepath = exePath
    subprocess.Popen([filepath] + list(args))

def  determineNumCPUs():
    """Taken from:
    http://stackoverflow.com/questions/1006289/how-to-find-out-the-number-of-cpus-in-python"""
    # Python 2.6+
    try:
        import multiprocessing
        return multiprocessing.cpu_count()
    except (ImportError,NotImplementedError): pass
    # POSIX
    try:
        res = int(os.sysconf('SC_NPROCESSORS_ONLN'))
        if res > 0: return res
    except (AttributeError,ValueError): pass
    # Windows
    try:
        res = int(os.environ['NUMBER_OF_PROCESSORS'])
        if res > 0: return res
    except (KeyError, ValueError): pass
    # jython
    try:
        from java.lang import Runtime
        runtime = Runtime.getRuntime()
        res = runtime.availableProcessors()
        if res > 0: return res
    except ImportError: pass
    # BSD
    try:
        sysctl = subprocess.Popen(['sysctl', '-n', 'hw.ncpu'],
                                      stdout=subprocess.PIPE)
        scStdout = sysctl.communicate()[0]
        res = int(scStdout)
        if res > 0: return res
    except (OSError, ValueError): pass
    # Linux
    try:
        res = open('/proc/cpuinfo').read().count('processor\t:')
        if res > 0: return res
    except IOError: pass
    # Solaris
    try:
        import re
        pseudoDevices = os.listdir('/devices/pseudo/')
        expr = re.compile('^cpuid@[0-9]+$')
        res = 0
        for pd in pseudoDevices:
            if expr.match(pd) != None:
                res += 1
        if res > 0: return res
    except OSError: pass
    # Other UNIXes (heuristic)
    try:
        try:
            dmesg = open('/var/run/dmesg.boot').read()
        except IOError:
            dmesgProcess = subprocess.Popen(['dmesg'], stdout=subprocess.PIPE)
            dmesg = dmesgProcess.communicate()[0]
        res = 0
        while '\ncpu' + str(res) + ':' in dmesg:
            res += 1
        if res > 0: return res
    except OSError: pass
    return 1


class OrderedDict(dict):
    def __init__(self, initVal=None):
        dict.__init__(self)
        self.__keys = []
        if initVal:
            try:
                if hasattr(initVal, 'values'): self.update(initVal)
                else: self.update_items(initVal)
            except AttributeError:
                raise AttributeError('OrderedDict expected a dictionary or sequence of (key,value) pairs as its argument.')
    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        if key not in self.__keys: self.__keys.append(key)
    def __delitem__(self, key):
        dict.__delitem__(self, key)
        self.__keys.remove(key)
    def __repr__(self):
        return '{%s}' % ', '.join('%r: %r' % item for item in self.iteritems())
    def __iter__(self, reverse=False):
        if reverse: return reversed(self.__keys)
        else: return iter(self.__keys)
    iterkeys = __iter__
    def itervalues(self, reverse=False):
        return (self[key] for key in self.iterkeys(reverse))
    def iteritems(self, reverse=False):
        return ((key, self[key]) for key in self.iterkeys(reverse))
    def keys(self, reverse=False): return list(self.iterkeys(reverse))
    def values(self, reverse=False): return list(self.itervalues(reverse))
    def items(self, reverse=False): return list(self.iteritems(reverse))
    def update(self, d): map(self.__setitem__, d, d.values())
    def update_items(self, seq): map(lambda pair: self.__setitem__(*pair), seq)
    def copy(self): return OrderedDict(self)
    def clear(self): map(self.__delitem__, self.keys())
    def sort(self, rev=False): self.__keys.sort(reverse=rev)
    def reverse(self): self.__keys.reverse()
    def popitem(self):
        k,v = dict.popitem(self)
        self.__delitem__(k)
        return k,v
    def pop(self, key):
        v = self[key]
        self.__delitem__(key)
        return v

class ScoresDict(dict):
    """ The scoreType() method indicates what type of score file has been opened, if
    known. Returns a string, one of 'docking', 'abinitio', 'homology', or
    'floppytail'. metrics()
    returns as a list of strings the scoring metrics in this score file."""
    def __init__(self, filepath):
        dict.__init__(self)
        self.filepath = self.__resolveFilepath(filepath)
        if not self.filepath: return
        self.dirpath = os.path.dirname(self.filepath)
        self.lines = {}
        self.__header = ''
        self.__metrics = []
        self.__rankingMetric = None
        self.__type = None
        self.__parseScoresFile()
        self.__determineType()
    def remove(self, name):
        if name in self:
            del self[name]
            del self.lines[name]
    def saveScores(self, filepath=None):
        if filepath == None: filepath = self.filepath
        buff = [self.__header]
        for name in sorted(self):
            buff.append(self.lines[name])
        buff.append('') # So the file can be appended to.
        f = open(filepath, 'wb')
        f.write('\n'.join(buff)); f.close()

    def sort(self, metric=None):
        if metric not in self.__metrics:
            metric = self.__rankingMetric
        return sorted(self, key=lambda name: self[name][metric])
    def score(self, filename):
        return self[filename][self.__rankingMetric]
    def scoreType(self):
        return self.__type
    def metrics(self):
        return self.__metrics

    def markDeleted(self, name, basename=None):
        if not basename: basename = name[:-9]
        nameTemplate = '%s_deleted_%04d.pdb'
        i = 1
        while nameTemplate % (basename, i) in self: i += 1
        self.rename(name, nameTemplate % (basename, i))
    def rename(self, old, new):
        if old == new: return
        oldScoreLine = self.lines.pop(old)
        newScoreLine = '%s %s' % (oldScoreLine.rpartition(' ')[0], new)
        self.lines[new] = newScoreLine
        self[new] = self.pop(old)
    def swap(self, nameA, nameB):
        if nameA == nameB: return
        tempName = nameA+'.temp'
        self.rename(nameA, tempName)
        self.rename(nameB, nameA)
        self.rename(tempName, nameB)

    # # # # #  File handling methods  # # # # #
    def orderDecoys(self, basename=None, rankingMetric=None):
        self.orderAndTrimDecoys(0, basename)
    def orderAndTrimDecoys(self, numToKeep, basename=None):
        files = os.listdir(os.path.dirname(self.filepath))
        newOrder = [fname for fname in self.sort() if fname in files]
        if numToKeep:
            toDelete = newOrder[numToKeep:]
            newOrder = newOrder[:numToKeep]
            for fname in toDelete:
                os.remove(os.path.join(self.dirpath, fname))
                self.markDeleted(fname, basename)
        self.__reorderDecoys(newOrder, basename)

    def renameFiles(self, old, new):
        if old == new: return
        self.rename(old, new)
        oldPath = os.path.join(self.dirpath, old)
        newPath = os.path.join(self.dirpath, new)
        os.rename(oldPath, newPath)
    def swapFiles(self, nameA, nameB):
        if nameA == nameB: return
        tempName = nameA+'.temp'
        self.renameFiles(nameA, tempName)
        self.renameFiles(nameB, nameA)
        self.renameFiles(tempName, nameB)

    # # # # #  Private Methods  # # # # #
    def __reorderDecoys(self, newOrder, basename=None):
        l = newOrder[:]
        if not basename: basename = l[0].rpartition('_')[0]
        for i, oldFilename in enumerate(l):
            l[i] = False
            filename = '%s_%.4d.pdb' % (basename, i+1)
            if oldFilename == filename: continue
            if filename in l: # Need to rename existing file
                l[l.index(filename)] = oldFilename
                self.swapFiles(oldFilename, filename)
            else:
                self.renameFiles(oldFilename, filename)
        self.saveScores()

    def __isScoreFile(self, filename):
        suffixes = ('fasc', 'fsc', 'sc')
        if filename.lower().split('.')[-1] in suffixes:
            return True
        return False
    def __resolveFilepath(self, filepath):
        if type(filepath) is not str: return False
        if os.path.isfile(filepath) and self.__isScoreFile(filepath):
            return os.path.realpath(filepath)
        elif os.path.isdir(filepath):
            for f in os.listdir(filepath):
                if self.__isScoreFile(f):
                    return os.path.realpath(os.path.join(filepath, f))
        return False
    def __parseScoresFile(self):
        f = open(self.filepath, 'r')
        def floatIfIs(score):
            try: return float(score)
            except ValueError: return score
        try:
            temp = f.readline()
            if temp.startswith('SEQUENCE:'): header = f.readline().strip()
            else: header, temp = temp.strip(), ''
            scoresHeader = [score.strip() for score in header.split()[1:]
                            if not score.isspace() and 'description' not in score]
            self.__metrics = scoresHeader
            self.__header = temp+header
            self.__rankingMetric = scoresHeader
            for line in f:
                resultDict = {}
                segs = [seg.strip() for seg in line.split()[1:]
                        if not seg.isspace()]
                name = os.path.basename(segs.pop().strip())
                if not name.endswith('.pdb'): name += '.pdb' # Results must be pdbs.
                resultDict = dict((metric, floatIfIs(score)) for metric, score in
                                  zip(scoresHeader, segs))
                self[name] = resultDict.copy()
                self.lines[name] = line.strip()
            return True
        except:
            print('\nError occured attempting to parse {}.\n'.format(self.filepath))
            self.__header = ''
            self.lines = {}
            self.clear()
            raise
        finally:
            f.close()
    def __determineType(self):
        suffix = self.filepath.rpartition('.')[-1]
        metrics = self.__metrics
        if 'total_score' in metrics:
            self.__rankingMetric = 'total_score'
            if 'I_sc' in metrics:
                self.__type = 'docking'
            elif 'aln_len' in metrics:
                self.__type = 'homology'
            elif suffix == 'sc':
                self.__type = 'floppytail'
        elif 'total_energy' in metrics:
            self.__type = 'loopmodel'
            self.__rankingMetric = 'total_energy'
        elif 'score' in metrics:
            self.__rankingMetric = 'score'
            if suffix == 'fsc':
                self.__type = 'abinitio'
        else:
            self.__type = 'unknown'
            self.__rankingMetric = metrics[0]
