
import os, subprocess
import tkMessageBox

def popupInfo(title, message):
    tkMessageBox.showinfo(title=title, message=message)
def popupError(title, message):
    tkMessageBox.showerror(title=title, icon=tkMessageBox.WARNING, message=message)
def startFiles(*filepaths):
    """Takes a sequence of filepath strings."""
    if not filepaths: return
    subprocess.call(('open',) + filepaths)
def residuesFromPdb(filepath):
    if not filepath or not os.path.exists(filepath): return False
    chains, residues = [], {}
    chain, residue, newResidue = '', '', ''
    f = open(filepath, 'rb')
    for line in f:
        if line.startswith('ATOM'):
            chain = line[21]
            if not chain in chains: chains.append(chain)
            newResidue = line[17:20]+line[22:26].strip()
            if newResidue != residue:
                residue = newResidue
                residues.setdefault(chain, []).append(residue)
    f.close()
    return [chains]+[residues[chain] for chain in chains]


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
            print '\nError occured attempting to parse %s.\n'%self.filepath
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


rosedockMainHelpMessage = '    Very important note: This program has been tested with PyMol for Mac version 1.3, and the rosetta suite versions 3.1 and 3.3. Other versions of rosetta may or may not work, as with other versions of PyMol. It is however known that PyMol version 0.99 DOES NOT WORK. Once rosetta is installed make sure that RoseDock knows where it is located. If it has not found it automatically, you can set the paths in the preferences menu.\n' + \
'    This program was written to run Rosetta protein docking, and to be able to parse the multitudes of results. To start a new run, create one pdb file containing the two proteins to be docked together. If possible, remove sections of the structures that you are confident are not involved in the docking, as this can save hours of run-time. The two proteins should be separated from each other by enough space to allow each to rotate. If you have a rough idea of how the proteins interact, rotate each so the appropriate face is towards its partner.\n' + \
'    The number of decoys depends on the type of run; for a global run (where both structures are randomized and no constraints are defined) it is recommended to have 10,000-100,000 decoys, while for a local run (where the proteins have already been roughly positioned and/or constraints have been entered) 1,000-5,000 decoys are recommended. This modelling can take a very long time; by way of example 3,000 decoys of 2 proteins 300 amino acids each generated by 7 cores at 2.8GHz takes approximately 4-5 hours, and up to 8-16 hours with constraints defined.\n' + \
'    For best performance the number of parallel processes should be equal to or one less than the number of processors in your computer. Randomize Structures should be checked if you want to do a gloabl docking run, in which case you should be generating 10,000 decoys or more. Leave it unchecked for a local docking run. The output folder specifies where to save all of the decoys, their scores, and the docking options. The number of top decoys to keep allows you to delete the low scoring decoys after the run is finished, as thousands of them can quickly consume a lot of hard-drive space. When they are deleted the score file will not be altered, so their information will show up in the results display.\n' + \
'    Docking constraints allow you to incorporate previous information about an amino acid into a docking run. If you have mutagenesis information, or can infer importance via homology, you can specify that certain amino acids must be within a certain distance (in Angstroms) of the other protein. The Edit button becomes active after choosing an input pdb, and you can specify amino acids from one or both proteins in the pdb. Remember to exit with the Done button once finished, so that the constraints are saved. After generating the previously specified number of decoys, all of those that do not meet every constraint are deleted, and docking is begun again in order to again reach that number. This continues until every decoy satisfies every constraint. Note: If the constraints are too restrictive this cycle could go on forever, or at least for a very very long time.\n' + \
'    Once options (including constraints) have been filled out, they can be saved for future reference. Likewise, options can be loaded from a .dockfile. The options for a run are automatically saved into the output folder when a run is begun, and the results will be opened once it has completed.\n' + \
'    When examining the results of a docking run, several scores are presented. Total Score is a combination of several scores for the complex, the lower the better. The interface Score is a separate measure, and describes the difference in score between each individual protein and that of their complex. Again the lower the better, with good models usually from -5 to -10. RMS is a measure of the difference of that decoy to the starting structure; its value has little meaning, but the greater the difference in RMS between 2 decoys, the more different they are from each other. The results can be sorted, and double-clicking a decoy will open it in your default pdb viewing program.\n' + \
'    The Tabulate Results button brings up a text box with the RMS and sum of the Total and Interface scores, separated by a tab. Pressing select all and then command-c will copy the data so you can paste it into your favourite graphing or spreadsheet program. This allows the creation of graphs overviewing the docking run, so funnel clusters can be shown.\n' + \
'    Since each decoy is generated independantly, many of them may be very similar. Clustering the results attempts to separate the decoys into clusters of similar structures, allowing a better understanding of the results. This process is computationally intensive, and increases semi-exponentially the more decoys that are considered. 30-200 is a reasonable number, taking a minute or two at most. The cluster radius allows you to define how similar structures must be (in Angstrom) to be considered a cluster.\n'
pathsHelpMessage = '    These paths allow RoseDock to find the Rosetta programs it requires to run.\n    Rosetta_Bundles should be the path to the main directory, and rosetta_database should be just within that folder.\n    The docking executable is usually in rosetta3.x_Bundles/rosetta_source/bin/docking_protocol.your_particular_version. If you cannot pick that alias, the true file should be found in rosetta_Bundles/rosetta_source/build/src.../docking_protocol.your_particular_version.\n\nRemember to press "Save Paths" when you are done.'
