#!/usr/bin/python
# Author: Dave Curran
# Date: July 2010
# # # # # # # # # # #
"""A module allowing a directory of pdb files to be filtered using distance constraints.

 Usage:   Instanciate. Call setDirectory(dirPath) with the folder containing all of the
            decoys to filter, as well as the scores.fasc file. Call setBasenames((name1, name2)),
            giving the names of the pdb from which the decoys were generated. These are not
            paths, as the files are not accessed, but just strings to identify the
            decoys. Use addDistConstraint(residue, toProtein, distance) to add as many
            constraints as desired. residue should be a string in the form 'PHE171B',
            containing the residue code, its number in the pdb, and the 1-letter chain
            identifying that peptide in the pdb. toProtein is the 1-letter chain
            identifying the peptide to be measured to, and distance is a string, int,
            or float describing the desired minimum distance (in Angstroms) between
            'residue' and the 'toChain' peptide. Call parseDecoys() when ready to filter.
 Objects: self.proteinAtoms - {'B': {'LYS89': [(x1,y1,z1), (x2,y2,z2)...]...}...}
            A dict, containing a dict for each chain where the whole protein must be
            measured. In the sub-dict(s) the keys are the codes for every residue, and
            the value is a list with the xyz coords for each atom.
          self.proteinNs - {'B': {'LYS89': (x,y,z)...}...}
            Similar to above, but instead of a list, there is only one tuple with the
            coordinates of the backbone nitrogen for that residue (this is the first
            atom listed for each residue).
          self.residueAtoms - {'PHE171B' : [(x1,y1,z1), (x2,y2,z2)...]...}
            A dict with keys as the residue code, number, and chain for each residue
            that is being measured to a protein. The value for each is a list of the
            xyz coordinates for every atom in the residue.
          self.residueNs - {'PHE171B' : (x,y,z)...}
            Similar to residueAtoms above, but instead of a list the values are a single
            tuple of the coords of the backbone nitrogen (the first listed atom) of that
            residue.
"""
import os, math
import util

class PdbFilter:
    def __init__(self):
        self.workingDir = None
        self.basenames = ()
        self.scoresDict = {}
        self.distConstraints = {}
        self.residuesToParse = {}
        self.proteinAtoms = {}
        self.proteinNs = {}
        self.residueAtoms = {}
        self.residueNs = {}

    def setDirectory(self, path):
        if os.path.isdir(path):
            self.workingDir = path
            return True
        else: return False

    def setBasenames(self, basenames):
        self.basenames = basenames
        return True

    def addDistConstraint(self, res, toChain, dist):
        """The residue code should be a string containing the amino acid code, the
        residue number, and the chainID of its protein. For example lysine 212 in chain
        T would be specified as 'LYS212T'. toChain is the code of the protein this
        residue should be close to, and dist is that distance in angstroms."""
        if type(toChain) is not str or type(res) is not str:
            return False
        dist = float(dist)
        res = res.strip().upper()
        if toChain not in self.distConstraints:
            self.distConstraints[toChain] = [(res, dist)]
        else:
            self.distConstraints[toChain].append((res, dist))
        if res[-1] not in self.residuesToParse:
            self.residuesToParse[res[-1]] = [res[:-1]]
        else:
            self.residuesToParse[res[-1]].append(res[:-1])
        return True

    def getDistConstraints(self):
        return self.distConstraints

    def parseDecoys(self):
        if not self.workingDir or not self.basenames: return False
        self.scoresDict = util.ScoresDict(self.workingDir)
        if not self.scoresDict:
            return False
        for chain in self.distConstraints:
            print('\nResidues constrained to be near chain {}:\n{}\n'.format(chain, '\n'.join('{}, {:.1f} Angstrom'.format(res, dist) for res, dist in self.distConstraints[chain])) )
        print('Starting to filter pdb files in {}.'.format(self.workingDir))
        deleted = 0
        filesToKeep = []
        basenames = tuple(self.basenames)
        for filename in os.listdir(self.workingDir):
            if filename.startswith(basenames) and filename[:-9] in basenames:
                filepath = os.path.join(self.workingDir, filename)
                if not self.__checkPdbConstraints(filepath):
                    os.remove(filepath)
                    self.scoresDict.remove(filename)
                    deleted += 1
                else:
                    filesToKeep.append(filename)
        self.scoresDict.saveScores()
        print('Filtering completed and score file updated.')
        print('{} files deleted, {} kept, out of {} total pdbs.\n'.format(deleted, len(filesToKeep), deleted+len(filesToKeep)))
        return len(filesToKeep)

    def orderDecoys(self):
        """Calls orderAndFilterDecoys on all files in the working directory."""
        return self.orderAndFilterDecoys(0)

    def orderAndFilterDecoys(self, numToKeep):
        """Takes an integer argument. Examines all of the pdb files in the working directory,
        sorts them by total_score, and deletes all but the top numToKeep. The data of the
        decoys deleted is not lost, but is changed within the score file. The name is changed
        to inputPdb_deleted_0001 with no file extension. This data is kept for graphing and
        other purposes."""
        if not self.workingDir or type(numToKeep) != int or numToKeep < 0:
            return False
        self.scoresDict = util.ScoresDict(self.workingDir)
        if not self.scoresDict:
            return False
        dirList = os.listdir(self.workingDir)
        fileList = sorted([filename for filename in dirList if filename in self.scoresDict],
                          key=lambda f: self.scoresDict[f]['total_score'])
        if numToKeep:
            deleteList = fileList[numToKeep:]
            fileList = fileList[:numToKeep]
            for filename in deleteList:
                os.remove(os.path.join(self.workingDir, filename))
                self.scoresDict.markDeleted(filename)
        self.__orderDecoys(fileList)
        return True

    # # # # # # # # # # # # # # # #  Private Functions  # # # # # # # # # # # # # # # #
    def __checkPdbConstraints(self, filepath):
        self.proteinAtoms, self.proteinNs = {}, {}
        self.residueAtoms, self.residueNs = {}, {}
        if not self.__parseCoordsFromPdb(filepath): # fills out the above 4 dicts
            return False
        return self.__testConstraints()

    def __orderDecoys(self, filesToKeep):
        l = filesToKeep[:]
        basename, _, s = l[0].rpartition('_')
        numDigits = len(s.partition('.')[0])
        for i, oldFilename in enumerate(l):
            l[i] = False
            filename = '%s_%0*d.pdb' % (basename, numDigits, i+1)
            if oldFilename == filename: continue
            if filename in l: # Need to rename existing file
                l[l.index(filename)] = oldFilename
                self.__swap(oldFilename, filename)
            else:
                self.__rename(oldFilename, filename)
        self.scoresDict.saveScores()

    def __parseCoordsFromPdb(self, filepath):
        dc = self.distConstraints; rtp = self.residuesToParse
        pa, pn = self.proteinAtoms, self.proteinNs
        ra, rn = self.residueAtoms, self.residueNs
        try:
            f = open(filepath, 'rb')
            for line in f:
                if not line.startswith('ATOM') or line[13] == 'H': continue
                chainID = line[21]
                if chainID in dc or chainID in rtp:
#                    resCodeNum ='%s%s'%(line[17:20],line[22:26].strip()) # Ex: 'LYS89'
                    resCodeNum =''.join((line[17:20],line[22:26].strip()))
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                if chainID in dc: # It's the protein to measure to.
                    if not pa.setdefault(chainID,{}).setdefault(
                        resCodeNum,[]):  # This is the first atom of this residue
                        pn.setdefault(chainID,{})[resCodeNum] = (x,y,z)
                    pa[chainID][resCodeNum].append((x,y,z))
                if chainID in rtp: # The input residues
                    resCodeNum = '%s%s' % (resCodeNum, chainID) # Ex: 'LYS89B'
                    if not ra.setdefault(resCodeNum,[]):
                        rn[resCodeNum] = (x,y,z)
                    ra[resCodeNum].append((x,y,z))
        finally:
            f.close()
        return True
    def __testConstraints(self):
        """For each constraint, starts iterating through every residue in the other
        protein, comparing the backbone N with that of the residue specified in the
        constraint. If the backbone atoms are close, each atom in both residues is
        compared. If one pair of atoms is found within the specified distance, the
        algorithm moves on. Comparisons are done by comparing X values, then Y, then
        Z; only if these are all close will the full distance calculation be done.
        Returns True if all constraints satisfied, else False"""
        for chainID, l in self.distConstraints.items():
            for res, dist in l:
                nearDist = max(20.0, dist+10.0)
                resNCoords = self.residueNs[res]
                res1Atoms = self.residueAtoms[res]
                for resCode, atomNCoords in self.proteinNs[chainID].items():
                    if abs(resNCoords[0] - atomNCoords[0]) <= nearDist and abs(resNCoords[1] - atomNCoords[1]) <= nearDist and abs(resNCoords[2] - atomNCoords[2]) <= nearDist:
                        res2Atoms = self.proteinAtoms[chainID][resCode]
                        if self.__checkAtomDistances(dist, res1Atoms, res2Atoms):
                            break
                else:             # Only entered if all atoms in 'other' protein iterated
                    return False  # through without satisfying a constraint
        return True
    def __checkAtomDistances(self, dist, coordsList1, coordsList2):
        for (x1,y1,z1), (x2,y2,z2) in ((coords1, coords2) for coords1 in coordsList1 for coords2 in coordsList2):
            if math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2) <= dist:
                return True

    def __rename(self, old, new):
        if old != new:
            oldPath = os.path.join(self.workingDir, old)
            if os.path.exists(oldPath):
                os.rename(oldPath, os.path.join(self.workingDir, new))
            self.scoresDict.rename(old, new)
    def __swap(self, nameA, nameB):
        if nameA == nameB: return
        pathA = os.path.join(self.workingDir, nameA)
        pathB = os.path.join(self.workingDir, nameB)
        if os.path.isfile(pathB):
            os.rename(pathA, pathA+'.temp')
            os.rename(pathB, pathA)
            os.rename(pathA+'.temp', pathB)
        else:
            os.rename(pathA, pathB)
        self.scoresDict.swap(nameA, nameB)
    def __renameAndUpdateFiles(self, filesToKeep):
        # Depreciated in favour of orderDecoys.
        filesToKeep.sort()
        basename, _, s = filesToKeep[0].rpartition('_')
        numDigits = len(s.partition('.')[0])
        for i in xrange(1, len(filesToKeep)+1):
            filename = '%s_%0*d.pdb' % (basename, numDigits, i)
            if filename not in filesToKeep:
                oldFilename = filesToKeep.pop()
                self.__rename(oldFilename, filename)
        self.scoresDict.saveScores()
