"""Classes to run rosetta applications.

These classes all inherit from baseClasses. This module defines the following:

  Docker -- Runs the rosetta docking program.
  Prepacker -- Runs the prepack program. This is used before running some other
programs, to ensure that the score is properly calibrated to the starting structure.

Notes:
-- If a program is run that uses SeedInputPdbs, and the renameFinalDecoys is set
to True (which is default), and the user wants to generate more decoys, they must
be aware. Say 2000 were created, but the user now wants 2500 total. In this case, the
number of new decoys should be passed as nstruct; so 500. If the program did not use
SeedInputPdbs, or if it did but the renameFinalDecoys was set to False, then nstruct
should be set to 2500.
-- If there are pdb files in the output directory, with names similar to the output
name, that are not recorded in the score file for that directory, they could be
overwritten and lost.
-- If the input pdb file is already present in the output directory, or if the
specified outputName+'.pdb' exists, Prepacker and some other programs may overwrite
it unless the outputName field is set to something else.
-- If one of these classes is used within a GUI, the two relevant functions to
overwrite are _reportProgress(self, numDone, numTotal) and
_reportCompletion(self, runTime). The number of completed structures will be counted
every 5 seconds, and if the number changes the reportProgress function will be called.
Once all of the structures have finished, the reportCompletion function will be called.
"""
from __future__ import with_statement # Needed for python 2.5
from baseClasses import *


class Docker(FilterPdbs, SeedInputPdbs, BaseClass):
    """
    Class to run rosetta's docking protocol.

    Attributes:
    outputName -- Name to be used for output. String.
    outputDir -- Directory to output all files from the run. Path.
    inputPdb -- Starting pdb structure. Path.
    numCPUs -- Number of processors to use. Integer, default 1.
    numStruct -- Number of structures to generate. Integer, default 1.
    randomize -- Randomize initial orientations. Boolean, default False.
    runPrepack -- Run rosettas prepack protocol first. Boolean, default True.
    renameFinalDecoys -- Rename and sort decoys after the run. Boolean, default True.
    keepTopDecoys -- Delete all but top scoring decoys. Integer, default 500.
    constraints -- List of constraints.
    """
    def __init__(self):
        BaseClass.__init__(self, 'docking_protocol')
        SeedInputPdbs.__init__(self)
        FilterPdbs.__init__(self)
        # # #  Default run options:
        scorefile = constants.docking_score_file
        self.update_items(
            [('-in:file:s',''), ('-nstruct','1'), ('-partners',''),
             ('-randomize1',''), ('-randomize2',''), ('-docking:spin','true'),
             ('-dock_pert','3 8'), ('-ex1','true'), ('-ex2aro','true'),
             ('-use_input_sc','true'), ('-out:file:fullatom','true'),
             ('-out:file:o',scorefile), ('-mute','core.util.prof')] )
        # # #  Run variables:
        self.randomize = False
        self.runPrepack = True

    # # # # #  Property Getters and Setters  # # # # #

    def _getRandomize(self): return self._randomize
    def _setRandomize(self, val):
        if val: self['-randomize1'] = self['-randomize2'] = 'true'
        else: self['-randomize1'] = self['-randomize2'] = ''
        self._randomize = val
    randomize = property(_getRandomize, _setRandomize)


class Relaxer(BaseClass):
    """
    Class to run rosetta's relax algorithm.

    Attributes:
    thorough -- Set True to enable thorough run, False to do fast. Default False.
    outputName -- Name to be used for output. String.
    outputDir -- Directory to output all files from the run. Path.
    inputPdb -- Starting pdb structure. Path.
    numCPUs -- Number of processors to use. Integer, default 1.
    numStruct -- Number of structures to generate. Integer, default 1.
    """
    def __init__(self):
        BaseClass.__init__(self, 'relax')
        self.update_items(
            [('-in:file:s',''), ('-nstruct','1'), ('-in:file:fullatom','true'),
             ('-relax:fast','true')] )
        self.thorough = False
        self.openScoreViewer = False

    def _getThorough(self):
        if '-relax:thorough' in self: return True
        return False
    def _setThorough(self, val=True):
        if val: toAdd, toRemove = '-relax:thorough', '-relax:fast'
        else: toAdd, toRemove = '-relax:fast', '-relax:thorough'
        if toRemove in self: self.pop(toRemove)
        self[toAdd] = 'true'
    thorough = property(_getThorough, _setThorough)

class Prepacker(Prepacker):
    """
    Class to run rosetta's prepack algorithm.

    Attributes:
    
    Methods:

    """
    pass # Code is located in baseClasses.py
