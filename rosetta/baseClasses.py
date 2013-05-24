"""The main classes used to run the various rosetta programs.

Each class is a subclass of dict, and their run options are accessed
directly in this way (example d = rosetta.Docker(); d['-in:file:s']).
If being set manually, the values must be strings.
The major options for each class should have their own set methods, and
these should be preferred over direct access to the options.

Defines:
  load_options -- A function to load a previously saved run. The state
of the object should be the same as just before the run started.
  BaseClass -- Main class that all other program objects inherit
from. Set up so that different programs may overwrite different
functions, allowing functionality with a minimum of overwriting.
  SeedInputPdbs -- A class acting as a patch over BaseClass. This
should be inherited from for those programs that do not respect
multiple processes writing to the same directory, with the same
base name. Beyond initiation, nothing further is required from
the child class.
  FilterPdbs -- A class acting as a patch over BaseClass. This
should be inherited from for programs wishing to use residue
distance constraints. Adds the self.constraints list, and re-
runs the main application untill all constraints in this list
are satisfied for every decoy.
  Prepacker -- A program object, child of BaseClass, that runs the
rosetta prepack algorithm. This is useful for many different runs,
allowing calibration of the score function. Is present here because
it is used by BaseClass itself.
"""
from __future__ import with_statement # Needed for python 2.5
from . import options, PathError # From __init__ file
import util, constants, pdbFilter
import os, subprocess, time, cPickle, datetime
try: import ScoreView
except ImportError:
    ScoreView = None


def load_options(filepath):
    return cPickle.load(open(filepath, 'rb'))

# Implement checks when an option is set. Make sure is str, allow options to
#   be accessed without the '-'; get rid of double ':'? etc.
class BaseClass(util.OrderedDict):
    def __init__(self, execName):
        util.OrderedDict.__init__(self)
        self.bundlePath = options['rosetta_bundle']
        self.dbPath = options['rosetta_database']
        self.execPath = options.paths.get(execName, '')
        self.execName = execName
        self.outputLog = constants.output_log_filename
        self.errorLog = constants.error_log_filename
        self.saveFile = self.execName + constants.saved_run_options_extension
        self.update_items([ ('-in:file:s',''), ('-nstruct','1') ])
        self._finalInputPdb = ''
        self._tempFiles = [] # To be deleted after a run.
        # # #  Common options:
        self.outputName = ''
        self.outputDir = constants.default_outputDir
        self.numStruct = 1
        self.numCPUs = constants.default_numCPUs
        self.runPrepack = False
        self.openScoreViewer = True
        self.renameFinalDecoys = True
        self.keepTopDecoys = constants.default_keepTopDecoys
        # # #  Uncommon options:
        self._runSilent = False
        self._saveState = True
        self._argsInSupplements = []
        
    # # # # #  Common Public Methods  # # # # #
    
    def run(self):
        startTime = time.time()
        try:
            if not self._preRunCheck(): return 3
            self._initialFileSetup()
            if self._saveState: self.save()
            self._preRun()
            args = self._generateArgs()
            supplements = self._generateSupplements()
            retcode = self._run(args, supplements, silent=self._runSilent)
            self._postRun()
        except Exception, excep:
            self._exceptionCallback(excep)
        else:
            self._successPostRun()
        finally:
            self._cleanUp()
        if self.renameFinalDecoys: self._renameDecoys()
        self._reportCompletion(time.time() - startTime)
        if self.openScoreViewer: self._startScoreViewer()
        return retcode
    def save(self, outPath=None):
        if not outPath: outPath = self.outputDir
        if not os.path.isdir(outPath): os.makedirs(outPath)
        with open(os.path.join(outPath, self.saveFile), 'wb') as f:
            cPickle.dump(self, f, 0)
    def cleanPdb(self, inFile, outFile):
        with open(inFile, 'rb') as f:
            buff = [line for line in f if line.startswith('ATOM')]
        with open(outFile, 'wb') as f: f.write(''.join(buff))
    def help(self):
        return self.__doc__

    # # # # #  Common Private Methods  # # # # #

    def _preRunCheck(self):
        if not os.path.isfile(self.inputPdb):
            return False
        return True
    def _initialFileSetup(self):
        outDir = self.outputDir
        if not os.path.isdir(outDir): os.makedirs(outDir)
        outPath = os.path.join(outDir, self.outputName+'.pdb')
        if os.path.isfile(outPath): outPath = os.path.join(
            outDir, self.outputName+'_clean.pdb')
        self.cleanPdb(self.inputPdb, outPath)
        self._finalInputPdb = outPath
        self._tempFiles.append(outPath)
        if self._runPrepack:
            ppk = Prepacker()
            ppk.inputPdb = outPath
            ppk.outputDir = self.outputDir
            ppk._saveState = False
            ppk.run()
        
    def _preRun(self):
        pass
    def _generateArgs(self):
        notThese = self._argsInSupplements[:]
        args = [self.execPath, '-database', self.dbPath]
        if '-in:file:s' not in notThese:
            inpath = self._finalInputPdb
            if ' ' in inpath: inpath = inpath.replace(' ', '\ ')
            args.extend(['-in:file:s', inpath])
            notThese.append('-in:file:s')
        for opt, val in self.items():
            if not val or opt in notThese: continue
            if val.lower() == 'true': args.append(opt)
            else: args.extend([opt] + val.split(' '))
        return args
    def _generateSupplements(self):
        return [[]] * self.numCPUs
    def _postRun(self):
        pass
    def _exceptionCallback(self, excep):
        raise excep
    def _successPostRun(self):
        pass
    def _cleanUp(self):
        for fpath in self._tempFiles[:]:
            if os.path.exists(fpath):
                os.remove(fpath)
                self._tempFiles.remove(fpath)
    def _renameDecoys(self):
        sd = util.ScoresDict(self.outputDir)
        if not sd: return False # Error here?
        sd.orderAndTrimDecoys(self.keepTopDecoys, self.outputName)
        
    def _startScoreViewer(self):
        if ScoreView:
            ScoreView.start(self.outputDir)
        elif not options['scoreView_path']:
            print 'Could not start ScoreView, as the application could not be located. Please update the path in the options.'
        else:
            print 'Starting ScoreView...'
            util.startWithArgs(options['scoreView_path'], [self.outputDir])
    def _run(self, args, supplements=[], numCPUs=None, silent=False):
        """Run the process described by the args list.

        --args should be a list formatted for subprocess.Popen().
        --supplements should be a list of lists, one for each process desired to run in
        parallel. One sublist will be appended to the args list for one subprocess call.
        --numCPUs should be an int, if desired to be different from the self attribute.
        --silent mode doesn't use _countCompleted or _reportProgress
        """
        # is the numCPUs arg needed here?
        curDir = os.getcwd()
        outPath = self.outputDir
        procs, retcode = [], 0
        if not numCPUs: numCPUs = self.numCPUs
        while len(supplements) < numCPUs:
            supplements.append([])
        os.chdir(outPath)
        with open(self.outputLog, 'ab') as outlog:
            with open(self.errorLog, 'ab') as errlog:
                for n, supp in zip(range(numCPUs), supplements):
                    argsList = args + supp
                    print 'Beginning %s process #%i with arguments:\n%s\n' % (
                        self.execName, n, ' '.join(argsList))
                    procs.append(subprocess.Popen(
                        argsList, stdout=outlog.fileno(),
                        stderr=errlog.fileno() ) )
                    time.sleep(2)
                if silent:
                    for p in procs: p.wait()
                else:
                    countInterval = constants.count_completed_interval
                    prevCount, total = -1, int(self['-nstruct'])
                    while None in (p.poll() for p in procs):
                        count = self._countCompleted()
                        if count != prevCount:
                            self._reportProgress(count, total)
                            prevCount = count
                        time.sleep(countInterval)
                    self._reportProgress(self._countCompleted(), total)
        if any(p.returncode for p in procs):
            retcode = 1
        if os.path.isfile(self.errorLog) and os.path.getsize(self.errorLog) == 0:
            os.remove(self.errorLog)
        os.chdir(curDir)
        return retcode

    def _getOutputBasenames(self):
        name = self.outputName or self['-in:file:s']
        if name.endswith('.pdb'): name = name[:-4]
        self.outputName = name
        return [name]

    def _countCompleted(self):
        # if slow can cache listdir contents
        # relies on the fact that output structures end with _0002.pdb. 4 digit number.
        count = 0
        basenames = tuple(self._getOutputBasenames())
        for fname in os.listdir(self.outputDir):
            if fname.startswith(basenames) and fname[:-9] in basenames:
                count += 1
        return count
            
    def _reportProgress(self, numDone, numTotal):
        now = datetime.datetime.now()
        timestamp = str(datetime.time(now.hour, now.minute, now.second))
        print '\t%s - %i of %i results completed.' % (timestamp, numDone, numTotal)
    def _reportCompletion(self, runTime):
        mins, secs = divmod(runTime, 60)
        hours, mins = divmod(mins, 60)
        timestr = ('%i hours, ' % hours if hours else '') +\
                  ('%i minutes, ' % mins if mins else '') +\
                  ('%.1f seconds' % secs)
        s = '%s run completed in %s.' % (self.execName, timestr)
        dashes = '-'*len(s)
        print '%s\n%s\n%s\n' % (dashes, s, dashes)
    
    # # # # #  Under-the-Hood Methods  # # # # #

    def _getInputPdb(self): return self['-in:file:s']
    def _setInputPdb(self, path):
        if not os.path.isfile(path):
            print '%s could not be found to set inputPdb' % path
            return
        self['-in:file:s'] = os.path.realpath(path)
    def _getOutputName(self):
        if not self._outputName:
            self.outputName = os.path.basename(self['-in:file:s'])
        return self._outputName
    def _setOutputName(self, name):
        if name.endswith('.pdb'): name = name[:-4]
        if ' ' in name:
            name = name.replace(' ', '_')
            print "Due to Rosetta's option handling, there can be no spaces in the file names. Any spaces have been converted to underscores."
        self._outputName = name
    def _getOutputDir(self): return self._outputDir
    def _setOutputDir(self, path):
        if not os.path.isabs(path):
            path = os.path.join(os.getcwd(), path)
        path = os.path.realpath(path)
        if ' ' in path: raise PathError("Due to Rosetta's option handling, the output path cannot contain any spaces.")
        self._outputDir = path
    def _getNumStruct(self): return int(self['-nstruct'])
    def _setNumStruct(self, num):
        num = int(num)
        self['-nstruct'] = str(num)
    def _getNumCPUs(self): return self._numCPUs
    def _setNumCPUs(self, val): self._numCPUs = int(val)
    def _getNumTopDecoys(self): return self._keepTopDecoys
    def _setNumTopDecoys(self, val): self._keepTopDecoys = int(val)
    def _getRunPrepack(self): return self._runPrepack
    def _setRunPrepack(self, val): self._runPrepack = val
    inputPdb = property(_getInputPdb, _setInputPdb)
    outputName = property(_getOutputName, _setOutputName)
    outputDir = property(_getOutputDir, _setOutputDir)
    numStruct = property(_getNumStruct, _setNumStruct)
    numCPUs = property(_getNumCPUs, _setNumCPUs)
    keepTopDecoys = property(_getNumTopDecoys, _setNumTopDecoys)
    runPrepack = property(_getRunPrepack, _setRunPrepack)

    def __repr__(self):
        return 'rosetta.%s object at %s' % (self.__class__.__name__,
                                            hex(id(self)))
    def __str__(self):
        opts = '{%s}' % ', '.join('%r: %r'%(k,v) for k,v in self.items())
        return "Wrapper for rosetta's %s. Options:\n%s" % (self.execName, opts)

class SeedInputPdbs:
    """A patch for BaseClass where the input pdb is seeded.

    The 1 input pdb is copied, one per process. This is for the many
    applications that would otherwise overwrite output files when run
    with multiple processes. The class is used in the exact same way."""
    def __init__(self):
        self._argsInSupplements = ['-in:file:s', '-nstruct']
    def _initialFileSetup(self):
        BaseClass._initialFileSetup(self)
        outDir = self.outputDir
        cleanPath = self._finalInputPdb
        data = open(cleanPath, 'rb').read()
        for name in self._getOutputBasenames():
            path = os.path.join(outDir, name + '.pdb')
            with open(path, 'wb') as f: f.write(data)
            self._tempFiles.append(path)
        self._tempFiles.append(cleanPath)
    def _generateSupplements(self):
        basenames = self._getOutputBasenames()
        numStruct, numCPUs = int(self['-nstruct']), float(self.numCPUs)
        supps = []
        for name in basenames:
            num = round(numStruct / numCPUs)
            numStruct -= num
            numCPUs -= 1
            supp = ['-in:file:s', name+'.pdb', '-nstruct', str(int(num))]
            supps.append(supp)
        return supps
    def _getOutputBasenames(self):
        name = self.outputName or self['-in:file:s']
        if name.endswith('.pdb'): name = name[:-4]
        self.outputName = name
        if self.numCPUs > 1:
            names = [name+'_%i'%i for i in range(self.numCPUs)]
        else: names = [name]
        return names


class FilterPdbs:
    """A patch for BaseClass allowing residue distance constraints.

    -- Adds the constraints variable.

    This will run the master application again and again untill all specified
    distance constraints are satisfied in every generated decoy. To use this,
    fill out the self.constraints list with as many constraints as desired.
    Each should consist of (resID, toChain, dist). resID is the residue code,
    the number for that residue, and the chain of that peptide ex: 'PHE171B'.
    toChain is the chainID of the other protein you want your residue to be
    close to, ex: 'B', and dist is an integer specifying the distance in
    Angstrom it must be within. So to ensure F171 and E112 from chain 'B' are
    within 10 A of protein 'T', self.constraints would look like:
    [('PHE171B', 'T', 10), ('GLU112B', 'T', 10)]."""
    def __init__(self):
        self.constraints = []
        self._pdbFilter = None
    def run(self):
        startTime = time.time()
        retcode = 0
        try:
            if not self._preRunCheck(): return 3
            self._initialFileSetup()
            if self._saveState: self.save()
            self._preRun()
            args = self._generateArgs()
            supplements = self._generateSupplements()
            if not self.constraints:
                retcode = self._run(args, supplements, silent=self._runSilent)
            else:
                self._setupPdbFilter()
                numFinished, numToMake = 0, self.numStruct
                while numFinished < numToMake:
                    retcode = self._run(args, supplements, silent=self._runSilent)
                    numFinished = self._pdbFilter.parseDecoys()
                    if numFinished is False: return 2
            self._postRun()
        except Exception, excep:
            self._exceptionCallback(excep)
        else:
            self._successPostRun()
        finally:
            self._cleanUp()
        if self.renameFinalDecoys: self._renameDecoys()
        self._reportCompletion(time.time() - startTime)
        if self.openScoreViewer: self._startScoreViewer()
        return retcode
    def _setupPdbFilter(self):
        self._pdbFilter = pdbFilter.PdbFilter()
        self._pdbFilter.setDirectory(self.outputDir)
        self._pdbFilter.setBasenames(self._getOutputBasenames())
        for res, toChain, dist in self.constraints:
            self._pdbFilter.addDistConstraint(res, toChain, dist)


class Prepacker(BaseClass):
    """
    Class to run rosetta's prepack algorithm.
    """
    def __init__(self):
        BaseClass.__init__(self, 'docking_prepack_protocol')
        self.renameFinalDecoys = False
        self.openScoreViewer = False
        self.numCPUs = 1
        self._runSilent = True
        if not self.execPath:
            self._version = 3.1
            self.execPath = options.paths['docking_protocol']
            self.update_items(
                [('-in:file:s',''), ('-nstruct','1'), ('-docking:dock_ppk','true'),
                 ('-run:no_scorefile','true'), ('-partners','')] )
            # don't know if no_scorefile is working here. if not, implement.
        else:
            self._version = 3.3
            self.update_items(
                [('-in:file:s',''), ('-nstruct','1'),
                 ('-run:no_scorefile','true'), ('-partners','')] )
    def _initialFileSetup(self):
        outDir = self.outputDir
        if not os.path.isdir(outDir): os.makedirs(outDir)
        outPath = os.path.join(outDir, self.outputName+'.pdb')
        if os.path.isfile(os.path.join(outDir, self.outputName+'_0001.pdb')):
            outPath = os.path.join(outDir, self.outputName+'_clean.pdb')
            self._tempFiles.append(outPath)
        self.cleanPdb(self.inputPdb, outPath)
        self._finalInputPdb = outPath
    def _cleanUp(self):
        cleanPath = self._finalInputPdb[:-4] + '_0001.pdb'
        cleanName = os.path.basename(cleanPath)
        if self._version == 3.3: #or higher? 3.2?
            for pref in ('away_', 'away_packed_', 'initial_', 'prepack_'):
                path = os.path.join(self.outputDir, pref+cleanName)
                self._tempFiles.append(path)
        else:
            pass # anything to do here?
        os.rename(cleanPath, os.path.join(self.outputDir, self.outputName+'.pdb'))
        BaseClass._cleanUp(self)
