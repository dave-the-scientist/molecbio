# Options database for the package. Reference to the dict is stored in the __init__ file
from __future__ import with_statement # Needed for python 2.5
import constants, util
from . import OptionWarning, DirectoryPathWarning, ExecutablePathWarning
import os

__saveDir__ = os.path.realpath(os.path.dirname(__file__))
if 'site-packages.zip' in __saveDir__ and 'Resources' in __saveDir__:
    __saveDir__ = __saveDir__[:__saveDir__.find('Resources') + 9] # If its in an app.
__default_general__ = util.OrderedDict(constants.default_general)
__default_paths__ = util.OrderedDict(constants.default_paths)

class OptionsDict(util.OrderedDict):
    """Options for the molecbio.rosetta package.

    Is itself an OrderedDict object, so can be accessed and modified
    using the standard dictionary interface. The paths to all of the
    rosetta executables are stored in another OrderedDict object, which
    can be accessed at self.paths. Changes will be good for one session
    only unless the save method is used.
    """
    def __init__(self):
        util.OrderedDict.__init__(self)
        self.paths = util.OrderedDict()
        self.options_filepath = os.path.join(
            __saveDir__, constants.options_filename)

        self.load()

    # # #  Public Methods  # # #
    def load(self): # Allow for different path to be specified?
        self.reset()
        if os.path.isfile(self.options_filepath):
            self.__parseOptions()
        self.validate()

    def save(self):
        def optsAsStr(section, opts):
            header = '[%s]' % section
            l = ['%s = %s' % (key, val) for key, val in opts.items()]
            return '\n'.join([header]+l)
        buff = [optsAsStr('general', self), optsAsStr('paths', self.paths)]
        with open(self.options_filepath, 'wb') as f: f.write('\n\n'.join(buff))

    def validate(self):
        self.__checkOptions()
        self.__checkPaths()

    def reset(self):
        self.paths.clear()
        self.paths.update(__default_paths__)
        self.clear()
        self.update(__default_general__)

    # # #  Private Methods  # # #
    def __parseOptions(self):
        with open(self.options_filepath, 'rb') as f:
            l = [['pre-section',[]]]
            for line in f:
                line = line.strip()
                if line.startswith('[') and line.endswith(']'):
                    section = line[1:-1]
                    l.append([section, []])
                elif line:
                    key, _, val = line.partition('=')
                    l[-1][1].append((key.strip(), val.strip()))
        for (section, opts) in l:
            if section == 'general':
                self.update(dict(opts))
            elif section == 'paths':
                self.paths.update(dict(opts))
    def __checkOptions(self):
        if not os.path.isdir(self['rosetta_bundle']):
            self['rosetta_bundle'] = self.__findRosetta()
        if not os.path.isdir(self['rosetta_database']):
            self['rosetta_database'] = self.__findDb()
        if not os.path.isdir(self['executables_dir']):
            self['executables_dir'] = self.__findExecDir()
        if not self['executables_suffix']:
            OptionWarning('executables_suffix')
        if not os.path.isfile(self['scoreView_path']):
            self['scoreView_path'] = self.__findScoreView()
    def __checkPaths(self):
        if not os.path.isdir(self['executables_dir']): return
        for name, path in self.paths.items():
            if os.path.isfile(path): continue
            fname = '%s.%s' % (name, self['executables_suffix'])
            fpath = os.path.join(self['executables_dir'], fname)
            if os.path.isfile(fpath): self.paths[name] = fpath
            else: ExecutablePathWarning(name, fpath)

    def __findRosetta(self):
        bundles = []
        for path in constants.searchPaths:
            if not os.path.isdir(path): continue
            for d in os.listdir(path):
                if d.lower().startswith('rosetta3') and d.lower().endswith('bundles'):
                    bundles.append(os.path.join(path, d))
        if bundles: return sorted(bundles)[-1]
        OptionWarning('rosetta_bundles')
        return ''
    def __findScoreView(self):
        for path in constants.searchPaths:
            if not os.path.isdir(path): continue
            if 'ScoreView.app' in os.listdir(path):
                return os.path.realpath(os.path.join(
                    path, 'ScoreView.app', 'Contents', 'MacOS', 'ScoreView'))
        return ''
    def __findDb(self):
        path = os.path.join(self['rosetta_bundle'], 'rosetta_database')
        if not os.path.isdir(path):
            DirectoryPathWarning('rosetta_database', path)
            return ''
        return path
    def __findExecDir(self):
        path = os.path.join(self['rosetta_bundle'], 'rosetta_source',
                            'build', 'src', 'release')
        if os.path.isdir(path):
            ostype = os.listdir(path)[-1]
            path = os.path.join(path, ostype) # os-specific
            path = os.path.join(path, os.listdir(path)[-1]) # os-version specific
            path = os.path.join(path, os.listdir(path)[-1]) # 32/64 bit specific
            path = os.path.join(path, os.listdir(path)[-1]) # cpu interface specific
            cCompiler = os.listdir(path)[-1]
            path = os.path.join(path, cCompiler) # c-compiler specific
        if not os.path.isdir(path):
            DirectoryPathWarning('executables', path)
            return ''
        self['executables_suffix'] = ostype+cCompiler+'release'
        return path

