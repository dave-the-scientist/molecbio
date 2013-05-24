"""Application and module for viewing rosetta scores

    This is a self-contained application that can handle many different kinds of rosetta
application scores, and should work even for those it has not been designed for, as it
defines a general case for scores. Has been tested with: docking, relax, ab initio,
homology modelling, loopmodel, flexdock, and floppytail.

    This can also be used as a module to import this functionality into some other
application, by using the ViewFrame class. This inherits from tkinter's Frame, and should
be treated as such, but brings with it everything necessary to function.
"""
import sys, os, subprocess, threading
from Tkinter import *
import ttk, tkFileDialog, tkMessageBox
import util
from molecbio import Cluster
if not Cluster:
    print "Could not import the Cluster module, so its functionality has been disabled. This is usually because the 'numpy' module could not be found." 


class ViewFrame(Frame):
    def __init__(self, parent, **args):
        Frame.__init__(self, parent, **args)
        self.parent = parent
        # # #  GUI Frames
        self.columnconfigure(0, weight=1)
        self.rowconfigure(2, weight=1)
        self.scoresLabel = Label(self, text='Rosetta Scores', font=('Verdana', '18'))
        self.scoresLabel.grid(row=0, column=0, sticky=EW)
        self.messageLabel = Label(self, text='', font=('Verdana','13'))
        self.messageLabel.grid(row=1, column=0, sticky=W)
        self.scoresFrame = Frame(self, bd=1, relief=SUNKEN)
        self.scoresFrame.grid(row=2, column=0, sticky=NSEW)
        self.scoresFrame.columnconfigure(0, weight=1)
        self.scoresFrame.rowconfigure(0, weight=1)
        self.buttonFrame = Frame(self, bd=0)
        self.buttonFrame.grid(row=3, column=0, pady=10, sticky=EW)
        self.buttonFrame.columnconfigure((0,1,2,3), weight=1)
        # # #  GUI code
        self.scoresListTree = ScoreTree(self.scoresFrame, height=10, selectmode=EXTENDED,
                                             padding=(-1,-1), dblcommand=self.treeDoubleClick)
        #self.scoresListTree.bind('<Return>', self.openSeleCommand)
        #self.scoresListTree.bind('<KP_Enter>', self.openSeleCommand)
        self.scoresListTree.grid(row=0, column=0, sticky=NSEW)
        self.scoresListScroll = Scrollbar(self.scoresFrame, orient=VERTICAL, 
                                          command=self.scoresListTree.yview)
        self.scoresListScroll.grid(row=0, column=1, sticky=NS)
        self.scoresListTree.configure(yscrollcommand=self.scoresListScroll.set)
        
        self.openScoreButton = Button(self.buttonFrame, text='Load Results', 
                                      command=self.openScoreCommand)
        self.openScoreButton.grid(row=0, column=0)
        self.tableButton = Button(self.buttonFrame, text='Results to Text', 
                                  state=DISABLED, command=self.tableResultsCommand)
        self.tableButton.grid(row=0, column=1)
        self.openSeleButton = Button(self.buttonFrame, text='Open Selected',
                                     state=DISABLED, command=self.openSeleCommand)
        self.openSeleButton.grid(row=0, column=2)
        self.clusterButton = Button(self.buttonFrame, text='Cluster Results',
                                    state=DISABLED, command=self.clusterResultsCommand)
        self.clusterButton.grid(row=0, column=3)
        

    def treeDoubleClick(self, treeData):
        names = zip(*treeData)[0]
        paths, notFound = [], []
        for name in names:
            path = os.path.join(self.scoresDir, name)
            if not os.path.isfile(path): notFound.append(name)
            else: paths.append(path)
        if notFound:
            self.messageLabel.configure(text='Could not locate %s.' % ', '.join(notFound))
        if paths: util.startFiles(*paths)

    def openSeleCommand(self, event=None):
        data = self.scoresListTree.getSelectedData()
        if not data: return
        self.treeDoubleClick(data)
    def openScoreCommand(self, event=None):
        directory = tkFileDialog.askdirectory(
            parent=self.parent, mustexist=True, title='Choose Results Folder',
            message='Choose a directory containing the score file and the decoys from some Rosetta run.')
        if not directory: return
        self._openScores(directory)
    def clusterResultsCommand(self):
        ClusterWindow(self)
    def tableResultsCommand(self):
        def floatIfIs(score):
            try: return '%.3f'%score
            except (ValueError, TypeError): return str(score)
        metrics = self.scoresListTree['columns']
        if not self.scoresDict or not metrics: return
        header = '\t'.join(metrics)
        buff = [header]
        for name, score in self.scoresDict.items():
            buff.append('\t'.join(floatIfIs(score.get(metric, '')) for metric in metrics))
        tblWidth = max(len(header.expandtabs()), len(buff[1].expandtabs())) + 2
        tblText = '\n'.join(buff)
        self.__showTabulateWindow(tblWidth, tblText)

    # # # # #  Private Methods  # # # # #
    def __fillTree(self):
        scoreType = self.scoresDict.scoreType()
        if scoreType == 'docking':
            headerList = [('rms','RMS'),('I_sc','Interface'),('total_score','Total Score')]
        elif scoreType == 'abinitio':
            headerList = [('rama','Ramachandran'),('omega','Omega'),('score','Score')]
        elif scoreType == 'floppytail':
            headerList = [('rama','Ramachandran'),('omega','Omega'),
                          ('total_score','Total Score')]
        elif scoreType == 'homology':
            headerList = [('rama','Ramachandran'),('omega','Omega'),
                          ('total_score','Total Score')]
        elif scoreType == 'loopmodel':
            headerList = [('total_energy','Total Energy')]
        else: headerList = [('score','Score'),('total_score','Total Score')]
        self.scoresListTree.setupColumns(headerList)
        self.scoresListTree.fillTree(self.scoresDict)
        self.scoresListTree.sortTree(len(headerList))
        if Cluster: self.clusterButton.config(state=NORMAL)
        self.tableButton.config(state=NORMAL)
        self.openSeleButton.config(state=NORMAL)
    def __clearTree(self):
        self.scoresListTree.clearTree()
        self.clusterButton.config(state=DISABLED)
        self.tableButton.config(state=DISABLED)
        self.openSeleButton.config(state=DISABLED)
    def __showTabulateWindow(self, tblWidth, tblText):
        dialog = Toplevel()
        dialog.transient(self.parent)
        dialog.minsize(291, 291)
        dialog.geometry('291x291+100+100')
        Message(dialog, aspect=600, text='Press Select All to highlight everything; you can then copy and paste as usual.').grid(row=0, column=0, columnspan=2, padx=15, pady=10, sticky=EW)
        tableText = Text(dialog, bg='white', width=tblWidth, height=10)
        tableText.insert(1.0, tblText)
        tableText.grid(row=1, column=0, columnspan=2, padx=15, pady=10, sticky=NS)
        tableText.focus_set(); tableText.config(state=DISABLED)
        Button(dialog, text='Select All', command=lambda:
               ( tableText.tag_remove(SEL, 1.0, END), tableText.tag_add(SEL, 1.0, END) )
               ).grid(row=2, column=0, pady=10)
        Button(dialog, text='Close', command=dialog.destroy).grid(row=2, column=1, pady=10)
        dialog.rowconfigure(1, weight=1)
        dialog.columnconfigure((0,1), weight=1)
        dialog.lift()
        self.tableButton.wait_window(dialog)
        
    def _openScores(self, directory):
        self.__clearTree()
        dirName = os.path.basename(directory)
        self.scoresLabel.configure(text="Loading %s..." % dirName)
        self.scoresDict = util.ScoresDict(directory)
        if not self.scoresDict:
            self.messageLabel.configure(text="Failed to parse %s."%directory)
            self.scoresDir = None
            return
        self.__fillTree()
        self.scoresDir = directory
        scoreType = self.scoresDict.scoreType()
        if scoreType == 'unknown': label = 'Scores for %s' % dirName
        else: label = '%s scores for %s' % (scoreType.capitalize(), dirName)
        self.scoresLabel.configure(text=label)
        if self.messageLabel.cget('text').startswith('Failed to parse'):
            self.messageLabel.configure(text='')


class ClusterWindow(Toplevel):
    def __init__(self, parent):
        Toplevel.__init__(self)
        self.parent = parent
        self.cluster = Cluster.Clusterer()
        self.scoresDict = parent.scoresDict
        self.scoresDir = parent.scoresDir
        self.clustered = []
        self.__lastSort = None

        #self.transient(self.parent.parent)
        self.columnconfigure(0, weight=1); self.rowconfigure(0, weight=1)
        self.minsize(365, 235); self.geometry('900x350+100+100')

        self.leftFrame = Frame(self)
        self.leftFrame.grid(row=0, column=0, padx=15, pady=15, sticky=NSEW)
        self.leftFrame.columnconfigure(0, weight=1)
        self.leftFrame.rowconfigure(0, weight=1)
        self.tree = ScoreTree(self.leftFrame, height=10, selectmode=EXTENDED, padding=(-1,-1),
                                   dblcommand=self.__treeDoubleClick)
        self.tree.grid(row=0,column=0, sticky=NSEW)
        self.tree.bind('<Return>', self.__openSeleCommand)
        self.tree.bind('<KP_Enter>', self.__openSeleCommand)
        self.treeScroll = Scrollbar(self.leftFrame, orient=VERTICAL, 
                                    command=self.tree.yview)
        self.treeScroll.grid(row=0, column=1, sticky=NS)
        self.tree.configure(yscrollcommand=self.treeScroll.set)

        self.rightFrame = Frame(self)
        self.rightFrame.grid(row=0, column=1, padx=(0,15), pady=15, sticky=NSEW)
        self.rightFrame.rowconfigure(5, weight=1)
        Label(self.rightFrame, text='Results Folder:').grid(row=0, column=0, pady=5, sticky=W)
        self.resultsDirButton = Button(self.rightFrame, text='Choose',
                                       command=self.__chooseResultsDir)
        self.resultsDirButton.grid(row=0, column=1, pady=5, sticky=E)
        self.resultsDirLabel = Label(self.rightFrame, text=self.scoresDir, relief=RIDGE,
                                     bd=1, width=30, anchor=E)
        self.resultsDirLabel.grid(row=1, column=0, columnspan=2, padx=(20,0),
                                  pady=(0,10), sticky=EW)
        Label(self.rightFrame, text='Top Decoys to Cluster:\n(0 to cluster them all)').grid(
            row=2, column=0, pady=5, sticky=W)
        self.numDecoysEntry = Entry(self.rightFrame, width=8, justify=RIGHT, bg='white',
                                    bd=2, relief=SUNKEN)
        self.numDecoysEntry.grid(row=2, column=1, pady=5, sticky=E)
        self.numDecoysEntry.insert(0, '50')
        Label(self.rightFrame, text='Cluster radius:').grid(row=3, column=0, sticky=W)
        self.radiusEntry = Entry(self.rightFrame, width=8, justify=RIGHT, bg='white',
                                 bd=2, relief=SUNKEN)
        self.radiusEntry.grid(row=3, column=1, sticky=E)
        self.radiusEntry.insert(0, '4.5')
        self.runButton = Button(self.rightFrame, text='Run Analysis', width=16,
                                command=self.__runClusteringCommand)
        self.runButton.grid(row=4, column=0, pady=15, columnspan=2)
        self.progBar = ttk.Progressbar(self.rightFrame)
        self.progBar.grid(row=5, column=0, columnspan=2, sticky=S+EW)
        self.progLabel = Label(self.rightFrame, font=('Helvetica', '12'))
        self.progLabel.grid(row=6, column=0, columnspan=2, pady=5, sticky=S+EW)
        self.progLabel.grid_remove(); self.progBar.grid_remove()
        butFrame = Frame(self.rightFrame); butFrame.grid(row=7, column=0, columnspan=2, sticky=S+EW)
        butFrame.columnconfigure((0,1), weight=1)
        Button(butFrame, text='Open Selected', command=self.__openSeleCommand).grid(row=0,column=0)
        Button(butFrame, text='Quit', command=self.destroy).grid(row=0, column=1)

    def __chooseResultsDir(self):
        # Adds 1-3 MB of ram every time is clicked, and when browsing directories.
        # Adds another 2MB or so when a dir is actually picked.
        directory = tkFileDialog.askdirectory(parent=self, title='Choose Results Folder',
                                              mustexist=True, initialdir=self.scoresDir,
                                              message='Choose a directory containing the .fasc or .sc file and the decoys from some Rosetta run.')
        if directory:
            scoresDict = util.ScoresDict(directory)
            if scoresDict:
                self.resultsDirLabel.config(text=directory)
                self.scoresDir = directory
                self.scoresDict = scoresDict
                self.__clearTree()
        
    def __runClusteringCommand(self):
        t = threading.Thread(target=self.__runClustering)
        t.start()

    def __runClustering(self):
        ### Do the progress bars using rescheduling events; root.after_idle or root.after.
        #   Should be a lot more efficient than causing the clusterer to keep track of counters.
        readProgStep = 1
        cmpProgStep = 50
        radius, directory, targets = self.__getInputs()
        if not radius: return False
        self.runButton.config(state=DISABLED)
        self.progLabel.config(text='Loading files...'); self.progLabel.grid()
        self.progBar.grid(); self.progBar.config(mode='indeterminate')
        self.progBar.start(); self.progBar.update_idletasks()
        self.cluster._readProgressStep = readProgStep
        self.cluster._readProgressFxn = self.__readProgressFxn(readProgStep).next
        self.cluster._cmpProgressStep = cmpProgStep
        self.cluster._cmpProgressFxn = self.__cmpProgressFxn(cmpProgStep).next

        clustered = self.cluster.cluster(radius, directory, targets)
        #clustered = self.cluster.cluster(radius, directory, targets, progStep,
        #                                 self.__progressFxn(progStep).next)
        if not clustered:
            util.popupError('Problem Clustering Files', 'Something went wrong trying to'+
                       ' cluster the specified files.')
        else:
            self.clustered = [sorted(cluster, key=lambda f: self.scoresDict.score(f))
                              if type(cluster) is list else cluster for cluster in clustered]
            self.__fillTree()
        self.progBar.grid_remove()
        self.progLabel.grid_remove()
        self.runButton.config(state=NORMAL)

    def __getInputs(self):
        try:
            radius = float(self.radiusEntry.get().strip())
            num = int(self.numDecoysEntry.get().strip())
        except ValueError:
            util.popupError('Problem Reading Options',
                       'One or more of the options were either left blank'+
                       ' or filled inappropriately.')
            return 0,0,0
        if num > 0 and num < len(self.scoresDict):
            targets = self.scoresDict.sort()[:num]
        else:
            targets = None
        directory = self.resultsDirLabel.cget('text')
        return radius, directory, targets

    def __readProgressFxn(self, progStep):
        numDone = 0
        num = len(self.scoresDict)
        num = min(int(self.numDecoysEntry.get().strip()), num) or num
        self.progBar.stop()
        self.progBar.config(mode='determinate', maximum=num)
        while numDone < num:
            numDone += progStep
            self.progBar.config(value=numDone)
            self.progLabel.config(text='%i of %i files loaded.' % (numDone, num))
            yield numDone
    def __cmpProgressFxn(self, progStep):
        """Little hack. Not used as a generator fxn, just a function with persistant
        memory."""
        total = 0
        num = len(self.scoresDict)
        num = min(int(self.numDecoysEntry.get().strip()), num) or num
        numCmps = num*(num-1)/2
        self.progBar.config(mode='determinate', maximum=numCmps)
        while total < numCmps:
            total += progStep
            self.progBar.config(value=total)
            self.progLabel.config(text='%d%% of %d calculations.' % (total*100/numCmps, numCmps))
            yield total
    
    def __openSeleCommand(self, event=None):
        data = self.tree.getSelectedData()
        if not data: return
        self.__treeDoubleClick(data)
    def __treeDoubleClick(self, treeData):
        names = zip(*treeData)[0]
        paths, notFound = [], []
        for name in names:
            path = os.path.join(self.scoresDir, name)
            if not os.path.isfile(path): notFound.append(name)
            else: paths.append(path)
        if notFound:
            self.messageLabel.configure(text='Could not locate %s.' % ', '.join(notFound))
        if paths: util.startFiles(*paths)

    def __fillTree(self):
        """self.clustered should be filled out, and each cluster sorted by score."""
        # Make this smarter. Automatically decide header, instead of hard-coding here.
        scoreType = self.scoresDict.scoreType()
        if scoreType == 'docking':
            headerList = [('I_sc','Interface'),('total_score','Total Score'),
                          ('clusterSize','Cluster Size')]
        elif scoreType == 'abinitio':
            headerList = [('rama','Ramachandran'),('score','Score'),
                          ('clusterSize','Cluster Size')]
        elif scoreType == 'floppytail':
            headerList = [('rama','Ramachandran'),('total_score','Total Score'),
                          ('clusterSize','Cluster Size')]
        elif scoreType == 'homology':
            headerList = [('rama','Ramachandran'),('total_score','Total Score'),
                          ('clusterSize','Cluster Size')]
        else:
            #headerList = [('score','Score'),('total_score','Total Score'),
            #                ('clusterSize','Cluster Size')]
            l = [('chainbreak','Chain Break'),('score','Score'),
                 ('total_score','Total Score'), ('total_energy','Energy')]
            headerList = [t for t in l if t[0] in self.scoresDict.metrics()] + [
                ('clusterSize','Cluster Size')]
        d = {}
        for node in self.clustered:
            if type(node) is list:
                for name in node:
                    d[name] = self.scoresDict[name].copy()
                d[node[0]]['clusterSize'] = len(node)
            else:
                d[node] = self.scoresDict[node].copy()
        self.tree.setupColumns(headerList)
        self.tree.fillTree(d, self.clustered)
        self.tree.sortTree(2)
        
    def __clearTree(self):
        self.tree.clearTree()
        self.progBar.grid_remove()
        self.runButton.config(state=NORMAL)

class ScoreTree(ttk.Treeview):
    """ A modified version of the ttk TreeView object. Implements a few methods and variables
    as well as introducing the option 'dblcommand' which should be a function that will be
    called when an entry is double clicked, with the entry's name as the argument. """
    def __init__(self, parent, **args):
        if 'dblcommand' in args:
            self.__dblClickFxn = args.pop('dblcommand')
        else: self.__dblClickFxn = None
        ttk.Treeview.__init__(self, parent, **args)
        self.bind('<Double-Button-1>', self.__dblClick)
        self.bind('<Return>', self.__dblClick)
        self.bind('<KP_Enter>', self.__dblClick)
        # # # # # # # # # #  Variables  # # # # # # # # # #
        self.__data = {}
        self.__lastSort = None
        self.__defaultWidth = 500

    def setupColumns(self, headerList, width=None):
        if not width: width = self.winfo_width()
        if width == 1: width = self.__defaultWidth
        unitWidth = width / (len(headerList)+2) - 15
        self['columns'] = [col[0] for col in headerList]
        widths = 0
        for i, (name,label) in enumerate(headerList):
            i += 1
            col = '#%i'%i
            minWidth = max(55, len(label)*10 - 30)
            colWidth = max(minWidth, unitWidth)
            self.heading(col, text=label, command=lambda ind=i: self.sortTree(ind))
            self.column(col, anchor=CENTER, width=colWidth, minwidth=minWidth)
            widths += colWidth
        filenameWidth = width - widths - 2
        self.heading('#0', text='Filename', command=lambda: self.sortTree(0))
        self.column('#0', anchor=CENTER, width=filenameWidth, minwidth=200)
          
    def fillTree(self, attDict, namesList=None):
        # Should be passed a list of filenames; if a node is to have children then
        # that entry should be a sublist of names instead, where the first is the parent.
        # If no namesList is provided, all entries in the attDict will be added without
        # children.
        if not namesList:
            namesList = list(attDict)
        idDict = {}
        attKeys = self['columns']
        def insertByName(name, parent=''):
            vals = [attDict[name].get(key, '') for key in attKeys]
            valStrs = ['%.3f'%v if type(v) == float else str(v) for v in vals]
            idNum = self.insert(parent,END, text=name, values=valStrs, tags='treeItem')
            idDict[idNum] = [name]+vals
            return idNum
        self.clearTree()
        for entry in namesList:
            if type(entry) is str:
                insertByName(entry)
            elif type(entry) is list:
                idNum = insertByName(entry[0])
                for name in entry[1:]:
                    insertByName(name, idNum)
        self.__data = idDict
        # self.tag_bind('treeItem', '<Double-Button-1>', self.__dblClick)

    def sortTree(self, valueIndex):
        if self.__lastSort == valueIndex:
            children = list(self.get_children())
            children.reverse()
        else:
            children = sorted(self.get_children(),
                              key=lambda idNum: self.__data[idNum][valueIndex])
            self.__lastSort = valueIndex
        self.set_children('', *children)
    def clearTree(self):
        self.delete(*self.get_children())
        self.__data = {}
        self.__lastSort = None

    def getSelectedData(self):
        sele = self.selection()
        if not sele: return False
        return map(self.__data.get, sele)

    # # # # # # # # # #  Private Functions  # # # # # # # # # #
    def __getitem__(self, key):
        if type(key) is int:
            return self.__data[self.get_children()[key]]
        else:
            return ttk.Treeview.__getitem__(self, key)
    def __dblClick(self, event=None):
        data = self.getSelectedData()
        if not callable(self.__dblClickFxn) or not data: return False
        self.__dblClickFxn(data)


class GUI(object):
    def __init__(self, root):
        self.root = root
        # # #  Menus
        self.menubar = Menu(self.root)
        helpMenu = Menu(self.menubar)
        helpMenu.add_command(label='ScoreView Help', command=self.helpMenuCommand)
        self.menubar.add_cascade(label='Help', menu=helpMenu)
        self.root.config(menu=self.menubar)

        self.mainScreen = Frame(self.root, bd=1, relief=RAISED, padx=0, pady=0)
        self.mainScreen.pack(expand=YES, fill=BOTH)
        self.mainScreen.columnconfigure(0, weight=1)
        self.mainScreen.rowconfigure(0, weight=1)

        self.viewFrame = ViewFrame(self.mainScreen)
        self.viewFrame.grid(row=0, column=0, padx=20, pady=10, sticky=NSEW)
    def openScores(self, resultsDir):
        if os.path.isdir(resultsdir):
            self.viewFrame._openScores(resultsDir)
    def helpMenuCommand(self):
        helpBox = Toplevel(padx=15, pady=15)
        helpBox.transient(self.root)
        helpBox.minsize(250,155); helpBox.geometry('775x650+10+30')
        helpBox.rowconfigure(1, weight=1)
        helpBox.columnconfigure(0, weight=1)
        Label(helpBox, text='ScoreView User Manual', font=('helvetica', 18)).grid(
            row=0, column=0, columnspan=2, pady=(0,10))
        helpText = Text(helpBox, bg='white', bd=2, relief=SUNKEN, highlightthickness=0,
                        font=('helvetica', 14), padx=10, spacing1=10, spacing2=2,
                        spacing3=10, wrap=WORD)
        helpText.grid(row=1, column=0, sticky=NSEW)
        helpText.insert(1.0, mainHelpMessage)
        helpText.config(state=DISABLED)
        scrollBar = Scrollbar(helpBox, orient=VERTICAL, command=helpText.yview)
        scrollBar.grid(row=1, column=1, sticky=NS)
        helpText.configure(yscrollcommand=scrollBar.set)
        helpBox.lift()


mainHelpMessage = '    Very important note: This program has been tested with PyMol for Mac version 1.3; other versions of PyMol may or may not work. It is however known that PyMol version 0.99 DOES NOT WORK.\n' + \
'    This program was written to parse the multitudes of results typically generated by some rosetta protein modelling application.\n' + \
'    When examining the results of some run, several scores are presented. Total Score is the metric used by the docking protocol, and is a combination of several scores for the complex, the lower the better. The interface Score is a separate measure, and describes the difference in score between each individual protein and that of their complex. Again the lower the better, with good models usually from -5 to -10. Other applications may have different scoring metrics, but at least one should be displayed. The results can be sorted, and double-clicking a decoy or clicking the "Open Selected" button will open it in your default pdb viewing program.\n' + \
'    The "Results to Text" button brings up a text box with the RMS value (if present) and the main scoring metric for that run separated by a tab, for each scored model. Pressing select all and then command-c will copy the data so you can paste it into your favourite graphing or spreadsheet program. This allows the creation of graphs overviewing the docking run, so funnel clusters can be viewed.\n' + \
'    Since each decoy is generated independantly, many of them may be very similar. Clustering the results attempts to separate the decoys into clusters of similar structures, allowing a better understanding of the results. This process is computationally intensive, and increases semi-exponentially the more decoys that are considered. Clustering between 50-200 is reasonable, and should take a minute or two at most to compute. However, if you are running the program on one computer, and opening results over a network that are stored somewhere else, loading the files can sometimes be very slow and can take between 0.5 and 5 seconds each. The cluster radius allows you to define how similar structures must be (in Angstrom) to be considered a cluster. Appropriate values depend on the size of the structures, and the magnitude by which the models are likely to be different. For docking proteins, 3-6 Angstrom is probably reasonable, while modelling a single loop might call for a radius of 0.5-2.5. Generally the size can be adjusted until a pattern emerges with a handful of well-populated clusters. Too many or too few clusters generally impart little information.\n'

def start(resultsDir=None):
    minWidth, minHeight = 500, 200
    startWidth, startHeight = 600, 500
    xCornerOffset, yCornerOffset = 0, 0
    
    root = Tk()
    root.title('Rosetta Score Viewer')
    root.minsize(minWidth, minHeight)
    root.geometry('%ix%i+%i+%i' % (
        startWidth, startHeight, xCornerOffset, yCornerOffset))
    root.config(bg='SystemButtonFace')
    root.option_add('*Background', 'SystemButtonFace')
    root.option_add('*Button*highlightBackground', 'SystemButtonFace')
    root.option_add('*Entry*highlightBackground', 'SystemButtonFace')
    root.option_add('*Text*Background', 'White')
    root.option_add('*Menubutton*Foreground', 'Black')
    root.option_add('*Menu*Foreground', 'Black')
    root.option_add('*tearOff', False)
    root.event_delete('<<PasteSelection>>', '<ButtonRelease-2>')
    program = GUI(root)

    if resultsDir and os.path.isdir(resultsDir):
        program.viewFrame._openScores(resultsDir)

    root.mainloop()

if __name__ == '__main__':
    if len(sys.argv) >= 2:
        scoresDir = sys.argv[1]
    else: scoresDir = None
    start(scoresDir)
