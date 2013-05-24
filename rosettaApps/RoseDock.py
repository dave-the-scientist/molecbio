#!/usr/bin/python
# Author: Dave Curran
# Date: July 2010
# Purpose: This is the GUI code that runs my rosetta framework to perform Rosetta
#   docking on protein structures. This program provides an easy-to-use interface to this
#   complicated procedure, as well as several features for parsing the result of a run.
#   It allows single amino acid constraints to be defined on either docking partner, and
#   ensures ad-hoc that no structure is kept that does not satisfy these constraints. Once
#   finished, resultant structures (decoys) are presented in a list, and can be examined by
#   the user's default pdb viewing program with a double-click. Decoys can be sorted by a
#   variety of metrics, and can also be clustered. This analysis shows how many unique
#   solutions have been found, and provides validation for certain conformations. More
#   detailed explanations can be found in the program's Help menu.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#### Look into a combo score, total + interface, or interface/2. test with known structs.
import os, threading, time
import rosedockScreens, util
import ScoreView
from molecbio import rosetta
from Tkinter import *
import ttk, tkFileDialog

class RoseDocker(rosetta.Docker):
    def __new__(cls, docker, gui):
        docker.__class__ = RoseDocker
        return docker
    def __init__(self, docker, gui):
        self.gui = gui
        self.openScoreViewer = False
        self._saveState = False
    def _initialFileSetup(self):
        self.gui.startRunButton.config(state=DISABLED)
        self.gui.dockProgBar.config(mode='indeterminate'); self.gui.dockProgBar.start()
        self.gui.progLabel.config(text='Initializing docking protocol...')
        self.gui.dockProgBar.grid(); self.gui.progLabel.grid()
        self.gui.dockProgBar.update_idletasks()
        rosetta.Docker._initialFileSetup(self)
    def _preRun(self):
        numDecoys = int(self['-nstruct'])
        numCPUs = self.numCPUs
        self.gui.dockProgBar.stop()
        self.gui.messageLabel.configure(text='%d processes currently running.' % numCPUs)
        self.gui.dockProgBar.config(mode='determinate', maximum=numDecoys, value=0)
        self.gui.progLabel.config(text='0 of %i completed.' % numDecoys)
    def _reportProgress(self, numDone, numTotal):
        self.gui.dockProgBar.config(value=numDone)
        self.gui.progLabel.config(text='%i of %i completed.' % (numDone, numTotal))
    def _reportCompletion(self, runTime):
        mins, sec = divmod(runTime, 60)
        hrs, mins = divmod(mins, 60)
        total = int(self['-nstruct'])
        self.gui.messageLabel.configure(text='%i total decoys generated in %i hours, %i minutes, %i seconds.' % (total, hrs, mins, sec))
        self.gui.scoreView._openScores(self.outputDir)
    def _cleanUp(self):
        rosetta.Docker._cleanUp(self)
        self.gui.startRunButton.config(state=NORMAL)
        self.gui.dockProgBar.grid_remove()
        self.gui.progLabel.grid_remove()
        d = rosetta.load_options(os.path.join(self.outputDir, self.saveFile))
        if not d or d.execName != 'docking_protocol':
            d = rosetta.Docker()
        self.gui.docker = d
        self.gui._displayDockingOptions()
    def _exceptionCallback(self, excep):
        self.gui.messageLabel.configure(text='Docking run encountered an error before completion.')
        raise excep


class GUI:
    def __init__(self, root):
        self.parent = root
        self.docker = rosetta.Docker()
        # # # # # # # # # # #   Variables   # # # # # # # # # # #
        self.inputResidues = []
        self.inputFilepathVar = StringVar()
        self.outputFilepathVar = StringVar()
        self.numProcessesVar = StringVar(value='1')
        self.randStrucVar = IntVar(value=0)

        # # # # # # # # # # #   GUI Steup   # # # # # # # # # # #
        # # #  Major Frames  # # #
        self.mainScreen = Frame(self.parent, bd=1, relief=RAISED, padx=10, pady=10)
        self.mainScreen.pack(expand=YES, fill=BOTH)
        self.mainScreen.columnconfigure(0, weight=1)
        self.mainScreen.rowconfigure(2, weight=1)
        self.titleFrame = Frame(self.mainScreen)
        self.titleFrame.grid(row=0, column=0, columnspan=3, padx=5, sticky=EW)
        self.titleFrame.columnconfigure(0, weight=1)
        
        Frame(self.mainScreen, width=4, bd=2, relief=SUNKEN).grid(
            row=2, column=1, sticky=NS, pady=(35, 43)) # The vertical divider frame

        self.rightFrame = Frame(self.mainScreen)
        self.rightFrame.grid(row=2, column=2, padx=15, sticky=NSEW)
        self.rightFrame.rowconfigure(2, weight=1)
        self.inputFrame = Frame(self.rightFrame, height=385, width=300)
        self.inputFrame.grid(row=1, column=0, columnspan=3, sticky=EW)
        self.inputFrame.columnconfigure((0,1), weight=1)
        self.inputFrame.grid_propagate(0)

        # # #  Misc Widgets  # # #
        self.title = Label(self.titleFrame, text='RoseDock', font=('Helvetica', '48'))
        self.title.grid(sticky=EW)
        self.messageLabel = Label(self.mainScreen, text='', anchor=W, 
                                  font=('Helvetica', '14', 'bold'))
        self.messageLabel.grid(row=1, column=0, columnspan=3, padx=10,
                               pady=(0,5), sticky=EW)

        # # #  Results Display  # # #
        self.scoreView = ScoreView.ViewFrame(self.mainScreen)
        self.scoreView.grid(row=2, column=0, padx=15, sticky=NSEW)

        # # #  Docking Run Inputs  # # #
        Label(self.rightFrame, text='New Docking Run', font=('Helvetica', '16')).grid(
            row=0, column=0, columnspan=3, pady=(0,10), sticky=EW)
        Label(self.inputFrame, text='Input PDB:').grid(row=0, column=0, 
                                                       pady=5, sticky=W)
        self.inputPdbButton=Button(self.inputFrame, text='Choose', 
                                   command=self.inputPdbButtonCommand)
        self.inputPdbButton.grid(row=0, column=1, pady=5, sticky=E)
        self.inputFilepathLabel = Label(self.inputFrame, relief=RIDGE, bd=1, width=10,
                                  textvariable=self.inputFilepathVar, anchor=E)
        self.inputFilepathLabel.grid(row=1, column=0, columnspan=2, 
                                     padx=(20,0), pady=(0,10), sticky=EW)
        Label(self.inputFrame, text='Decoys to Generate:').grid(
            row=2, column=0, pady=5, sticky=W)
        self.decoysEntry = Entry(self.inputFrame, width=8, justify=RIGHT, 
                                 bg='white', bd=2, relief=SUNKEN)
        self.decoysEntry.grid(row=2, column=1, pady=5, sticky=E)
        Label(self.inputFrame, text='Parallel Processes:').grid(
            row=3, column=0, pady=5, sticky=W)
        self.numProcsEntry = Entry(self.inputFrame, width=4, justify=RIGHT, 
                                 bg='white', bd=2, relief=SUNKEN)
        self.numProcsEntry.grid(row=3, column=1, pady=5, sticky=E)
        Label(self.inputFrame, text='Randomize Structures?').grid(
            row=4, column=0, pady=5, sticky=W)
        self.randStrucCheckbutton = Checkbutton(self.inputFrame, 
                                                variable=self.randStrucVar)
        self.randStrucCheckbutton.grid(row=4, column=1, pady=5,  sticky=E)
        Label(self.inputFrame, text='Output to Folder:').grid(
            row=5, column=0, pady=5, sticky=W)
        self.outputFolderButton = Button(self.inputFrame, text='Choose', 
                                      command=self.outputFolderButtonCommand)
        self.outputFolderButton.grid(row=5, column=1, pady=5, sticky=E)
        self.outputFilepathLabel = Label(self.inputFrame, relief=RIDGE, bd=1, width=10,
                               textvariable=self.outputFilepathVar, anchor=E)
        self.outputFilepathLabel.grid(row=6, column=0, columnspan=2,
                                      padx=(20,0), pady=(0,10), sticky=EW)
        Label(self.inputFrame, text='Top decoys to keep:\n(0 to keep them all).',
              anchor=W).grid(row=7, column=0, pady=5, sticky=W)
        self.filterEntry = Entry(self.inputFrame, width=6, justify=RIGHT, bg='white',
                                 bd=2, relief=SUNKEN)
        self.filterEntry.grid(row=7, column=1, sticky=E)
        self.constLabel = Label(self.inputFrame)
        self.constLabel.grid(row=8, column=0, pady=10, sticky=W)
        self._updateConstraintsView()
        self.editConstButton = Button(self.inputFrame, text='Edit', state=DISABLED,
                                      command=self.editConstraintsCommand)
        self.editConstButton.grid(row=8, column=1, sticky=E)
        self.dockProgBar = ttk.Progressbar(self.inputFrame)
        self.dockProgBar.grid(row=9, column=0, columnspan=2, pady=(10,0), sticky=N+EW)
        self.progLabel = Label(self.inputFrame, font=('Helvetica', '12'), text='testing')
        self.progLabel.grid(row=10, column=0, columnspan=2, pady=0, sticky=N+EW)
        self.dockProgBar.grid_remove()
        self.progLabel.grid_remove()
        self.loadDockfileButton = Button(self.rightFrame, text='Load Options', 
                                         command=self.loadDockfileCommand)
        self.loadDockfileButton.grid(row=2, column=0, sticky=S, pady=5)
        self.saveDockfileButton = Button(self.rightFrame, text='Save Options',
                                         command=self.saveDockfileCommand)
        self.saveDockfileButton.grid(row=2, column=1, sticky=S, pady=5)
        self.resetButton = Button(self.rightFrame, text='Reset', 
                                  command=self.resetButtonCommand)
        self.resetButton.grid(row=2, column=2, sticky=S, pady=5)
        self.startRunButton = Button(self.rightFrame, text='Start Run', width=17,
                                     command=self.startRunButtonCommand)
        self.startRunButton.grid(row=3, column=0, columnspan=3, sticky=S, pady=5)
        self.resetButtonCommand()

        # # #  Top Menu Setup  # # #
        self.menubar = Menu(self.parent)
        appMenu = Menu(self.menubar, name='apple') #Application menu.
        fileMenu = Menu(self.menubar)
        editMenu = Menu(self.menubar)
        helpMenu = Menu(self.menubar)
        for menu,label,cmnd,accel in (
            (fileMenu,'Open .dockfile...',self.loadDockfileCommand,'Command-o'),
            (fileMenu,'Open Results Folder...',self.scoreView.openScoreCommand,'Command-Shift-O'),
            (fileMenu,'Save .dockfile...',self.saveDockfileCommand,'Command-s'),
            (editMenu,'Cut',lambda e=None: self.__focusEvent('<<Cut>>'),'Command-x'),
            (editMenu,'Copy',lambda e=None: self.__focusEvent('<<Copy>>'),'Command-c'),
            (editMenu,'Paste',lambda e=None: self.__focusEvent('<<Paste>>'),'Command-v'),
            (helpMenu,'RoseDock Help',self.helpMenuCommand,'') ):
            bindStr = '<%s>' % accel
            menu.add_command(label=label, command=cmnd, accelerator=accel)
            if menu == fileMenu: self.parent.bind_all(bindStr, cmnd)
        for label, menu in (('', appMenu), ('File', fileMenu), ('Edit', editMenu),
                            ('Help', helpMenu)):
            self.menubar.add_cascade(label=label, menu=menu)
        self.parent.config(menu=self.menubar)
        self.parent.createcommand('::tk::mac::ShowPreferences', lambda:rosedockScreens.PrefsWindow(self))

    # # # # # # # # # # # #   GUI Actions/Commands   # # # # # # # # # # # #
    def inputPdbButtonCommand(self):
        filepath = tkFileDialog.askopenfilename(parent=self.parent, filetypes=[('pdb Structure', '.pdb')], title='Select pdb file.', message='Please select a pdb file containing the two structures you wish to dock.')
        if not filepath: return
        self.__setInputPdb(filepath)

    def editConstraintsCommand(self):
        rosedockScreens.ConstraintsWindow(self)

    def outputFolderButtonCommand(self):
        dirPath = tkFileDialog.askdirectory(parent=self.parent, title='Select a Directory',
                message='Please choose a directory in which to save all of the decoy files.')
        if not dirPath: return
        if ' ' in dirPath:
            util.popupError('Problem wiith the path',
                            "Due to Rosetta's option handling, there can be no spaces in the output directory path.")
            return
        self.outputFilepathVar.set(dirPath)

    def startRunButtonCommand(self):
        opts = rosetta.options
        if not opts['rosetta_bundle'] or not opts['rosetta_database'] or not opts.paths['docking_protocol']:
            rosedockScreens.PrefsWindow(self)
            return
        if not self.__setDockingOptions():
            return
        self.docker.save()
        d = RoseDocker(self.docker, self)
        t = threading.Thread(target=d.run)
        t.start()

    def saveDockfileCommand(self, event=None):
        if not rosetta.options['rosetta_bundle']:
            rosedockScreens.PrefsWindow(self)
            return
        if not self.__setDockingOptions():
            return False
        self.__saveOptions()

    def loadDockfileCommand(self, event=None):
        filepath = tkFileDialog.askopenfilename(parent=self.parent,
                                                title='Choose an options file to load.',
                                                message='Choose an options file to load.')
        if not filepath: return
        d = rosetta.load_options(filepath)
        if not d or d.execName != 'docking_protocol': return False
        self.docker = d
        self._displayDockingOptions()
        self.messageLabel.configure(text='Options loaded from %s.'%os.path.basename(filepath))

    def resetButtonCommand(self):
        self.inputFilepathVar.set('')
        self.decoysEntry.delete(0,END); self.decoysEntry.insert(END, '1')
        self.numProcsEntry.delete(0,END)
        self.numProcsEntry.insert(END, str(rosetta.constants.default_numCPUs))
        self.filterEntry.delete(0,END)
        self.filterEntry.insert(END, '500')
        self.randStrucVar.set(0)
        self.outputFilepathVar.set('')
        for i in range(len(self.docker.constraints)): self.docker.constraints.pop()
        self._updateConstraintsView()

    def helpMenuCommand(self):
        helpBox = Toplevel(padx=15, pady=15)
        helpBox.transient(self.parent)
        helpBox.minsize(250,155); helpBox.geometry('750x625+10+30')
        helpBox.rowconfigure(1, weight=1)
        helpBox.columnconfigure(0, weight=1)
        Label(helpBox, text='RoseDock User Manual', font=('helvetica', 18)).grid(
            row=0, column=0, columnspan=2, pady=(0,10))
        helpText = Text(helpBox, bg='white', bd=2, relief=SUNKEN, highlightthickness=0,
                        font=('helvetica', 14), padx=10, spacing1=10, spacing2=2,
                        spacing3=10, wrap=WORD)
        helpText.grid(row=1, column=0, sticky=NSEW)
        helpText.insert(1.0, util.rosedockMainHelpMessage)
        helpText.config(state=DISABLED)
        scrollBar = Scrollbar(helpBox, orient=VERTICAL, command=helpText.yview)
        scrollBar.grid(row=1, column=1, sticky=NS)
        helpText.configure(yscrollcommand=scrollBar.set)
        helpBox.lift()

    # # # # # # # # # # # # #   Docking Functions   # # # # # # # # # # # # #
    def __setDockingOptions(self):
        problem = None
        inputPdb = self.inputFilepathVar.get()
        decoys = self.decoysEntry.get().strip()
        procs = self.numProcsEntry.get().strip()
        outputDir = self.outputFilepathVar.get()
        toFilter = self.filterEntry.get().strip()
        if not os.path.isfile(inputPdb): problem = 'Input PDB'
        elif not decoys.isdigit(): problem = 'Number of decoys to generate'
        elif not procs.isdigit(): problem = 'Number of processes'
        elif not os.path.isdir(outputDir): problem = 'Output to folder'
        elif not toFilter.isdigit(): problem = 'Number of decoys to keep'
        if problem:
            util.popupError('Problem wiith the options',
                            "The run was not started because of a problem with the '%s' option." % problem)
            return False
        self.docker.inputPdb = inputPdb
        self.docker.numStruct = decoys
        self.docker.numCPUs = procs
        self.docker.randomize = self.randStrucVar.get()
        self.docker.outputDir = outputDir
        self.docker.keepTopDecoys = toFilter
        return True
    def __saveOptions(self):
        filepath = tkFileDialog.askdirectory(parent=self.parent, title='Select a Directory',
                message='Please choose a directory in which to save the docking options.',
                initialdir=self.docker.outputDir)
        if not filepath: return
        self.docker.save(filepath)
        self.messageLabel.configure(text='Options saved to %s.' % os.path.basename(filepath))
    def _displayDockingOptions(self):
        self.__setInputPdb(self.docker.inputPdb)
        self.decoysEntry.delete(0,END); self.decoysEntry.insert(END,self.docker['-nstruct'])
        self.numProcsEntry.delete(0,END); self.numProcsEntry.insert(END,str(self.docker.numCPUs))
        self.randStrucVar.set(self.docker.randomize)
        self.outputFilepathVar.set(self.docker.outputDir)
        self.filterEntry.delete(0,END); self.filterEntry.insert(END,str(self.docker.keepTopDecoys))
    def __setInputPdb(self, filepath):
        if filepath and os.path.isfile(filepath):
            state = NORMAL
            self.inputResidues = util.residuesFromPdb(filepath)
            self.__checkConstraints()
        else:
            filepath = ''
            state = DISABLED
        self.inputFilepathVar.set(filepath)
        self.editConstButton.config(state=state)
        self._updateConstraintsView()
    def __checkConstraints(self):
        chains, residues = self.inputResidues[0], self.inputResidues[1:]
        for const in self.docker.constraints[:]:
            res, fc, tc, dist = const[0][:-1], const[0][-1], const[1], const[2]
            if fc not in chains or tc not in chains or res not in residues[chains.index(fc)]:
                self.docker.constraints.remove(const)
    def _updateConstraintsView(self):
        num = str(len(self.docker.constraints) or 'No')
        self.constLabel.configure(text=num + ' docking constraints defined.')

    def _startDockingProcesses_OLD(self):
        num = int(self.controller.userOptions['numProcesses'])
        self.messageLabel.configure(text='%d processes currently running.' % num)
        self.startRunButton.config(state=DISABLED)
        self.dockProgBar.config(mode='indeterminate'); self.dockProgBar.start()
        self.progLabel.config(text='Initializing docking protocol...')
        self.dockProgBar.grid(); self.progLabel.grid()
        self.dockProgBar.update_idletasks()
        self.__dockingRunning = True
        t = threading.Thread(target=self.__monitorDockingProgress)
        t.start()
        startTime = time.time()
        try:
            runSuccess = self.controller.startDockingRun()
        finally:
            self.__dockingRunning = False
            self.messageLabel.configure(text='Docking run encountered an error before completion.')
            self.startRunButton.config(state=NORMAL)
            self.dockProgBar.grid_remove()
            self.progLabel.grid_remove()
        if runSuccess:
            runTime = time.time()-startTime
            mins, sec = divmod(runTime, 60)
            hrs, mins = divmod(mins, 60)
            self.messageLabel.configure(text='%i total decoys generated in %i hours, %i minutes, %i seconds.' % (runSuccess, hrs, mins, sec))
            self.scoreView._openScores(self.controller.userOptions['outputLocation'])
    def __monitorDockingProgress_OLD(self):
        dirPath = self.controller.userOptions['outputLocation']
        name = os.path.basename(self.controller.userOptions['inputPdb'])[:-4]
        numDecoys = int(self.controller.userOptions['numDecoys'])
        numPdbs, numFinished = 0, 0
        self.progLabel.config(text='Running docking protocol...')
        while self.__dockingRunning:
            time.sleep(1)
            if os.path.isdir(dirPath):
                numPdbs = len(filter(lambda f: f.startswith(name) and f.endswith('.pdb'),
                                         os.listdir(dirPath)))
            if numPdbs != numFinished:
                if numFinished == 0:
                    self.dockProgBar.stop()
                    self.dockProgBar.config(mode='determinate', maximum=numDecoys)
                numFinished = numPdbs
                self.dockProgBar.config(value=numFinished)
                self.progLabel.config(text = '%i of %i completed.' % (numFinished, numDecoys))

    # # # # # # # # # #  Misc Functions  # # # # # # # # # #
    def __focusEvent(self, event):
        w = self.parent.focus_get()
        if w:
            w.event_generate(event)


root = Tk()
root.title('RoseDock - Protein Docking using Rosetta')
root.minsize(856, 605); root.geometry('920x605+0+0')
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

root.mainloop()
