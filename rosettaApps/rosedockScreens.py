#!/usr/bin/python
# Author: Dave Curran
# Date: October 2011
# Purpose: All of the non-main GUI screens used by RoseDock. Includes the Prefences
#   and Constraints windows.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
from Tkinter import *
import ttk
import tkFileDialog
import util
from molecbio import rosetta
import threading, os

class PrefsWindow(Toplevel):
    def __init__(self, parent):
        Toplevel.__init__(self)
        self.parent = parent
        self.options = rosetta._options.OptionsDict()
        self.transient(self.parent.parent)
        self.minsize(335, 350); self.maxsize(1600,350)
        self.geometry('415x350+100+100')
        self.title('Set paths to Rosetta')

        Button(self, text='?', command=self.__helpCmnd).grid(
            row=0, column=0, padx=15, pady=(10,0), sticky=E)
        f = LabelFrame(self, text='Configure Rosetta file paths', labelanchor=NW, padx=20, pady=10)
        f.grid(row=1,column=0, padx=15, pady=0, sticky=EW)
        f.columnconfigure(0, weight=1)
        self.columnconfigure(0, weight=1)

        labGrid = {'column':0, 'sticky':W}
        butGrid = {'column':1, 'sticky':E}
        for i, (text, cmnd) in enumerate((('Path to Rosetta_bundles', self.__bundleCmnd),
                                          ('Path to rosetta database', self.__dbCmnd),
                                          ('Path to docking executable', self.__dockingCmnd))):
            Label(f, text=text).grid(row=i*2, **labGrid)
            Button(f, text='Change', command=cmnd).grid(row=i*2, **butGrid)
        pathLabOpts = {'bd':1, 'relief':RIDGE, 'width':40, 'anchor':E}
        pathLabGridOpts = {'column':0, 'columnspan':2, 'sticky':EW, 'padx':(3,0), 'pady':(5,20)}
        self.bundleLabel = Label(f, text=self.options['rosetta_bundle'], **pathLabOpts)
        self.bundleLabel.grid(row=1, **pathLabGridOpts)
        self.dbLabel = Label(f, text=self.options['rosetta_database'], **pathLabOpts)
        self.dbLabel.grid(row=3, **pathLabGridOpts)
        self.dockingLabel = Label(f, text=self.options.paths['docking_protocol'], **pathLabOpts)
        self.dockingLabel.grid(row=5, **pathLabGridOpts)
        Button(self, text='Save Paths', command=self.__okCmnd, width=10
               ).grid(row=2, column=0, padx=20, pady=5, sticky=E)

    def __dockingCmnd(self):
        filepath = self.__findPath(self.dockingLabel, 'Find docking executable',
                                   'Find the appropriate docking executable.', False)
    def __dbCmnd(self):
        dirPath = self.__findPath(self.dbLabel, 'Find database folder',
                                  'Find the rosetta_database directory.')
    def __bundleCmnd(self):
        dirPath = self.__findPath(self.bundleLabel, 'Find Rosetta folder',
                                  'Find the main Rosetta_Bundles directory.')
        if not dirPath: return
        self.options['rosetta_bundle'] = dirPath
        self.options.validate()
        self.dbLabel.config(text=self.options['rosetta_database'])
        self.dockingLabel.config(text=self.options.paths['docking_protocol'])

    def __okCmnd(self):
        self.options.validate()
        self.options.save()
        rosetta.options.load()
        self.parent.docker = rosetta.Docker()
        self.parent._displayDockingOptions()
        self.destroy()
        
    def __helpCmnd(self):
        util.popupInfo('Help Information', util.pathsHelpMessage)
    def __findPath(self, lbl, title, message, isDirectory=True):
        if isDirectory:
            filepath = tkFileDialog.askdirectory(parent=self, title=title,
                                                 message=message, mustexist=True)
        else:
            filepath = tkFileDialog.askopenfilename(parent=self, title=title, message=message)
        if not filepath: return ''
        lbl.config(text=filepath)
        return filepath

class ConstraintsWindow(Toplevel):
    def __init__(self, parent):
        Toplevel.__init__(self)
        self.parent = parent
        self.transient(self.parent.parent)
        self.protocol("WM_DELETE_WINDOW", self.__closeConstWindow)
        self.columnconfigure(0, weight=1); self.rowconfigure(1, weight=1)
        self.minsize(365, 235); self.geometry('415x270+100+100')
        
        self.f = Frame(self, bd=2, padx=10, pady=10, relief=GROOVE)
        self.f.grid(row=0, column=0, columnspan=2, padx=15, pady=15, sticky=EW)
        self.f.columnconfigure((0,1,2,3), weight=1)
        self.fcVar,self.resVar = StringVar(value='Chain'),StringVar(value='Residue')
        self.tcVar = StringVar(value='Chain')
        for i, text in enumerate(('From Chain:', 'Residue:', 'To Chain:', 'Distance:')):
            Label(self.f, text=text).grid(row=0,column=i, sticky=W)
        self.fromChain = OptionMenu(self.f, self.fcVar,'')
        self.fromChain.configure(width=9)
        self.fromChain.grid(row=1, column=0, sticky=W)
        self.resOM = OptionMenu(self.f, self.resVar, '')
        self.resOM.configure(width=10, state=DISABLED)
        self.resOM.grid(row=1, column=1, sticky=W)
        self.toChain = OptionMenu(self.f, self.tcVar, '')
        self.toChain.config(width=8, state=DISABLED)
        self.toChain.grid(row=1, column=2, sticky=W)
        self.distEnt = Entry(self.f, width=7, bg='white', bd=2, relief=SUNKEN)
        self.distEnt.insert(END,'10'); self.distEnt.grid(row=1, column=3, sticky=W)

        self.fromChain['menu'].delete(0,END)
        for i, chain in enumerate(self.parent.inputResidues[0]):
            self.fromChain['menu'].add_command(label=chain,
                                    command=lambda ind=i,ch=chain: (self.fcVar.set(ch),
                                    self.__fcCommand(ind)) )
        
        self.constTree = ttk.Treeview(self, columns=('fromChain', 'distance','toChain'),
                                 height=5, selectmode=BROWSE, padding=(-1,-1))
        self.constTree.grid(row=1, column=0, padx=(15,0), sticky=NSEW)
        self.constTree.heading('#0', text='Residue')
        self.constTree.heading('fromChain', text='From Chain')
        self.constTree.heading('distance', text='Distance')
        self.constTree.heading('toChain', text='To Chain')
        self.constTree.column('#0', width=80)
        self.constTree.column('fromChain', anchor=CENTER, width=70)
        self.constTree.column('distance', anchor=CENTER, width=60)
        self.constTree.column('toChain', anchor=CENTER, width=60)
        self.constScroll = Scrollbar(self, orient=VERTICAL, 
                                     command=self.constTree.yview)
        self.constScroll.grid(row=1, column=1, padx=(0,15), sticky=NS)
        self.constTree.configure(yscrollcommand=self.constScroll.set)
        self.__fillConstTree()
        self.bf = Frame(self)
        self.bf.grid(row=2, column=0, columnspan=2, padx=15, pady=15, sticky=EW)
        self.bf.columnconfigure((0,1,2), weight=1)
        Button(self.bf, text='Add Constraint', command=self.__addConstCommand).grid(row=0,column=0)
        Button(self.bf, text='Remove Constraint', command=self.__rmvConstCommand).grid(row=0,column=1)
        Button(self.bf, text='Done', width=10, command=self.__closeConstWindow).grid(row=0,column=2)
    def __fillConstTree(self):
        children = self.constTree.get_children()
        if children: self.constTree.delete(*children)
        for (res, tochain, dist) in self.parent.docker.constraints:
            fromchain, res = res.strip()[-1], res.strip()[:-1]
            self.constTree.insert('',END, text=res, values=(fromchain,dist,tochain))
        
    def __fcCommand(self, chainIndex):
        resLabels = self.parent.inputResidues[chainIndex+1]
        tcLabels = self.parent.inputResidues[0][:]
        tcLabels.remove(self.parent.inputResidues[0][chainIndex])
        for om, omVar, labels in ((self.resOM,self.resVar,resLabels),
                                  (self.toChain,self.tcVar,tcLabels)):
            om.config(state=NORMAL)
            om['menu'].delete(0,END)
            for name in labels:
                om['menu'].add_command(label=name, command=lambda var=omVar,nm=name: var.set(nm))
        self.resVar.set('Residue')
        self.tcVar.set(tcLabels[0])
    def __addConstCommand(self):
        res, fc, tc = self.resVar.get(), self.fcVar.get(), self.tcVar.get()
        dist = self.distEnt.get().strip()
        constraint = [res+fc, tc, dist]
        message = False
        if fc == 'Chain' or res == 'Residue' or not dist or not dist.isdigit():
            message = 'One or more of the options was filled out incorrectly.'
        elif constraint in self.parent.docker.constraints:
            message = 'That constraint has already been added.'
        if message:
            util.popupError('Problem Saving Constraint', message)
            return
        self.parent.docker.constraints.append(constraint)
        self.__fillConstTree()
    def __rmvConstCommand(self):
        sele = self.constTree.selection()
        if not sele: return
        self.parent.docker.constraints.pop(self.constTree.index(sele[0]))
        self.constTree.delete(sele)
    def __closeConstWindow(self):
        self.parent._updateConstraintsView()
        self.destroy()
