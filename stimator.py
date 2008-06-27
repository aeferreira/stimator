#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2008 Antonio Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""
stimatorVersion = "0.601"
stimatorDate = "27 June 2008"

import sys
import os
import os.path
import re
import thread
#import time
import DESolver
from numpy import *
import wx
import wx.lib.newevent
from   stimatorwidgts import SDLeditor, ParamGrid, TCGrid, readTCinfo
import modelparser
import stimator_timecourse
from scipy import integrate
import pylab as p

ABOUT_TEXT = __doc__ + "\n\nVersion %s, %s" % (stimatorVersion, stimatorDate)

demoText = """\
Write your model here...

"""

##------------- Globals for DeOdeSolver
timecoursedata = []
besttimecoursedata = []
m_Parameters = None

##------------- Computing thread class

# NewEvent objects and a EVT binder functions
(UpdateGenerationEvent, EVT_UPDATE_GENERATION) = wx.lib.newevent.NewEvent()
(EndComputationEvent, EVT_END_COMPUTATION) = wx.lib.newevent.NewEvent()


class DeOdeSolver(DESolver.DESolver):

    def setup(self, parser, calcThread):
        global m_Parameters, timecoursedata
        self.parser = parser
        self.calcThread = calcThread
        
        # cutoffEnergy is 0.1% of deviation from data
        self.cutoffEnergy =  0.000001*sum([nansum(abs(tc[:,1:])) for tc in timecoursedata])
        
        # scale times to maximum time in data
        scale = float(max([ (tc[-1,0]-tc[0,0]) for tc in timecoursedata]))
        
        # generate function calcDerivs (with scale)
        sss = parser.ODEcalcString(scale = scale)
        cc = compile(sss, 'bof.log','exec')
        exec cc
        self.calcDerivs = calcDerivs

        # store initial values and (scaled) time points
        self.X0 = []
        self.times = []
        for data in timecoursedata:
            y0 = copy(data[0, 1:]) # variables are in columns 1 to end
            self.X0.append(y0)
            
            t  = data[:, 0]        # times are in columns 0
            t0 = t[0]
            times = (t-t0)/scale+t0  # this scales time points
            self.times.append(times)
        
        # dump messages to files instead of standard output streams
        self.msgfile = open("msgs.txt", 'w')
        self.msgfileerr = open("msgserr.txt", 'w')
        self.oldstdout = sys.stdout
        self.oldstderr = sys.stderr
        sys.stdout = self.msgfile
        sys.stderr = self.msgfileerr
        
    
    def externalEnergyFunction(self,trial):
        global m_Parameters, timecoursedata
        for par in range(self.parameterCount):
            if trial[par] > self.maxInitialValue[par] or trial[par] < self.minInitialValue[par]:
                return 1.0E300
        
        m_Parameters = trial
        
        timecourse_scores = zeros(len(timecoursedata))

        for (i,data) in enumerate(timecoursedata):
            y0 = copy(self.X0[i])
            t  = self.times[i]

            Y, infodict = integrate.odeint(self.calcDerivs, y0, t, full_output=True, printmessg=False)

            if infodict['message'] != 'Integration successful.':
                return (1.0E300)
            S = (Y- data[:, 1:])**2
            score = nansum(S)
            timecourse_scores[i]=score
        
        gscore = timecourse_scores.sum()
        return gscore

    def reportGeneration (self):
        evt = UpdateGenerationEvent(generation = self.generation, energy = float(self.bestEnergy))
        wx.PostEvent(self.calcThread.win, evt)

    def reportFinal (self):
        global besttimecoursedata
        if self.exitCode==0: outCode = -1 
        else: 
            outCode = self.exitCode
            #generate best time-courses
            besttimecoursedata = []
            for i in range(len(timecoursedata)):
                y0 = copy(self.X0[i])
                t  = self.times[i]
                Y, infodict = integrate.odeint(self.calcDerivs, y0, t, full_output=True, printmessg=False)
                besttimecoursedata.append(Y)

        self.msgfile.close()
        self.msgfileerr.close()
        sys.stdout = self.oldstdout
        sys.stderr = self.oldstderr
        evt = EndComputationEvent(exitCode = outCode)
        wx.PostEvent(self.calcThread.win, evt)

class CalcOptmThread:
    def __init__(self, win):
        self.win = win
        #self.genNum = 0
        self.computationEnded = False

    def Start(self, parser):
        mins = array([k[1] for k in parser.parameters])
        maxs = array([k[2] for k in parser.parameters])
        
        self.solver = DeOdeSolver(len(parser.parameters), # number of parameters
                                 int(parser.genomesize),  # genome size
                                 int(parser.generations), # max number of generations
                                 mins, maxs,              # min and max parameter values
                                 "Best2Exp",              # DE strategy
                                 0.7, 0.6, 0.0,           # DiffScale, Crossover Prob, Cut off Energy
                                 False)                   # use class random number methods

        self.solver.setup(parser,self)
        
        self.keepGoing = self.running = True
        thread.start_new_thread(self.Run, ())

    def Stop(self):
        self.keepGoing = False

    def IsRunning(self):
        return self.running

    def Run(self):
        while self.keepGoing:
            self.solver.computeGeneration()
            if self.solver.exitCode !=0: self.Stop()

        self.solver.finalize()
        self.running = False


##------------- Fonts to be used.
if wx.Platform == '__WXMSW__':
    face1 = 'Arial'
    face2 = 'Times New Roman'
    face3 = 'Courier New'
    pb = 10
else:
    face1 = 'Helvetica'
    face2 = 'Times'
    face3 = 'Courier'
    pb = 12

##------------- Log class

class MyLog(wx.PyLog):
    def __init__(self, textCtrl, logTime=0):
        wx.PyLog.__init__(self)
        self.tc = textCtrl
        self.logTime = logTime

    def DoLogString(self, message, timeStamp):
        #print message, timeStamp
        #if self.logTime:
        #    message = time.strftime("%X", time.localtime(timeStamp)) + \
        #              ": " + message
        if self.tc:
            self.tc.AppendText(message + '\n')
            self.tc.GotoLine(self.tc.GetLineCount())


##------------- Main Window

class stimatorMainFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)

        self.mainwindow = wx.SplitterWindow(self, -1, style=wx.SP_3D|wx.SP_BORDER)

        self.top_pane = wx.Panel(self.mainwindow, -1)
        self.bottom_pane = wx.Panel(self.mainwindow, -1)

        self.bottomwindow = wx.SplitterWindow(self.bottom_pane, -1, style=wx.SP_3D|wx.SP_BORDER)
        self.bottom_left_pane = wx.Panel(self.bottomwindow, -1)
        self.bottom_right_pane = wx.Panel(self.bottomwindow, -1)

        self.resNotebook = wx.Notebook(self.bottom_right_pane, -1, style=0)

        self.InitVariables()

        self.MakeMenus()
        self.MakeToolbar()
        self.MakeStimatorWidgets()
        self.mainstatusbar = self.CreateStatusBar(1, wx.ST_SIZEGRIP)

        #self.Bind(wx.EVT_IDLE, self.OnIdle)
        self.Bind(EVT_UPDATE_GENERATION, self.OnUpdateGeneration)
        self.Bind(EVT_END_COMPUTATION, self.OnEndComputation)

        self.__set_properties()
        self.__do_layout()

        # Set the wxWindows log target to be the LogText control
        # using our own wx.Log class
        wx.Log_SetActiveTarget(MyLog(self.LogText))
        # for serious debugging
        #wx.Log_SetActiveTarget(wx.LogStderr())
        #wx.Log_SetTraceMask(wx.TraceMessages)

##------------- Initialization and __del__ functions

    def InitVariables(self):
        self.fileName = None
        self.TCpaths = []
        self.optimizerThread = None
        self.needRefreshParamsGrid = False
        #self.pid = None

    def __del__(self):
        pass
        #~ if self.optimizerThread is not None:
            #~ self.optimizerThread.Detach()
            #~ self.optimizerThread.CloseOutput()
            #~ self.optimizerThread = None

    def __set_properties(self):
        # main window configuration
        self.SetTitle("S-timator [untitled]")
        self.SetSize((1024, 768))
        self.SetBackgroundColour(wx.Colour(229, 229, 229))

        # statusbar configuration
        self.mainstatusbar.SetStatusWidths([-1])
        mainstatusbar_fields = ["S-timator %s"%(stimatorVersion)]
        for i in range(len(mainstatusbar_fields)):
            self.mainstatusbar.SetStatusText(mainstatusbar_fields[i], i)
        self.maintoolbar.Realize()

        # timecourse and parameter grids configuration
        #self.timecoursegrid.CreateGrid(5, 3)
        #self.timecoursegrid.SetSelectionMode(wx.grid.Grid.wxGridSelectRows)
        #self.timecoursegrid.SetColLabelValue(0, "File")
        #self.timecoursegrid.SetColSize(0, 100)
        #self.timecoursegrid.SetColLabelValue(1, "Comment")
        #self.timecoursegrid.SetColSize(1, 180)
        #self.timecoursegrid.SetColLabelValue(2, "Path")

        self.parametergrid.CreateGrid(5, 3)
        #self.parametergrid.EnableEditing(False)
        self.parametergrid.SetColLabelValue(0, "Name")
        self.parametergrid.SetColLabelValue(1, "Value")
        self.parametergrid.SetColLabelValue(2, "SE")

        # STC widgets configuration
        # ModelEditor
        self.ModelEditor.SetText(demoText)
        self.ModelEditor.EmptyUndoBuffer()

        self.ModelEditor.StyleSetSpec(wx.stc.STC_STYLE_DEFAULT, "size:%d,face:%s" % (pb, face3))
        self.ModelEditor.StyleClearAll()

        # line numbers in the margin
        self.ModelEditor.SetMarginType(0, wx.stc.STC_MARGIN_NUMBER)
        self.ModelEditor.SetMarginWidth(0, 22)
        self.ModelEditor.StyleSetSpec(wx.stc.STC_STYLE_LINENUMBER, "size:%d,face:%s" % (pb, face1))

        # LogText
        self.LogText.SetText(u"")
        self.LogText.EmptyUndoBuffer()

        self.LogText.StyleSetSpec(wx.stc.STC_STYLE_DEFAULT, "size:%d,face:%s" % (pb, face3))
        self.LogText.StyleClearAll()


    def __do_layout(self):
        sizer_1 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_2 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_5 = wx.BoxSizer(wx.HORIZONTAL)
        sizer_bottom = wx.BoxSizer(wx.HORIZONTAL)
        sizer_top = wx.BoxSizer(wx.HORIZONTAL)

        #self.resNotebook.AddPage(self.timecoursegrid, "Time courses")
        self.resNotebook.AddPage(self.parametergrid, "Parameters")

        sizer_5.Add(self.resNotebook, 1, wx.EXPAND, 0)
        self.bottom_right_pane.SetAutoLayout(True)
        self.bottom_right_pane.SetSizer(sizer_5)
        sizer_5.Fit(self.bottom_right_pane)
        sizer_5.SetSizeHints(self.bottom_right_pane)

        sizer_bottom.Add(self.LogText, 1, wx.EXPAND, 0)
        self.bottom_left_pane.SetAutoLayout(True)
        self.bottom_left_pane.SetSizer(sizer_bottom)
        sizer_bottom.Fit(self.bottom_left_pane)
        sizer_bottom.SetSizeHints(self.bottom_left_pane)


        self.bottomwindow.SplitVertically(self.bottom_left_pane, self.bottom_right_pane, sashPosition = 800)

        sizer_2.Add(self.bottomwindow, 1, wx.EXPAND, 0)
        self.bottom_pane.SetAutoLayout(True)
        self.bottom_pane.SetSizer(sizer_2)
        sizer_2.Fit(self.bottom_pane)
        sizer_2.SetSizeHints(self.bottom_pane)

        sizer_top.Add(self.ModelEditor, 1, wx.EXPAND, 0)
        self.top_pane.SetAutoLayout(True)
        self.top_pane.SetSizer(sizer_top)
        sizer_top.Fit(self.top_pane)
        sizer_top.SetSizeHints(self.top_pane)

        self.mainwindow.SplitHorizontally(self.top_pane, self.bottom_pane, sashPosition = 450)
        #self.mainwindow.SplitHorizontally(self.top_pane, self.LogText)

        sizer_1.Add(self.mainwindow, 1, wx.EXPAND, 0)
        self.SetAutoLayout(True)
        self.SetSizer(sizer_1)
        self.Layout()

        self.mainwindow.SetSashGravity(1.0)
        self.bottomwindow.SetSashGravity(0.5)

##------------- Init Subwindows

    def MakeStimatorWidgets(self):
        global ID_LT; ID_LT = wx.NewId()
        self.LogText = SDLeditor(self.bottom_left_pane, ID_LT, self)

        global ID_ME; ID_ME = wx.NewId()
        self.ModelEditor = SDLeditor(self.top_pane, ID_ME, self)
        #global ID_TCGRID; ID_TCGRID = wx.NewId()
        #self.timecoursegrid = TCGrid(self.resNotebook, ID_TCGRID, self)
        global ID_PARGRID; ID_PARGRID = wx.NewId()
        self.parametergrid = ParamGrid(self.resNotebook, ID_PARGRID, self)

    def MakeMenus(self):
        self.mainmenu = wx.MenuBar()
        self.AddMenus(self.mainmenu)
        self.SetMenuBar(self.mainmenu)

    def MakeToolbar(self):
        self.maintoolbar = wx.ToolBar(self, -1, style=wx.TB_HORIZONTAL|wx.TB_TEXT|wx.TB_NOICONS)
        self.SetToolBar(self.maintoolbar)
        buttonId = wx.NewId()
        self.maintoolbar.AddLabelTool(buttonId, "Compute", wx.NullBitmap, wx.NullBitmap, wx.ITEM_NORMAL, "", "")
        self.Bind(wx.EVT_TOOL, self.OnComputeButton, id=buttonId)
        buttonId = wx.NewId()
        self.maintoolbar.AddLabelTool(buttonId, "Abort", wx.NullBitmap, wx.NullBitmap, wx.ITEM_NORMAL, "", "")
        self.Bind(wx.EVT_TOOL, self.OnAbortButton, id=buttonId)

    def AddMenuItem(self, menu, itemText, itemDescription, itemHandler):
        menuId = wx.NewId()
        menu.Append(menuId, itemText, itemDescription)
        self.Bind(wx.EVT_MENU, itemHandler, id=menuId)
        return menuId

    def AddMenus(self, menu):
        # File menu
        fileMenu = wx.Menu()
        self.AddMenuItem(fileMenu, '&New\tCtrl-N', 'New File', self.OnNewMenu)
        self.AddMenuItem(fileMenu, '&Open\tCtrl-O', 'Open File', self.OnOpenMenu)
        self.AddMenuItem(fileMenu, '&Save\tCtrl-S', 'Save File', self.OnSaveMenu)
        self.AddMenuItem(fileMenu, 'Save &As\tCtrl-A', 'Save File As',self.OnSaveAsMenu)
        self.AddMenuItem(fileMenu, 'E&xit\tAlt-X', 'Exit', self.OnExitMenu)
        menu.Append(fileMenu, 'File')

        # Edit menu
        editMenu = wx.Menu()
        self.AddMenuItem(editMenu, 'Cut\tCtrl-X', 'Cut', self.OnCutSelection)
        self.AddMenuItem(editMenu, '&Copy\tCtrl-C', 'Copy', self.OnCopySelection)
        self.AddMenuItem(editMenu, 'Paste\tCtrl-V', 'Paste', self.OnPaste)
        menu.Append(editMenu, 'Edit')

        # Timecourse menu
        TCMenu = wx.Menu()
        self.AddMenuItem(TCMenu, 'Add...', 'Add time courses', self.OnLoadTC)
        #self.AddMenuItem(TCMenu, 'Remove selected', 'Remove selected time courses', self.OnUnloadTC)
        #self.AddMenuItem(TCMenu, 'Remove All', 'Remove all time courses', self.OnUnloadTCAll)
        menu.Append(TCMenu, 'Time Courses')

        # Results menu
        ResultsMenu = wx.Menu()
        self.AddMenuItem(ResultsMenu, 'Save Results As...', 'Save results file as', self.OnSaveResults)
        if wx.Platform == '__WXMSW__':
              self.AddMenuItem(ResultsMenu, 'Generate XLS ', 'Generate Excel file form results', self.OnGenXLS)
        menu.Append(ResultsMenu, 'Results')

        # Settings menu
        #fileMenu = wx.Menu()
        #menu.Append(fileMenu, 'Settings')

        # Help menu
        helpMenu = wx.Menu()
        self.AddMenuItem(helpMenu, 'About', 'About the program', self.OnAboutMenu)
        menu.Append(helpMenu, 'Help')

##------------- Write funcs

    def WriteText(self, text):
        if text[-1:] == '\n':
            text = text[:-1]
        wx.LogMessage(text)

    def write(self, txt):
        self.WriteText(txt)

##-------------- Dialogs

    def MessageDialog(self, text, title):
        messageDialog = wx.MessageDialog(self, text, title, wx.OK | wx.ICON_INFORMATION)
        messageDialog.ShowModal()
        messageDialog.Destroy()

    def OkCancelDialog(self, text, title):
        dialog = wx.MessageDialog(self, text, title, wx.OK | wx.CANCEL | wx.ICON_INFORMATION)
        result = dialog.ShowModal()
        dialog.Destroy()
        if result == wx.ID_OK:
            return True
        else:
            return False

    def SelectFileDialog(self, IsOpen=True, defaultDir=None, defaultFile=None, wildCard=None):
        if defaultDir == None:
            defaultDir = "."
        if defaultFile == None:
            defaultFile = ""
        if wildCard == None:
            wildCard = "*.*"
        osflag = wx.SAVE
        if IsOpen: osflag = wx.OPEN

        fileName = None
        fileDialog = wx.FileDialog(self, "Choose a file", defaultDir, defaultFile, wildCard, osflag|wx.FILE_MUST_EXIST)
        result = fileDialog.ShowModal()
        if result == wx.ID_OK:
            fileName = fileDialog.GetPath()
        fileDialog.Destroy()
        return fileName

    def SelectFilesDialog(self, IsOpen=True, defaultDir=None, defaultFile=None, wildCard=None):
        if defaultDir == None:
            defaultDir = "."
        if defaultFile == None:
            defaultFile = ""
        if wildCard == None:
            wildCard = "*.*"
        fileName = None
        fileDialog = wx.FileDialog(self, "Choose some files", defaultDir, defaultFile, wildCard, wx.OPEN|wx.MULTIPLE|wx.FILE_MUST_EXIST)
        result = fileDialog.ShowModal()
        if result == wx.ID_OK:
            fileNames = fileDialog.GetPaths()
        fileDialog.Destroy()
        return fileNames

    def OpenFileError(self, fileName):
        wx.LogMessage('Open file error.')
        self.MessageDialog("Error opening file '%s'!" % fileName, "Error")

    def SaveFileError(self, fileName):
        wx.LogMessage('Save file error.')
        self.MessageDialog("Error saving file '%s'!" % fileName, "Error")


##---------------- Utility functions

    def GetFileDir(self):
        if self.fileName is not None:
            return os.path.split(self.fileName)[0]
        return "."

    def GetFileName(self):
        if self.fileName is not None:
            return os.path.split(self.fileName)[1]
        return ""

    def NewFile(self):
        self.ModelEditor.SetText("")
        self.fileName = None
        self.SetTitle("S-timator [untitled]")

    def SaveFile(self, fileName):
        sucess = self.ModelEditor.SaveFile(fileName)
        if sucess:
             self.SetTitle("S-timator [%s]" % self.GetFileName())
        return sucess

    def OpenFile(self, fileName):
        sucess = self.ModelEditor.LoadFile(fileName)
        if sucess:
             self.fileName = fileName
             self.SetTitle("S-timator [%s]" % self.GetFileName())
        return sucess

##---------------- Event handlers

    def OnNewMenu(self, event):
        if self.ModelEditor.GetModify():
            if not self.OkCancelDialog("New file - abandon changes?", "New File"):
                return
        self.NewFile()
        self.ModelEditor.SetFocus()

    def OnOpenMenu(self, event):
        if self.ModelEditor.GetModify():
            if not self.OkCancelDialog("Open file - abandon changes?", "Open File"):
                return
        fileName = self.SelectFileDialog(True, self.GetFileDir())
        if fileName is not None:
            if self.OpenFile(fileName) is False:
                self.OpenFileError(fileName)
        self.ModelEditor.SetFocus()

    def OnSaveMenu(self, event):
        if self.fileName is None:
            return self.OnSaveAsMenu(event)
        wx.LogMessage("Saving %s..." % self.fileName)
        if self.SaveFile(self.fileName) is not True:
            self.SaveFileError(self.fileName)
        self.ModelEditor.SetFocus()

    def OnSaveAsMenu(self, event):
        fileName = self.SelectFileDialog(False, self.GetFileDir(),self.GetFileName())
        if fileName is not None:
            self.fileName = fileName
            wx.LogMessage("Saving %s..." % self.fileName)
            if self.SaveFile(self.fileName) is not True:
                self.SaveFileError(self.fileName)
        self.ModelEditor.SetFocus()

    def OnExitMenu(self, event):
        self.Close()

    def OnCutSelection(self, event):
        self.ModelEditor.Cut()

    def OnCopySelection(self, event):
        self.ModelEditor.Copy()

    def OnPaste(self, event):
        self.ModelEditor.Paste()

    def OnLoadTC(self, event):
        fileNames = self.SelectFilesDialog(True, self.GetFileDir())
        if len(fileNames) > 0:
            self.ModelEditor.Home()
            for name in fileNames:
                d,n = os.path.split(name)
                if d == self.GetFileDir():
                    name = n
                self.ModelEditor.AddText("timecourse %s"%(name))
                self.ModelEditor.NewLine()

    def OnSaveResults(self, event):
        self.write("'SaveResults' not implemented!")
        event.Skip()

    def OnAboutMenu(self, event):
        self.MessageDialog(ABOUT_TEXT, "About S-timator")
        pass

    def IndicateError(self, parser):
       self.MessageDialog("The model description contains errors.\nThe computation was aborted.", "Error")
       msg = "ERROR in line %d:" % (parser.errorline)
       msg = msg +"\n" +parser.errorlinetext
       caretline = [" "]*(len(parser.errorlinetext)+1)
       caretline[parser.errorstart] = "^"
       caretline[parser.errorend] = "^"
       caretline = "".join(caretline)
       msg = msg +"\n" + caretline
       msg = msg +"\n" + parser.error
       self.write(msg)

    def OnAbortButton(self, event):
        if self.optimizerThread is None:
           self.MessageDialog("S-timator is NOT performing a computation!", "Error")
           return

        self.optimizerThread.Stop()

    def OnComputeButton(self, event):
        global timecoursedata
        if self.optimizerThread is not None:
           self.MessageDialog("S-timator is performing a computation!\nPlease wait.", "Error")
           return
        self.LogText.Clear()
        self.LogText.Refresh()
        self.parametergrid.ClearGrid()
        self.parametergrid.ClearSelection()
        self.parametergrid.Refresh()

        parser = modelparser.StimatorParser()
        sepre = re.compile(r"\r\n|\n|\r")
        textlines = sepre.split(self.ModelEditor.GetText())
        
        parser.parse(textlines)
        if parser.error:
           self.IndicateError(parser)
           return

        if len(parser.timecourses) == 0 :
           self.MessageDialog("No time courses to load!\nPlease indicate some time courses with 'timecourse <filename>'", "Error")
           return
        os.chdir(self.GetFileDir())
        pathlist = [os.path.abspath(k) for k in parser.timecourses]
        #pathlist = [k.replace("\\","\\\\") for k in pathlist] #no longer needed?
        for filename in pathlist:
            if not os.path.exists(filename) or not os.path.isfile(filename):
                self.MessageDialog("Time course file \n%s\ndoes not exist"% filename, "Error")
                return

        #TODO: the following two lists should be owned by the computation object
        self.write("-------------------------------------------------------")
        self.timecourseheaders = []
        timecoursedata = []
        for filename in pathlist:
            h,d = stimator_timecourse.readTimeCourseFromFile(filename)
            if d.shape == (0,0):
                self.write("file %s does not contain valid time-course data"%(filename))
            else:
                self.write("%d time points for %d variables read from file %s" % (d.shape[0], d.shape[1], filename))
                self.timecourseheaders.append(h)
                timecoursedata.append(d)
        
        
        self.write("-------------------------------------------------------")
        self.write("Solving %s..."%self.GetFileName())

        self.optimizerThread=CalcOptmThread(self)
        self.optimizerThread.Start(parser)
        
    def OnUpdateGeneration(self, evt):
        self.write("%-4d: %f" % (evt.generation, evt.energy))

    def OnEndComputation(self, evt):
        if evt.exitCode == -1:
            self.write("\nOptimization aborted by user!")
        else:
            self.write(self.optimizerThread.solver.reportFinalString())
            self.needRefreshParamsGrid = True
            self.PostProcessEnded()
        self.optimizerThread = None

    #def OnIdle(self, evt):
        #~ if self.needRefreshParamsGrid:
           #~ self.needRefreshParamsGrid = False
           #~ self.PostProcessEnded()

    def PostProcessEnded(self):
        global timecoursedata, besttimecoursedata
        solver = self.optimizerThread.solver
        
        #write parameters to grid
        self.parametergrid.ClearGrid()
        self.parametergrid.ClearSelection()
        solver = self.optimizerThread.solver
        np = solver.parameterCount
        nr = self.parametergrid.GetNumberRows()
        if np > nr:
                self.parametergrid.AppendRows(np-nr)
        for i in range(np):
            self.parametergrid.SetCellValue(i,0,solver.parser.parameters[i][0])
            self.parametergrid.SetCellValue(i,1,"%f" % solver.bestSolution[i])
            self.parametergrid.SetCellValue(i,2,'N/A')
        self.parametergrid.AutoSizeColumns()
 
        SDLTSH, HTA = besttimecoursedata[0].T
        t = timecoursedata[0][:,0]
        f1 = p.figure()
        p.plot(t, SDLTSH, 'r-', label='SDLTSH')
        p.plot(t, HTA  , 'b-', label='HTA')
        p.grid()
        p.legend(loc='best')
        p.xlabel('time (s)')
        p.ylabel('concentrations (mM)')
        p.title('Example from time course TSH2a.txt')
        p.show()
        

    def OnGenXLS(self, event):
        os.chdir(self.GetFileDir())
        if not os.path.exists("best.dat"):
           self.MessageDialog("best.dat was not found", "Error")
           return
        try:
          import best2xls
          best2xls.genXLS("best.dat")
        except:
           self.MessageDialog("An error ocurred while generating the Excel file", "Error")
        pass


# end of class stimatorMainFrame


class stimatorApp(wx.App):
    def OnInit(self):
        wx.InitAllImageHandlers()
        stimatormain = stimatorMainFrame(None, -1, "")
        self.SetTopWindow(stimatormain)
        stimatormain.Show()
        return 1

# end of class stimatorApp

if __name__ == "__main__":
    stimator = stimatorApp(0)
    stimator.MainLoop()
