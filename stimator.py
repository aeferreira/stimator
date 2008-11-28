#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2008 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""
stimatorVersion = "0.70"
stimatorDate = "31 July 2008"

import sys
import os
import os.path
import re
import thread
import time
import DESolver
from numpy import *
import wx
import wx.lib.newevent
from   stimatorwidgts import SDLeditor, ParamGrid, TCGrid, readTCinfo
import modelparser
import resultsframe
from scipy import integrate

ABOUT_TEXT = __doc__ + "\n\nVersion %s, %s" % (stimatorVersion, stimatorDate)

demoText = """\
Write your model here...

"""

##------------- Globals for DeOdeSolver
timecoursedata = []
m_Parameters = None
NT = None
ratebytecode = None
nrates = None
v = None

##------------- Computing thread class

# NewEvent objects and a EVT binder functions
(UpdateGenerationEvent, EVT_UPDATE_GENERATION) = wx.lib.newevent.NewEvent()
(EndComputationEvent, EVT_END_COMPUTATION) = wx.lib.newevent.NewEvent()


class DeOdeSolver(DESolver.DESolver):

    def calcDerivs(self, variables, t):
        global v, nrates, ratebytecode, NT, m_Parameters
        for i in nrates:
            v[i] = eval(ratebytecode[i])
        return dot(v,NT)
        
    def setup(self, parser, calcThread):
        global m_Parameters, timecoursedata, v, nrates, ratebytecode, NT
        self.parser = parser
        self.calcThread = calcThread
        
        # cutoffEnergy is 0.1% of deviation from data
        self.cutoffEnergy =  1.0e-6*sum([nansum(abs(tc[:,1:])) for tc in timecoursedata])
        
        # scale times to maximum time in data
        scale = float(max([ (tc[-1,0]-tc[0,0]) for tc in timecoursedata]))
        
        # compute stoichiometry matrix and transpose
        N = zeros((len(parser.variables),len(parser.rates)), dtype=float)
        for m, srow in enumerate(parser.stoichmatrixrows):
            for i,k in enumerate(parser.rates):
                if srow.has_key(k['name']):
                    N[m,i] = scale *srow[k['name']]
        NT = N.transpose()

        #compile rate laws
        ratebytecode = [compile(parser.rateCalcString(k['rate']), 'bof.log','eval') for k in parser.rates]
        
        # other globals
        nrates = range(len(ratebytecode))
        v = empty(len(parser.rates))

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

        self.timecourse_scores = empty(len(timecoursedata))
        
        # dump messages to files instead of standard output streams
        #~ self.msgfile = open("msgs.txt", 'w')
        #~ self.msgfileerr = open("msgserr.txt", 'w')
        #~ self.oldstdout = sys.stdout
        #~ self.oldstderr = sys.stderr
        #~ sys.stdout = self.msgfile
        #~ sys.stderr = self.msgfileerr
        
    
    def externalEnergyFunction(self,trial):
        global m_Parameters, timecoursedata
        #~ if (trial > self.maxInitialValue).any():
            #~ return 1.0E300
        #~ if (trial < self.minInitialValue).any():
            #~ return 1.0E300
        for par in range(self.parameterCount):
            if trial[par] > self.maxInitialValue[par] or trial[par] < self.minInitialValue[par]:
                return 1.0E300
        
        m_Parameters = trial
        salg=integrate._odepack.odeint

        for i in range(len(timecoursedata)):
            y0 = copy(self.X0[i])
            t  = copy(self.times[i])

#            Y, infodict = integrate.odeint(self.calcDerivs, y0, t, full_output=True, printmessg=False)
            output = salg(self.calcDerivs, y0, t, (), None, 0, -1, -1, 0, None, 
                            None, None, 0.0, 0.0, 0.0, 0, 0, 0, 12, 5)
            if output[-1] < 0: return (1.0E300)
            #~ if infodict['message'] != 'Integration successful.':
                #~ return (1.0E300)
            Y = output[0]
            S = (Y- timecoursedata[i][:, 1:])**2
            score = nansum(S)
            self.timecourse_scores[i]=score
        
        gscore = self.timecourse_scores.sum()
        return gscore

    def reportGeneration (self):
        evt = UpdateGenerationEvent(generation = self.generation, energy = float(self.bestEnergy))
        wx.PostEvent(self.calcThread.win, evt)
    
    def reportFinal (self):
        if self.exitCode==0: outCode = -1 
        else: 
            outCode = self.exitCode
            self.calcThread.win.bestData = self.generateBestData() # bestData is own by main window
        #~ self.msgfile.close()
        #~ self.msgfileerr.close()
        #~ sys.stdout = self.oldstdout
        #~ sys.stderr = self.oldstderr
        evt = EndComputationEvent(exitCode = outCode)
        wx.PostEvent(self.calcThread.win, evt)

    def generateBestData (self):
        global timecoursedata
        bestData = [{'section':"PARAMETERS"}, 
                    {'section':"D.E. OPTIMIZATION"}, 
                    {'section':"TIMECOURSES"},
                    {'section':"best timecourses"}]
        bestData[0]['data'] = [(self.parser.parameters[i][0], "%g"%value) for (i,value) in enumerate(self.bestSolution)]
        bestData[0]['format'] = "%s\t%s"
        bestData[0]['header'] = None
        
        bestData[1]['data'] = [('Final Score', "%g"% self.bestEnergy),
                                ('Generations', "%d"% self.generation),
                                ('Exit by    ', "%s"% self.exitCodeStrings[self.exitCode])]
        bestData[1]['format'] = "%s\t%s"
        bestData[1]['header'] = None

        #TODO: Store initial solver parameters?

        #generate best time-courses
        bestData[2]['data'] = []
        bestData[3]['data'] = []
        besttimecoursedata = bestData[3]['data']
        for (i,data) in enumerate(timecoursedata):
            y0 = copy(self.X0[i])
            t  = self.times[i]
            Y, infodict = integrate.odeint(self.calcDerivs, y0, t, full_output=True, printmessg=False)
            besttimecoursedata.append(Y)
            if infodict['message'] != 'Integration successful.':
                bestData[2]['data'].append((self.parser.timecourseshortnames[i], self.parser.timecourseshapes[i][0], 1.0E300))
            else:
                S = (Y- data[:, 1:])**2
                score = nansum(S)
                bestData[2]['data'].append((self.parser.timecourseshortnames[i], self.parser.timecourseshapes[i][0], score))
        bestData[2]['format'] = "%s\t%d\t%g"
        bestData[2]['header'] = ['Name', 'Points', 'Score']
        
        return bestData
        

class CalcOptmThread:
    def __init__(self, win):
        self.win = win
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
                                 True)                   # use class random number methods

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
        self.parser = modelparser.StimatorParser()


    def __del__(self):
        pass

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

        # STC widgets configuration
        # ModelEditor
        self.ModelEditor.SetText(demoText)
        self.ModelEditor.EmptyUndoBuffer()
        # Set up the numbers in the margin for margin #1
        self.ModelEditor.SetMarginType(1, wx.stc.STC_MARGIN_NUMBER)
        # Reasonable value for, say, 4-5 digits using a mono font (40 pix)
        self.ModelEditor.SetMarginWidth(1, 40)

        # LogText
        self.LogText.SetText(u"")
        self.LogText.EmptyUndoBuffer()
        if wx.Platform == '__WXMSW__':
            face = 'Courier New'
            pb = 10
        else:
            face = 'Courier'
            pb = 12

        self.LogText.StyleSetSpec(wx.stc.STC_STYLE_DEFAULT, "size:%d,face:%s" % (pb, face))
        self.LogText.StyleClearAll()


    def __do_layout(self):
        sizer_main = wx.BoxSizer(wx.HORIZONTAL)
        sizer_bottom = wx.BoxSizer(wx.HORIZONTAL)
        sizer_top = wx.BoxSizer(wx.HORIZONTAL)


        sizer_bottom.Add(self.LogText, 1, wx.EXPAND, 0)
        self.bottom_pane.SetAutoLayout(True)
        self.bottom_pane.SetSizer(sizer_bottom)
        sizer_bottom.Fit(self.bottom_pane)
        sizer_bottom.SetSizeHints(self.bottom_pane)

        sizer_top.Add(self.ModelEditor, 1, wx.EXPAND, 0)
        self.top_pane.SetAutoLayout(True)
        self.top_pane.SetSizer(sizer_top)
        sizer_top.Fit(self.top_pane)
        sizer_top.SetSizeHints(self.top_pane)

        self.mainwindow.SplitHorizontally(self.top_pane, self.bottom_pane, sashPosition = 450)

        sizer_main.Add(self.mainwindow, 1, wx.EXPAND, 0)
        self.SetAutoLayout(True)
        self.SetSizer(sizer_main)
        self.Layout()

        self.mainwindow.SetSashGravity(1.0)

##------------- Init Subwindows

    def MakeStimatorWidgets(self):
        global ID_LT; ID_LT = wx.NewId()
        self.LogText = SDLeditor(self.bottom_pane, ID_LT, self)

        global ID_ME; ID_ME = wx.NewId()
        self.ModelEditor = SDLeditor(self.top_pane, ID_ME, self)

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
        if IsOpen: osflag = wx.OPEN|wx.FILE_MUST_EXIST

        fileName = None
        fileDialog = wx.FileDialog(self, "Choose a file", defaultDir, defaultFile, wildCard, osflag)
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
        fileNames = None
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
        #wx.LogMessage("Saving %s..." % self.fileName)
        if self.SaveFile(self.fileName) is not True:
            self.SaveFileError(self.fileName)
        self.ModelEditor.SetFocus()

    def OnSaveAsMenu(self, event):
        fileName = self.SelectFileDialog(False, self.GetFileDir(),self.GetFileName())
        if fileName is not None:
            self.fileName = fileName
            #wx.LogMessage("Saving %s..." % self.fileName)
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

        sepre = re.compile(r"\r\n|\n|\r")
        textlines = sepre.split(self.ModelEditor.GetText())
        
        self.parser.parse(textlines)
        if self.parser.error:
           self.IndicateError(self.parser)
           return

        if len(self.parser.timecourses) == 0 :
           self.MessageDialog("No time courses to load!\nPlease indicate some time courses with 'timecourse <filename>'", "Error")
           return
        os.chdir(self.GetFileDir())
        pathlist = [os.path.abspath(k) for k in self.parser.timecourses]
        #pathlist = [k.replace("\\","\\\\") for k in pathlist] #no longer needed?
        for filename in pathlist:
            if not os.path.exists(filename) or not os.path.isfile(filename):
                self.MessageDialog("Time course file \n%s\ndoes not exist"% filename, "Error")
                return

        self.write("-------------------------------------------------------")
        self.parser.timecourseheaders = []
        self.parser.timecoursenans = []
        timecoursedata = []
        for filename in pathlist:
            h,d = modelparser.readTimeCourseFromFile(filename)
            if d.shape == (0,0):
                self.MessageDialog("File\n%s\ndoes not contain valid time-course data"% filename, "Error")
                return
            else:
                self.write("%d time points for %d variables read from file %s" % (d.shape[0], d.shape[1]-1, filename))
                self.parser.timecourseheaders.append(h)
                timecoursedata.append(d)
        
        for i,d in enumerate(timecoursedata):
            if d.shape[1] != len(self.parser.variables)+1:
                self.MessageDialog("There are %i initial values in time course %s but model has %i variables"%(d.shape[1]-1,
                                   self.parser.timecourses[i],len(self.parser.variables)),"Error in data")
                return

            
        
        self.parser.timecourseshapes = [i.shape for i in timecoursedata]
        self.parser.timecourseshortnames = [os.path.split(filename)[1] for filename in pathlist]
        self.parser.problemname = self.GetFileName()
        
        
        self.write("-------------------------------------------------------")
        self.write("Solving %s..."%self.GetFileName())
        self.time0 = time.time()

        self.optimizerThread=CalcOptmThread(self)
        self.optimizerThread.Start(self.parser)
        
    def OnUpdateGeneration(self, evt):
        self.write("%-4d: %f" % (evt.generation, evt.energy))

    def OnEndComputation(self, evt):
        if evt.exitCode == -1:
            self.write("\nOptimization aborted by user!")
        else:
            self.write(self.optimizerThread.solver.reportFinalString())
            self.write("Optimization took %f s"% (time.time()-self.time0))
            self.PostProcessEnded()
        self.optimizerThread = None

    def PostProcessEnded(self):
        global timecoursedata
        solver = self.optimizerThread.solver        
        win = resultsframe.resultsFrame(self, -1, "Results", size=(350, 200), style = wx.DEFAULT_FRAME_STYLE)
        win.loadBestData(self.parser, self.bestData, timecoursedata)
        win.Show(True)


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
