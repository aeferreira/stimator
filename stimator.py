#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : A Python front-end to AGEDO.

Copyright 2005-2006 Antonio Ferreira
S-timator uses Python, wxPython, and wxWindows."""
stimatorVersion = "0.20"
stimatorDate = "7 July 2006"

import os
import os.path
import re
import wx
from   stimatorwidgts import SDLeditor, ParamGrid, TCGrid, readTCinfo
import modelparser
import stimator_drivergen

ABOUT_TEXT = __doc__ + "\n\nVersion %s, %s" % (stimatorVersion, stimatorDate)

demoText = """\
//Write your model here...

"""

# fonts to be used.
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

        self.Bind(wx.EVT_IDLE, self.OnIdle)

        # We can either derive from wx.Process and override OnTerminate
        # or we can let wx.Process send this window an event that is
        # caught in the normal way...
        self.Bind(wx.EVT_END_PROCESS, self.OnProcessEnded)

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
        self.AGEDOprocess = None
        self.needRefreshParamsGrid = False
        self.pid = None

    def __del__(self):
        if self.AGEDOprocess is not None:
            self.AGEDOprocess.Detach()
            self.AGEDOprocess.CloseOutput()
            self.AGEDOprocess = None

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
#        self.LogText = wx.TextCtrl(self.bottom_pane, ID_LT, "", style=wx.TE_MULTILINE|wx.TE_READONLY|wx.HSCROLL)
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
        self.AddMenuItem(TCMenu, 'Remove selected', 'Remove selected time courses', self.OnUnloadTC)
        self.AddMenuItem(TCMenu, 'Remove All', 'Remove all time courses', self.OnUnloadTCAll)
        menu.Append(TCMenu, 'Time Courses')

        # Results menu
        ResultsMenu = wx.Menu()
        self.AddMenuItem(ResultsMenu, 'Save Results As...', 'Save results file as', self.OnSaveResults)
        if wx.Platform == '__WXMSW__':
              self.AddMenuItem(ResultsMenu, 'Generate XLS ', 'Generate Excel file form results', self.OnGenXLS)
        menu.Append(ResultsMenu, 'Results')

        # Settings menu
        fileMenu = wx.Menu()
        menu.Append(fileMenu, 'Settings')

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

    def GetCurrentDir(self):
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

##---------------- Grid functions

    def RefreshTCGrid(self):
        self.timecoursegrid.ClearGrid()
        self.timecoursegrid.ClearSelection()
        self.timecoursegrid.EnableEditing(0)
        nr = self.timecoursegrid.GetNumberRows()
        np = len(self.TCpaths)
        if np > nr:
                self.timecoursegrid.AppendRows(np-nr)
        for i in range(np):
            info = readTCinfo(self.TCpaths[i])
            self.timecoursegrid.SetCellValue(i,2,info["fullpath"])
            self.timecoursegrid.SetCellValue(i,0,info["filename"])
        self.timecoursegrid.AutoSizeColumns()

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
        fileName = self.SelectFileDialog(True, self.GetCurrentDir())
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
        fileName = self.SelectFileDialog(False, self.GetCurrentDir(),self.GetFileName())
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
        fileNames = self.SelectFilesDialog(True, self.GetCurrentDir())
        if len(fileNames) > 0:
            for n in fileNames:
              if self.TCpaths.count(n)==0:
                 self.TCpaths.append(n)
            self.RefreshTCGrid()

    def OnUnloadTC(self, event):
        tl = self.timecoursegrid.GetSelectionBlockTopLeft()
        br = self.timecoursegrid.GetSelectionBlockBottomRight()
        toRemove=[]
        for i in range(len(tl)) :
              toRemove.extend(range(tl[i][0],br[i][0]+1))
        toRemove.reverse()
        for i in toRemove:
             del self.TCpaths[i]
        self.RefreshTCGrid()

    def OnUnloadTCAll(self, event):
        del self.TCpaths[:]
        self.RefreshTCGrid()

    def OnSaveResults(self, event):
        self.write("`SaveResults' not implemented!")
        event.Skip()

    def OnAboutMenu(self, event):
        self.MessageDialog(ABOUT_TEXT, "About S-timator")
        pass

    def OnComputeButton(self, event):
        if self.AGEDOprocess is not None:
           self.MessageDialog("Stimator is performing a computation!\nPlease wait.", "Error")
           return
        #if len(self.TCpaths) == 0 :
           #self.MessageDialog("No time courses loaded!\nPlease load some time courses with menu Time Courses -> Add...", "Error")
           #return
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
           return

        if len(parser.timecourses) == 0 :
           self.MessageDialog("No time courses to load!\nPlease indicate some time courses with 'timecourse <filename>'", "Error")
           return

        #generate C++ driver
        currDir = self.GetCurrentDir()
        os.chdir(currDir)
        templ = stimator_drivergen.stimatorDriverCPPtemplate
        templ = templ.replace("#####verbose_level#####",       str(1) )
        templ = templ.replace("#####nTCs#####",                str(len(parser.timecourses)) )
        #templ = templ.replace("#####nTCs#####",                str(len(self.TCpaths)) )
        templ = templ.replace("#####npars#####",               str(len(parser.parameters)) )
        templ = templ.replace("#####nvars#####",               str(len(parser.variables)) )
        templ = templ.replace("#####individuals#####",         str(parser.genomesize))
        templ = templ.replace("#####generations#####",         str(parser.generations))

        #pathlist = [k.replace("\\","\\\\") for k in self.TCpaths]
        #templ = templ.replace("#####filelist#####",            '"%s"'% ('",\n"'.join(pathlist)))
        #pathlist = [os.path.split(k)[1] for k in self.TCpaths]
        #templ = templ.replace("#####filenameslist#####",       '"%s"'% ('",\n"'.join(pathlist)))

        pathlist = parser.timecourses
        pathlist = [os.path.abspath(k) for k in pathlist]
        pathlist = [k.replace("\\","\\\\") for k in pathlist]
        templ = templ.replace("#####filelist#####",            '"%s"'% ('",\n"'.join(pathlist)))
        pathlist = [os.path.split(k)[1] for k in pathlist]
        templ = templ.replace("#####filenameslist#####",       '"%s"'% ('",\n"'.join(pathlist)))

        templ = templ.replace("#####parameternames#####",      '"%s"'% ('","'.join([k[0] for k in parser.parameters])))
        templ = templ.replace("#####varnames#####",            '"%s"'% ('","'.join(parser.variables)))
        templ = templ.replace("#####constraints_min#####",     ','.join([str(k[1]) for k in parser.parameters]))
        templ = templ.replace("#####constraints_max#####",     ','.join([str(k[2]) for k in parser.parameters]))
        templ = templ.replace("#####CalculateDerivatives#####", stimator_drivergen.writeAGEDOCalculateDerivatives(parser))

        self.write("writing C++ driver file...")
        out = open("stimator_driver.cpp", 'w')
        out.write(templ)
        out.close()
        self.write("OK")
        #return
        #compile C++ driver

        self.write(" ")
        self.write("compiling C++ driver file...")
        mypath = os.path.split(__file__)[0]
        cmd = "bcc32 "
        cmd = cmd + '-I"%s/agedo" ' % mypath
        cmd = cmd + '-L"%s" ' % mypath
        cmd = cmd + "-DNDEBUG -O2 -Oi -Ov -Oc -w-8027 "
        cmd = cmd + "stimator_driver.cpp AGEDOLib.lib"

        fo = os.popen(cmd)
        for ll in fo:
            self.write(ll)
        ret = fo.close()
        if ret :
            self.write("Compilation FAILED!\n")
            return
        self.write("compilation sucessful!\n")
        os.remove("stimator_driver.obj")
        os.remove("stimator_driver.tds")
        #os.remove("stimator_driver.cpp") # this should be optional
        
        #execute driver

        if os.path.exists("best.dat"):
           os.remove("best.dat")
        self.write(" ")
        self.write("Running driver executable...")
        self.AGEDOprocess = wx.Process(self)
        self.AGEDOprocess.Redirect();
        cmd = "stimator_driver.exe"
        self.pid = wx.Execute(cmd, wx.EXEC_ASYNC, self.AGEDOprocess)
        if not self.pid:
            self.write("Could not start %s"%cmd)
            self.AGEDOprocess.Destroy()
            self.AGEDOprocess = None
            return
        self.write('process started: "%s" pid: %s\n' % (cmd, self.pid))


    def OnIdle(self, evt):
        if self.AGEDOprocess is not None:
            stream = self.AGEDOprocess.GetInputStream()
            if stream.CanRead():
                text = stream.read()
                self.write(text)
            return
        if self.needRefreshParamsGrid:
           self.needRefreshParamsGrid = False
           self.PostProcessEnded()

    def OnProcessEnded(self, evt):
        self.write('AGEDO ended, pid:%s,  exitCode: %s\n' %
                       (evt.GetPid(), evt.GetExitCode()))

        stream = self.AGEDOprocess.GetInputStream()

        if stream.CanRead():
            text = stream.read()
            self.write(text)

        self.AGEDOprocess.Destroy()
        self.AGEDOprocess = None
        self.needRefreshParamsGrid = True

    def PostProcessEnded(self):
        self.write("Reading results...")
        os.chdir(self.GetCurrentDir())
        if not os.path.exists("best.dat"):
           self.write("Results not available!")
           return
        #if os.path.exists("generations.dat"):
           #os.remove("generations.dat")
        if os.path.exists("stimator_driver.exe"): #this should be optional
           os.remove("stimator_driver.exe")
        best = open ("best.dat", 'r')
        line = best.readline()
        while not line.startswith("Parameters"):
              line = best.readline()
        parlist = []
        line = best.readline()
        line = best.readline()
        tokens = line.split()
        while len(tokens)>0:
           parlist.append(tokens)
           line = best.readline()
           tokens = line.split()
        best.close()
        self.parametergrid.ClearGrid()
        self.parametergrid.ClearSelection()
        #self.parametergrid.EnableEditing(0)
        np = len(parlist)
        nr = self.parametergrid.GetNumberRows()
        if np > nr:
                self.parametergrid.AppendRows(np-nr)
        for i in range(np):
            self.parametergrid.SetCellValue(i,0,parlist[i][0])
            self.parametergrid.SetCellValue(i,1,parlist[i][1])
            self.parametergrid.SetCellValue(i,2,parlist[i][2])
        self.parametergrid.AutoSizeColumns()
        self.write("Results read.")
        #self.resNotebook.SetSelection(1)

    def OnAbortButton(self, event):
        if self.AGEDOprocess is None:
           self.MessageDialog("Stimator is not performing a computation!", "Error")
           return

        wx.Kill(self.pid, wx.SIGKILL)

    def OnGenXLS(self, event):
        os.chdir(self.GetCurrentDir())
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
