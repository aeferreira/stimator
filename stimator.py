#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2009 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""
stimatorVersion = "0.86"
stimatorDate = "17 Jun 2009"

import sys
import os
import os.path
import re
import time
import thread
import wx
import wx.lib.newevent
import wx.stc  as  stc
import resultsframe
import modelparser
import deode

ABOUT_TEXT = __doc__ + "\n\nVersion %s, %s" % (stimatorVersion, stimatorDate)

demoText = """\
Write your model here...

"""

##------------- NewEvent objects and a EVT binder functions

(UpdateGenerationEvent, EVT_UPDATE_GENERATION) = wx.lib.newevent.NewEvent()
(EndComputationEvent, EVT_END_COMPUTATION) = wx.lib.newevent.NewEvent()

debug = 1

##------------- Fonts to be used.
if wx.Platform == '__WXMSW__':
    faces = { 'times': 'Times New Roman',
              'mono' : 'Courier New',
              'helv' : 'Arial',
              'other': 'Comic Sans MS',
              'size' : 10,
              'size2': 8,
             }
else:
    faces = { 'times': 'Times',
              'mono' : 'Courier',
              'helv' : 'Helvetica',
              'other': 'new century schoolbook',
              'size' : 12,
              'size2': 10,
             }

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

##------------- Editor class

class SDLeditor(stc.StyledTextCtrl):
    def __init__(self, parent, ID, log):
        stc.StyledTextCtrl.__init__(self, parent, ID)
        self.log = log
        self.SetLexer(stc.STC_LEX_PYTHON)
        self.SetKeyWords(0, "variables find timecourse rate generations genomesize in reaction title")

        # Highlight tab/space mixing (shouldn't be any)
        self.SetProperty("tab.timmy.whinge.level", "1")

        # Set left and right margins
        self.SetMargins(2,2)

        # Indentation and tab stuff
        self.SetIndent(4)               # Prescribed indent size for wx
        self.SetIndentationGuides(True) # Show indent guides
        self.SetBackSpaceUnIndents(True)# Backspace unindents rather than delete 1 space
        self.SetTabIndents(True)        # Tab key indents
        self.SetTabWidth(4)             # Proscribed tab size for wx
        self.SetUseTabs(False)          # Use spaces rather than tabs, or
                                        # TabTimmy will complain!    
        # White space
        self.SetViewWhiteSpace(False)   # Don't view white space

        #self.SetBufferedDraw(False)
        #self.SetViewEOL(True)
        #self.SetEOLMode(stc.STC_EOL_CRLF)
        #self.SetUseAntiAliasing(True)
        
        # No right-edge mode indicator
        self.SetEdgeMode(stc.STC_EDGE_NONE)
        #self.SetEdgeMode(stc.STC_EDGE_BACKGROUND)
        #self.SetEdgeColumn(78)

        self.Bind(stc.EVT_STC_UPDATEUI, self.OnUpdateUI)
        #self.Bind(stc.EVT_STC_MARGINCLICK, self.OnMarginClick)
        #self.Bind(wx.EVT_KEY_DOWN, self.OnKeyPressed)
        self.Bind(stc.EVT_STC_DO_DROP, self.OnDoDrop)
        self.Bind(stc.EVT_STC_DRAG_OVER, self.OnDragOver)
        self.Bind(stc.EVT_STC_START_DRAG, self.OnStartDrag)

        # Make some styles,  The lexer defines what each style is used for, we
        # just have to define what each style looks like.  This set is adapted from
        # Scintilla sample property files.

        # Global default styles for all languages
        self.StyleSetSpec(stc.STC_STYLE_DEFAULT,     "face:%(mono)s,size:%(size)d" % faces)
        self.StyleClearAll()  # Reset all to be like the default

        # Line numbers in margin
        self.StyleSetSpec(wx.stc.STC_STYLE_LINENUMBER,'fore:#000000,back:#99A9C2')    
        # Highlighted brace
        self.StyleSetSpec(wx.stc.STC_STYLE_BRACELIGHT,'fore:#00009D,back:#FFFF00')
        # Unmatched brace
        self.StyleSetSpec(wx.stc.STC_STYLE_BRACEBAD,'fore:#00009D,back:#FF0000')
        # Indentation guide
        self.StyleSetSpec(wx.stc.STC_STYLE_INDENTGUIDE, "fore:#CDCDCD")

        self.StyleSetSpec(stc.STC_STYLE_DEFAULT,     "face:%(mono)s,size:%(size)d" % faces)
        self.StyleSetSpec(stc.STC_STYLE_CONTROLCHAR, "face:%(other)s" % faces)

        # Python styles
        # Default 
        self.StyleSetSpec(stc.STC_P_DEFAULT, "fore:#000000,face:%(mono)s,size:%(size)d" % faces)
        # Comments
        self.StyleSetSpec(stc.STC_P_COMMENTLINE, "fore:#7F7F7F,face:%(helv)s,size:%(size)d" % faces)
        # Number
        self.StyleSetSpec(stc.STC_P_NUMBER, "fore:#007F7F,size:%(size)d" % faces)
        # String
        self.StyleSetSpec(stc.STC_P_STRING, "fore:#7F007F,face:%(mono)s,size:%(size)d" % faces)
        # Single quoted string
        self.StyleSetSpec(stc.STC_P_CHARACTER, "fore:#7F007F,face:%(mono)s,size:%(size)d" % faces)
        # Keyword
        self.StyleSetSpec(stc.STC_P_WORD, "fore:#00007F,bold,size:%(size)d" % faces)
        # Triple quotes
        self.StyleSetSpec(stc.STC_P_TRIPLE, "fore:#7F0000,size:%(size)d" % faces)
        # Triple double quotes
        self.StyleSetSpec(stc.STC_P_TRIPLEDOUBLE, "fore:#7F0000,size:%(size)d" % faces)
        # Class name definition
        self.StyleSetSpec(stc.STC_P_CLASSNAME, "fore:#0000FF,bold,underline,size:%(size)d" % faces)
        # Function or method name definition
        self.StyleSetSpec(stc.STC_P_DEFNAME, "fore:#007F7F,bold,size:%(size)d" % faces)
        # Operators
        self.StyleSetSpec(stc.STC_P_OPERATOR, "bold,size:%(size)d" % faces)
        # Identifiers
        self.StyleSetSpec(stc.STC_P_IDENTIFIER, "fore:#000000,face:%(mono)s,size:%(size)d" % faces)
        # Comment-blocks
        self.StyleSetSpec(stc.STC_P_COMMENTBLOCK, "fore:#7F7F7F,size:%(size)d" % faces)
        # End of line where string is not closed
        self.StyleSetSpec(stc.STC_P_STRINGEOL, "fore:#000000,face:%(mono)s,back:#E0C0E0,eol,size:%(size)d" % faces)

        self.SetCaretForeground("BLUE")


    #def OnDestroy(self, evt):
        # This is how the clipboard contents can be preserved after
        # the app has exited.
        #wx.TheClipboard.Flush()
        #evt.Skip()

    def IsModified(self):
        return self.GetModify()

    def Clear(self):
        self.ClearAll()

    def SetInsertionPoint(self, pos):
        self.SetCurrentPos(pos)
        self.SetAnchor(pos)

    def ShowPosition(self, pos):
        line = self.LineFromPosition(pos)
        #self.EnsureVisible(line)
        self.GotoLine(line)

    def GetLastPosition(self):
        return self.GetLength()

    def GetPositionFromLine(self, line):
        return self.PositionFromLine(line)

    def GetRange(self, start, end):
        return self.GetTextRange(start, end)

    def GetSelection(self):
        return self.GetAnchor(), self.GetCurrentPos()

    def SetSelection(self, start, end):
        self.SetSelectionStart(start)
        self.SetSelectionEnd(end)

    def SelectLine(self, line):
        start = self.PositionFromLine(line)
        end = self.GetLineEndPosition(line)
        self.SetSelection(start, end)

    def OnUpdateUI(self, evt):
        # check for matching braces
        braceAtCaret = -1
        braceOpposite = -1
        charBefore = None
        caretPos = self.GetCurrentPos()

        if caretPos > 0:
            charBefore = self.GetCharAt(caretPos - 1)
            styleBefore = self.GetStyleAt(caretPos - 1)

        # check before
        if charBefore and chr(charBefore) in "[]{}()" and styleBefore == stc.STC_P_OPERATOR:
            braceAtCaret = caretPos - 1

        # check after
        if braceAtCaret < 0:
            charAfter = self.GetCharAt(caretPos)
            styleAfter = self.GetStyleAt(caretPos)

            if charAfter and chr(charAfter) in "[]{}()" and styleAfter == stc.STC_P_OPERATOR:
                braceAtCaret = caretPos

        if braceAtCaret >= 0:
            braceOpposite = self.BraceMatch(braceAtCaret)

        if braceAtCaret != -1  and braceOpposite == -1:
            self.BraceBadLight(braceAtCaret)
        else:
            self.BraceHighlight(braceAtCaret, braceOpposite)
            #pt = self.PointFromPosition(braceOpposite)
            #self.Refresh(True, wxRect(pt.x, pt.y, 5,5))
            #print pt
            #self.Refresh(False)


    def OnStartDrag(self, evt):
        #self.log.write("OnStartDrag: %d, %s\n"
        #               % (evt.GetDragAllowMove(), evt.GetDragText()))

        if debug and evt.GetPosition() < 250:
            evt.SetDragAllowMove(False)     # you can prevent moving of text (only copy)
            evt.SetDragText("DRAGGED TEXT") # you can change what is dragged
            #evt.SetDragText("")             # or prevent the drag with empty text


    def OnDragOver(self, evt):
        #self.log.write(
        #    "OnDragOver: x,y=(%d, %d)  pos: %d  DragResult: %d\n"
        #    % (evt.GetX(), evt.GetY(), evt.GetPosition(), evt.GetDragResult())
        #    )

        if debug and evt.GetPosition() < 250:
            evt.SetDragResult(wx.DragNone)   # prevent dropping at the beginning of the buffer


    def OnDoDrop(self, evt):
        #self.log.write("OnDoDrop: x,y=(%d, %d)  pos: %d  DragResult: %d\n"
        #               "\ttext: %s\n"
        #               % (evt.GetX(), evt.GetY(), evt.GetPosition(), evt.GetDragResult(),
        #                  evt.GetDragText()))

        if debug and evt.GetPosition() < 500:
            evt.SetDragText("DROPPED TEXT")  # Can change text if needed
            #evt.SetDragResult(wx.DragNone)  # Can also change the drag operation, but it
                                             # is probably better to do it in OnDragOver so
                                             # there is visual feedback

            #evt.SetPosition(25)             # Can also change position, but I'm not sure why
                                             # you would want to...

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
        self.optimizerThread = None

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
        self.LogText.SetIndentationGuides(False)
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
        self.maintoolbar = wx.ToolBar(self, -1, style=wx.TB_HORIZONTAL|wx.TB_TEXT|wx.TB_NOICONS|wx.TB_FLAT)
        self.SetToolBar(self.maintoolbar)
        buttonId = wx.NewId()
        self.maintoolbar.AddLabelTool(buttonId, "Compute", wx.NullBitmap, wx.NullBitmap, wx.ITEM_NORMAL, "", "")
        self.Bind(wx.EVT_TOOL, self.OnComputeButton, id=buttonId)
        buttonId = wx.NewId()
        self.maintoolbar.AddLabelTool(buttonId, "Abort", wx.NullBitmap, wx.NullBitmap, wx.ITEM_NORMAL, "", "")
        self.Bind(wx.EVT_TOOL, self.OnAbortButton, id=buttonId)
        buttonId = wx.NewId()
        self.maintoolbar.AddLabelTool(buttonId, "Example", wx.NullBitmap, wx.NullBitmap, wx.ITEM_NORMAL, "", "")
        self.Bind(wx.EVT_TOOL, self.OnExampleButton, id=buttonId)

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

    def OnExampleButton(self, event):
        if self.ModelEditor.GetModify():
            if not self.OkCancelDialog("Open file - abandon changes?", "Open File"):
                return
        fileName = os.path.join(os.path.dirname(__file__),'examples','glxs_hta.mdl')
        if not os.path.exists(fileName) or not os.path.isfile(fileName):
            self.MessageDialog("File \n%s\ndoes not exist"% fileName, "Error")
            return
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

    def IndicateError(self, error, errorloc):
       #self.MessageDialog("The model description contains errors.\nThe computation was aborted.", "Error")
       if errorloc['line'] != -1:
           msg = "ERROR in line %d:" % (errorloc['line']+1)
           msg = msg +"\n" +errorloc['linetext']
           caretline = [" "]*(len(errorloc['linetext'])+1)
           caretline[errorloc['start']] = "^"
           caretline[errorloc['end']] = "^"
           caretline = "".join(caretline)
           msg = msg +"\n" + caretline
       else:
           msg = ""
       msg = msg +"\n" + error
       self.write(msg)

    def OnAbortButton(self, event):
        if self.optimizerThread is None:
           self.MessageDialog("S-timator is NOT performing a computation!", "Error")
           return

        self.optimizerThread.Stop()

    def OnComputeButton(self, event):
        if self.optimizerThread is not None:
           self.MessageDialog("S-timator is performing a computation!\nPlease wait.", "Error")
           return
        self.LogText.Clear()
        self.LogText.Refresh()

        sepre = re.compile(r"\r\n|\n|\r")
        textlines = sepre.split(self.ModelEditor.GetText())
        parser = modelparser.StimatorParser()
        
        oldout = sys.stdout #parser may need to print messages
        sys.stdout = self

        parser.parse(textlines)
        if parser.error:
           self.IndicateError(parser.error, parser.errorLoc)
           sys.stdout = oldout
           return

        self.model       = parser.model
        self.optSettings = parser.optSettings
        self.tc          = parser.tc

        ntcread = self.tc.loadTimeCourses (self.GetFileDir())
        if ntcread == 0:
           #self.IndicateError(parser.error, parser.errorLoc)
           sys.stdout = oldout
           return
        
        sys.stdout = oldout

        os.chdir(self.GetFileDir())
        
        if self.model.title == "":
            self.model.title   = self.GetFileName()

        self.write("-------------------------------------------------------")
        self.write("Solving %s..."%self.model.title)
        self.time0 = time.time()

        self.optimizerThread=CalcOptmThread()
        self.optimizerThread.Start(self.model, self.optSettings, self.tc, self.generationTick, self.finalTick)
        
    def OnUpdateGeneration(self, evt):
        self.write("%-4d: %f" % (evt.generation, evt.energy))

    def OnEndComputation(self, evt):
        if evt.exitCode == -1:
            self.write("\nOptimization aborted by user!")
        else:
            self.write(self.optimizerThread.solver.reportFinalString())
            #print >> self, "Optimization took %f s"% (time.time()-self.time0) #this works too!
            self.write("Optimization took %f s"% (time.time()-self.time0))
            #self.bestData = evt.bestData
            self.PostProcessEnded()
        self.optimizerThread = None

    def PostProcessEnded(self):
        solver = self.optimizerThread.solver        
        win = resultsframe.resultsFrame(self, -1, "Results", size=(350, 200), style = wx.DEFAULT_FRAME_STYLE)
        win.loadBestData(self.model, solver.optimum, solver.timecoursedata)
        win.Show(True)
    
    def generationTick(self, generation, energy):
        evt = UpdateGenerationEvent(generation = generation, energy = energy)
        wx.PostEvent(self, evt)
    
    def finalTick(self, exitCode):
        evt = EndComputationEvent(exitCode = exitCode)
        wx.PostEvent(self, evt)
# end of class stimatorMainFrame


class CalcOptmThread:

    def Start(self, model, optSettings, timecoursecollection, aGenerationTicker, anEndComputationTicker):
        self.solver = deode.DeODESolver(model,optSettings, timecoursecollection, aGenerationTicker, anEndComputationTicker)
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
