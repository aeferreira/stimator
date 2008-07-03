#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

import sys
import os
import os.path
import wx
from   stimatorwidgts import SDLeditor
import modelparser
from matplotlib.numerix import arange, sin, pi, cos

import matplotlib

# uncomment the following to use wx rather than wxagg
#matplotlib.use('WX')
#from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas

# comment out the following to use wx rather than wxagg
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import _load_bitmap
from matplotlib.figure import Figure
from matplotlib.numerix.mlab import rand


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


##------------- Results Frame

class resultsFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)

        self.mainwindow = wx.SplitterWindow(self, -1, style=wx.SP_3D|wx.SP_BORDER)

        self.top_pane = wx.Panel(self.mainwindow, -1)
        self.bottom_pane = wx.Panel(self.mainwindow, -1)

        self.InitVariables()
        self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)

        #self.MakeMenus()
        #self.MakeToolbar()
        global ID_RT; ID_RT = wx.NewId()
        self.ReportEditor = SDLeditor(self.top_pane, ID_RT, self)
        self.mainstatusbar = self.CreateStatusBar(1, wx.ST_SIZEGRIP)
        
        
        self.figure = Figure(figsize=(4,4), dpi=100)
        self.axes = self.figure.add_subplot(211)
        t = arange(0.0,3.0,0.01)
        s = sin(2*pi*t)
        
        self.axes.plot(t,s)

        self.axes = self.figure.add_subplot(212)
        t = arange(0.0,3.0,0.01)
        s = cos(2*pi*t)
        
        self.axes.plot(t,s)
        self.canvas = FigureCanvas(self.bottom_pane, -1, self.figure)

        wx.EVT_PAINT(self, self.OnPaint)        

        self.__set_properties()
        self.__do_layout()

    def OnPaint(self, event):
        self.canvas.draw()
        event.Skip()


##------------- Initialization and __del__ functions

    def InitVariables(self):
        pass

    def __del__(self):
        pass

    def loadBestData(self, parser, bestData, besttimecoursedata):
        """Main initialization function.
        
        Should be called after __init__() but before Show()."""
        nameOfProblem = '...'
        self.SetTitle("Results for %s" % nameOfProblem)

        reportText = ""
        # generate report Text
        for section in bestData:
            reportText += "%-20s --------------------------------\n" % section['section']
            if section['header']:
                reportText += '\t'.join(section['header'])+'\n'
            reportText += "\n".join([section['format'] % i for i in section['data']])
            reportText += '\n\n'
        
        self.ReportEditor.SetText(reportText)
        self.ReportEditor.EmptyUndoBuffer()
       
        
    def __set_properties(self):
        # main window configuration
        self.SetTitle("Results")
        self.SetSize((1024, 768))
        self.SetBackgroundColour(wx.Colour(229, 229, 229))

        # statusbar configuration
        self.mainstatusbar.SetStatusWidths([-1])
        mainstatusbar_fields = ["Results frame"]
        for i in range(len(mainstatusbar_fields)):
            self.mainstatusbar.SetStatusText(mainstatusbar_fields[i], i)
        #self.maintoolbar.Realize()


    def __do_layout(self):
        sizer_main = wx.BoxSizer(wx.HORIZONTAL)
        sizer_bottom = wx.BoxSizer(wx.HORIZONTAL)
        sizer_top = wx.BoxSizer(wx.HORIZONTAL)

        sizer_bottom.Add(self.canvas, 1, wx.TOP | wx.LEFT | wx.EXPAND)
        self.bottom_pane.SetAutoLayout(True)
        self.bottom_pane.SetSizer(sizer_bottom)
        sizer_bottom.Fit(self.bottom_pane)
        sizer_bottom.SetSizeHints(self.bottom_pane)


        sizer_top.Add(self.ReportEditor, 1, wx.EXPAND, 0)
        self.top_pane.SetAutoLayout(True)
        self.top_pane.SetSizer(sizer_top)
        sizer_top.Fit(self.top_pane)
        sizer_top.SetSizeHints(self.top_pane)

        self.mainwindow.SplitVertically(self.top_pane, self.bottom_pane, sashPosition = -400)

        sizer_main.Add(self.mainwindow, 1, wx.EXPAND, 0)
        self.SetAutoLayout(True)
        self.SetSizer(sizer_main)

        self.Layout()

        self.mainwindow.SetSashGravity(1.0)

##------------- Init Subwindows


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

        # Results menu
        ResultsMenu = wx.Menu()
        self.AddMenuItem(ResultsMenu, 'Save Results As...', 'Save results file as', self.OnSaveResults)
        if wx.Platform == '__WXMSW__':
              self.AddMenuItem(ResultsMenu, 'Generate XLS ', 'Generate Excel file form results', self.OnGenXLS)
        menu.Append(ResultsMenu, 'Results')

        # Help menu
        helpMenu = wx.Menu()
        self.AddMenuItem(helpMenu, 'About', 'About the program', self.OnAboutMenu)
        menu.Append(helpMenu, 'Help')



##---------------- Event handlers

    def OnCloseWindow(self, event):
        self.Destroy()

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

    def OnSaveResults(self, event):
        self.write("'SaveResults' not implemented!")
        event.Skip()

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


# end of class resultsFrame


