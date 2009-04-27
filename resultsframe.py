#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

import sys
import os
import os.path
import math
import wx
import  wx.stc  as  stc
import modelparser
from matplotlib.numerix import arange, sin, pi, cos, isnan

import matplotlib
matplotlib.interactive(True)
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure

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

##------------- Editor class

class SDLeditor(stc.StyledTextCtrl):
    def __init__(self, parent, ID, log):
        stc.StyledTextCtrl.__init__(self, parent, ID)
        self.log = log
        self.SetLexer(stc.STC_LEX_PYTHON)
        self.SetKeyWords(0, "TIMECOURSES OPTIMIZATION PARAMETERS")

        # Highlight tab/space mixing (shouldn't be any)
        self.SetProperty("tab.timmy.whinge.level", "1")

        # Set left and right margins
        self.SetMargins(2,2)

        # Indentation and tab stuff
        self.SetIndent(4)                # Proscribed indent size for wx
        self.SetIndentationGuides(False) # Show indent guides
        self.SetBackSpaceUnIndents(True) # Backspace unindents rather than delete 1 space
        self.SetTabIndents(True)         # Tab key indents
        self.SetTabWidth(4)              # Proscribed tab size for wx
        self.SetUseTabs(False)           # Use spaces rather than tabs, or
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


    def OnDestroy(self, evt):
        # This is how the clipboard contents can be preserved after
        # the app has exited.
        wx.TheClipboard.Flush()
        evt.Skip()

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

##------------- Supporting classes for a "no flicker"-wxPanel to draw matplotlib plots
#~ """
#~ A demonstration of creating a matlibplot window from within wx.
#~ A resize only causes a single redraw of the panel.
#~ The WXAgg backend is used as it is quicker.

#~ Edward Abraham, Datamine, April, 2006
#~ (works with wxPython 2.6.1, Matplotlib 0.87 and Python 2.4)
#~ """




#~ class NoRepaintCanvas(FigureCanvas):
    #~ """We subclass FigureCanvasWxAgg, overriding the _onPaint method, so that
    #~ the draw method is only called for the first two paint events. After that,
    #~ the canvas will only be redrawn when it is resized.
    #~ """
    #~ def __init__(self, *args, **kwargs):
        #~ FigureCanvas.__init__(self, *args, **kwargs)
        #~ self._drawn = 0
        #~ self._isRealized = False

    #~ def _onPaint(self, evt):
        #~ """
        #~ Called when wxPaintEvt is generated
        #~ """
        #~ if not self._isRealized:
            #~ self.realize()
        #~ if self._drawn < 2:
            #~ self.draw(repaint = False)
            #~ self._drawn += 1
        #~ self.gui_repaint(drawDC=wx.PaintDC(self))

#~ class PlotPanel(wx.Panel):
    #~ """
    #~ The PlotPanel has a Figure and a Canvas. OnSize events simply set a 
    #~ flag, and the actually redrawing of the
    #~ figure is triggered by an Idle event.
    #~ """
    #~ def __init__(self, parent, id = -1, color = None,\
        #~ dpi = None, style = wx.NO_FULL_REPAINT_ON_RESIZE, **kwargs):
        #~ wx.Panel.__init__(self, parent, id = id, style = style, **kwargs)
        #~ self.figure = Figure(None, dpi)
        #~ self.canvas = NoRepaintCanvas(self, -1, self.figure)
        #~ self.SetColor(color)
        #~ self.Bind(wx.EVT_IDLE, self._onIdle)
        #~ self.Bind(wx.EVT_SIZE, self._onSize)
        #~ self._resizeflag = True
        #~ self._SetSize()

    #~ def SetColor(self, rgbtuple):
        #~ """Set figure and canvas colours to be the same"""
        #~ if not rgbtuple:
            #~ rgbtuple = wx.SystemSettings.GetColour(wx.SYS_COLOUR_BTNFACE).Get()
        #~ col = [c/255.0 for c in rgbtuple]
        #~ self.figure.set_facecolor(col)
        #~ self.figure.set_edgecolor(col)
        #~ self.canvas.SetBackgroundColour(wx.Colour(*rgbtuple))

    #~ def _onSize(self, event):
        #~ self._resizeflag = True

    #~ def _onIdle(self, evt):
        #~ if self._resizeflag:
            #~ self._resizeflag = False
            #~ self._SetSize()
            #~ self.draw()

    #~ def _SetSize(self, pixels = None):
        #~ """
        #~ This method can be called to force the Plot to be a desired size, which defaults to
        #~ the ClientSize of the panel
        #~ """
        #~ if not pixels:
            #~ pixels = self.GetClientSize()
        #~ self.canvas.SetSize(pixels)
        #~ self.figure.set_size_inches(pixels[0]/self.figure.get_dpi(),
                                    #~ pixels[1]/self.figure.get_dpi())

    #~ def draw(self):
        #~ """Where the actual drawing happens"""
        #~ pass

class PlotPanel (wx.Panel):
    """The PlotPanel has a Figure and a Canvas. OnSize events simply set a 
flag, and the actual resizing of the figure is triggered by an Idle event."""
    def __init__( self, parent, color=None, dpi=None, **kwargs ):
        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
        from matplotlib.figure import Figure

        # initialize Panel
        if 'id' not in kwargs.keys():
            kwargs['id'] = wx.ID_ANY
        if 'style' not in kwargs.keys():
            kwargs['style'] = wx.NO_FULL_REPAINT_ON_RESIZE
        wx.Panel.__init__( self, parent, **kwargs )

        # initialize matplotlib stuff
        self.figure = Figure( None, dpi )
        self.canvas = FigureCanvasWxAgg( self, -1, self.figure )
        self.SetColor( color )

        self._SetSize()
#        self.draw()

        self._resizeflag = False

        self.Bind(wx.EVT_IDLE, self._onIdle)
        self.Bind(wx.EVT_SIZE, self._onSize)

    def SetColor( self, rgbtuple=None ):
        """Set figure and canvas colours to be the same."""
        if rgbtuple is None:
            rgbtuple = wx.SystemSettings.GetColour( wx.SYS_COLOUR_BTNFACE ).Get()
        clr = [c/255. for c in rgbtuple]
        self.figure.set_facecolor( clr )
        self.figure.set_edgecolor( clr )
        self.canvas.SetBackgroundColour( wx.Colour( *rgbtuple ) )

    def _onSize( self, event ):
        self._resizeflag = True

    def _onIdle( self, evt ):
        if self._resizeflag:
            self._resizeflag = False
            self._SetSize()
            self.draw()

    def _SetSize( self ):
        pixels = tuple( self.GetClientSize() )
        #pixels = tuple( self.parent.GetClientSize() )
        self.SetSize( pixels )
        self.canvas.SetSize( pixels )
        self.figure.set_size_inches( float( pixels[0] )/self.figure.get_dpi(),
                                     float( pixels[1] )/self.figure.get_dpi() )

    def draw(self): pass # abstract, to be overridden by child classes

class DemoPlotPanel(PlotPanel):
    """An example plotting panel. The only method that needs 
    overriding is the draw method"""
    def draw(self):
        #self._SetSize()    #?????
        if not hasattr(self, 'subplot'):
            self.subplot = self.figure.add_subplot(111)
        theta = arange(0, 45*2*pi, 0.02)
        rad = (0.8*theta/(2*pi)+1)
        r = rad*(8 + sin(theta*7+rad/1.8))
        x = r*cos(theta)
        y = r*sin(theta)
        #Now draw it
        self.subplot.plot(x,y, '-b')
        #Set some plot attributes
        #self.subplot.set_title("A polar flower (%s points)" % len(x), fontsize = 12)
        self.subplot.set_title("Results for %s" % self.parser.problemname)
        self.subplot.set_xlabel("Flower is from  http://www.physics.emory.edu/~weeks/ideas/rose.html")
        self.subplot.set_xlim([-400, 400])
        self.subplot.set_ylim([-400, 400])

class BestPlotPanel(PlotPanel):
    """Plots best data."""
    #TODO: implement graph settings
    def draw(self):
        self._SetSize()    #?????
        self.figure.clear()
        #if not hasattr(self, 'tcsubplots'):
        self.tcsubplots = []
        ntc = len(self.bestData[3]['data'])
        besttimecoursedata = self.bestData[3]['data']
        timecoursedata = self.timecoursedata
        ncols = int(math.ceil(math.sqrt(ntc)))
        nrows = int(math.ceil(float(ntc)/ncols))
        #print "ncols=",ncols,'nrows=',nrows
        for i in range(ntc):
            self.tcsubplots.append(self.figure.add_subplot(nrows,ncols,i+1))

        for i in range(ntc):
            subplot = self.tcsubplots[i]
            #subplot.set_xlabel("time")
            subplot.set_title("%s (%d)%g"% self.bestData[2]['data'][i], fontsize = 12)
            x = timecoursedata[i][:,0]
            nx = len(x)
            for line in range(1, timecoursedata[i].shape[1]):
                #count NaN
                yexp = timecoursedata[i][:,line]
                nnan = len(yexp[isnan(yexp)])
                if nnan >= nx-1: continue
                ysim = besttimecoursedata[i][:,line-1]
                subplot.plot(x,yexp, '-b')
                subplot.plot(x,ysim, '-r')
            #yexp = timecoursedata[i][:,1:]
            #ysim = besttimecoursedata[i]
            #Set some plot attributes
            #self.subplot.set_title("Results for %s" % self.parser.problemname)
            #subplot.set_xlim([-400, 400])
            #subplot.set_ylim([-400, 400])




##------------- Results Frame

class resultsFrame(wx.Frame):
    def __init__(self, *args, **kwds):
        kwds["style"] = wx.DEFAULT_FRAME_STYLE
        wx.Frame.__init__(self, *args, **kwds)

        self.notebook = wx.Notebook(self, -1, style=0)

        global ID_RT; ID_RT = wx.NewId()
        self.ReportEditor = SDLeditor(self.notebook, ID_RT, self)
        
        self.plotpanel = BestPlotPanel(self.notebook, color=[255.0]*3)
        
        self.notebook.AddPage(self.ReportEditor, "Results")
        self.notebook.AddPage(self.plotpanel, "Plots") 

        self.Bind(wx.EVT_CLOSE, self.OnCloseWindow)

        #self.MakeMenus()
        #self.MakeToolbar()

        #self.SetTitle("Results")
        self.SetSize((1024, 768))
        #self.SetBackgroundColour(wx.Colour(229, 229, 229))
        

        # statusbar configuration
        self.mainstatusbar = self.CreateStatusBar(1, wx.ST_SIZEGRIP)
        self.mainstatusbar.SetStatusWidths([-1])
        mainstatusbar_fields = ["Results frame"]
        for i in range(len(mainstatusbar_fields)):
            self.mainstatusbar.SetStatusText(mainstatusbar_fields[i], i)
        #self.maintoolbar.Realize()


##------------- Initialization and __del__ functions

    def __del__(self):
        pass
    
    def loadBestData(self, parser, bestData, timecoursedata):
        """Main initialization function.
        
        Should be called after __init__() but before Show()."""

        self.plotpanel.parser = parser
        self.plotpanel.bestData = bestData
        self.plotpanel.timecoursedata = timecoursedata

        self.parser = parser
        self.bestData = bestData
        self.SetTitle("Results for %s" % parser.problemname)

        # generate report
        reportText = ""
        for section in bestData:
            if section['section'] =="best timecourses": continue
            reportText += "%-20s --------------------------------\n" % section['section']
            if section['header']:
                reportText += '\t'.join(section['header'])+'\n'
            reportText += "\n".join([section['format'] % i for i in section['data']])
            reportText += '\n\n'
        
        self.ReportEditor.SetText(reportText)
        self.ReportEditor.EmptyUndoBuffer()
               

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

    def OnSaveMenu(self, event):
        if self.fileName is None:
            return self.OnSaveAsMenu(event)
        #wx.LogMessage("Saving %s..." % self.fileName)
        if self.SaveFile(self.fileName) is not True:
            self.SaveFileError(self.fileName)
        self.ReportEditor.SetFocus()

    def OnSaveAsMenu(self, event):
        fileName = self.SelectFileDialog(False, self.GetFileDir(),self.GetFileName())
        if fileName is not None:
            self.fileName = fileName
            #wx.LogMessage("Saving %s..." % self.fileName)
            if self.SaveFile(self.fileName) is not True:
                self.SaveFileError(self.fileName)
        self.ReportEditor.SetFocus()

    def OnExitMenu(self, event):
        self.Close()

    def OnCutSelection(self, event):
        self.ReportEditor.Cut()

    def OnCopySelection(self, event):
        self.ReportEditor.Copy()

    def OnPaste(self, event):
        self.ReportEditor.Paste()

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

