#!/usr/bin/env python
# -*- coding: latin1 -*-
"""S-timator : Time-course parameter estimation using Differential Evolution.

Copyright 2005-2010 António Ferreira
S-timator uses Python, SciPy, NumPy, matplotlib, wxPython, and wxWindows."""
stimatorVersion = "0.91"
stimatorDate = "22 Apr 2010"

import wx
import wx.aui
import sys
import os
import os.path
import time
import thread
import wx.lib.newevent
import resultsframe
import stimator.modelparser
import stimator.deode
import images


ABOUT_TEXT = __doc__ + "\n\nVersion %s, %s" % (stimatorVersion, stimatorDate)

demoText = """\
Write your model here...

"""
##------------- NewEvent objects and a EVT binder functions


ID_File_New = wx.NewId()
ID_File_Open = wx.NewId()
ID_File_Save_As = wx.NewId()

ID_CreatePerspective = wx.NewId()
ID_CopyPerspective = wx.NewId()

ID_TransparentHint = wx.NewId()
ID_VenetianBlindsHint = wx.NewId()
ID_RectangleHint = wx.NewId()
ID_NoHint = wx.NewId()
ID_HintFade = wx.NewId()
ID_AllowFloating = wx.NewId()
ID_NoVenetianFade = wx.NewId()
ID_TransparentDrag = wx.NewId()
ID_AllowActivePane = wx.NewId()
ID_NoGradient = wx.NewId()
ID_VerticalGradient = wx.NewId()
ID_HorizontalGradient = wx.NewId()

ID_Settings = wx.NewId()
ID_About = wx.NewId()
ID_FirstPerspective = ID_CreatePerspective+1000

(UpdateGenerationEvent, EVT_UPDATE_GENERATION) = wx.lib.newevent.NewEvent()
(EndComputationEvent, EVT_END_COMPUTATION) = wx.lib.newevent.NewEvent()

debug = 1

##------------- Main frame class

class MyFrame(wx.Frame):

    def __init__(self, parent, id=-1, title="S-timator [untitled]", 
                pos=wx.DefaultPosition,
                size=(1024, 768), style= wx.DEFAULT_FRAME_STYLE |
                                        wx.SUNKEN_BORDER |
                                        wx.CLIP_CHILDREN):

        wx.Frame.__init__(self, parent, id, title, pos, size, style)

        self._mgr = wx.aui.AuiManager(self)
        
        # init variables
        self._perspectives = []
        self.n = 0
        self.x = 0
        self.fileName = None
        self.optimizerThread = None
        
        self.SetIcon(images.getMondrianIcon())
        self.SetBackgroundColour(wx.Colour(229, 229, 229))

        # create menu
        mb = wx.MenuBar()

        # File menu
        file_menu = wx.Menu()

        file_menu.Append(ID_File_New, '&New\tCtrl-N', 'New File')
        file_menu.Append(ID_File_Open, '&Open\tCtrl-O', 'Open File')
        file_menu.Append(wx.ID_SAVE, '&Save\tCtrl-S', 'Save File')
        file_menu.Append(ID_File_Save_As, 'Save &As\tCtrl-A', 'Save File As')
        file_menu.Append(wx.ID_EXIT, 'E&xit\tAlt-X', 'Exit')

        # Edit menu
        edit_menu = wx.Menu()
        edit_menu.Append(wx.ID_UNDO, 'Undo\tCtrl-Z', 'Undo')
        edit_menu.Append(wx.ID_REDO, 'Redo\tCtrl-Y', 'Redo')
        edit_menu.Append(wx.ID_CUT, 'Cut\tCtrl-X', 'Cut')
        edit_menu.Append(wx.ID_COPY, '&Copy\tCtrl-C', 'Copy')
        edit_menu.Append(wx.ID_PASTE, 'Paste\tCtrl-V', 'Paste')

        options_menu = wx.Menu()
        options_menu.AppendRadioItem(ID_TransparentHint, "Transparent Hint")
        options_menu.AppendRadioItem(ID_VenetianBlindsHint, "Venetian Blinds Hint")
        options_menu.AppendRadioItem(ID_RectangleHint, "Rectangle Hint")
        options_menu.AppendRadioItem(ID_NoHint, "No Hint")
        options_menu.AppendSeparator();
        options_menu.AppendCheckItem(ID_HintFade, "Hint Fade-in")
        options_menu.AppendCheckItem(ID_AllowFloating, "Allow Floating")
        options_menu.AppendCheckItem(ID_NoVenetianFade, "Disable Venetian Blinds Hint Fade-in")
        options_menu.AppendCheckItem(ID_TransparentDrag, "Transparent Drag")
        options_menu.AppendCheckItem(ID_AllowActivePane, "Allow Active Pane")
        options_menu.AppendSeparator();
        options_menu.AppendRadioItem(ID_NoGradient, "No Caption Gradient")
        options_menu.AppendRadioItem(ID_VerticalGradient, "Vertical Caption Gradient")
        options_menu.AppendRadioItem(ID_HorizontalGradient, "Horizontal Caption Gradient")
##         options_menu.AppendSeparator();
##         options_menu.Append(ID_Settings, "Settings Pane")

        self._perspectives_menu = wx.Menu()
        self._perspectives_menu.Append(ID_CreatePerspective, "Create Perspective")
        self._perspectives_menu.Append(ID_CopyPerspective, "Copy Perspective Data To Clipboard")
        self._perspectives_menu.AppendSeparator()
        self._perspectives_menu.Append(ID_FirstPerspective+0, "Default Startup")
        self._perspectives_menu.Append(ID_FirstPerspective+1, "All Panes")

        help_menu = wx.Menu()
        help_menu.Append(ID_About, "About...")
        
        mb.Append(file_menu, "File")
        mb.Append(edit_menu, "Edit")
        mb.Append(self._perspectives_menu, "Perspectives")
        mb.Append(options_menu, "Options")
        mb.Append(help_menu, "Help")
        self.mb = mb
        
        self.SetMenuBar(mb)

        # statusbar configuration
##         self.mainstatusbar = self.CreateStatusBar(1, wx.ST_SIZEGRIP)
##         self.mainstatusbar.SetStatusWidths([-1])
##         mainstatusbar_fields = ["S-timator %s"%(stimatorVersion)]
##         for i in range(len(mainstatusbar_fields)):
##             self.mainstatusbar.SetStatusText(mainstatusbar_fields[i], i)
        self.mainstatusbar = wx.StatusBar(self, -1)
        self.mainstatusbar.SetFieldsCount(2)
        self.mainstatusbar.SetStatusText("S-timator %s"%(stimatorVersion), 0)
        self.mainstatusbar.SetStatusText("", 1)
        self.SetStatusBar(self.mainstatusbar)


        # min size for the frame itself isn't completely done.
        # see the end up FrameManager::Update() for the test
        # code. For now, just hard code a frame minimum size
        self.SetMinSize(wx.Size(400, 300))

        # create some toolbars

        tb2 = wx.ToolBar(self, -1, wx.DefaultPosition, wx.DefaultSize,
                         wx.TB_FLAT | wx.TB_NODIVIDER)
        tb2.SetToolBitmapSize(wx.Size(20,20))
        tb2_bmp1 = wx.ArtProvider_GetBitmap(wx.ART_QUESTION, wx.ART_OTHER, wx.Size(20, 20))
        tb2.AddTool(ID_File_Open, images.get_rt_openBitmap(), shortHelpString="Open")
        tb2.AddTool(wx.ID_SAVE, images.get_rt_saveBitmap(), shortHelpString="Save")
        tb2.AddSeparator()
        tb2.AddTool(wx.ID_CUT, images.get_rt_cutBitmap(), shortHelpString="Cut")
        tb2.AddTool(wx.ID_COPY, images.get_rt_copyBitmap(), shortHelpString="Copy")
        tb2.AddTool(wx.ID_PASTE, images.get_rt_pasteBitmap(), shortHelpString="Paste")
        tb2.AddSeparator()
        tb2.AddTool(wx.ID_UNDO, images.get_rt_undoBitmap(), shortHelpString="Undo")
        tb2.AddTool(wx.ID_REDO, images.get_rt_redoBitmap(), shortHelpString="Redo")
        tb2.AddSeparator()
        buttonId = wx.NewId()
        b = wx.Button(tb2, buttonId, "Compute", (20, 20), style=wx.NO_BORDER)
        tb2.AddControl(b)
        self.Bind(wx.EVT_BUTTON, self.OnComputeButton, b)
        #~ b.SetDefault()
        b.SetSize(b.GetBestSize())
        
        buttonId = wx.NewId()
        b = wx.Button(tb2, buttonId, "Example", (20, 20), style=wx.NO_BORDER|wx.BU_EXACTFIT )
        tb2.AddControl(b)
        self.Bind(wx.EVT_BUTTON, self.OnExampleButton, b)
        tb2.Realize()
        self.tb2 = tb2
       
        tb3 = wx.ToolBar(self, -1, wx.DefaultPosition, wx.DefaultSize,
                         wx.TB_FLAT | wx.TB_NODIVIDER)
        tb3.SetToolBitmapSize(wx.Size(16,16))
        tb3_bmp1 = wx.ArtProvider_GetBitmap(wx.ART_FOLDER, wx.ART_OTHER, wx.Size(16, 16))
        tb3.AddLabelTool(101, "Test", tb3_bmp1)
        tb3.AddLabelTool(101, "Test", tb3_bmp1)
        tb3.AddLabelTool(101, "Test", tb3_bmp1)
        tb3.AddLabelTool(101, "Test", tb3_bmp1)
        tb3.AddSeparator()
        tb3.AddLabelTool(101, "Test", tb3_bmp1)
        tb3.AddLabelTool(101, "Test", tb3_bmp1)
        tb3.Realize()

        # add a bunch of panes
                      
        self._mgr.AddPane(self.CreateTreeCtrl(), wx.aui.AuiPaneInfo().
                          Name("test8").Caption("Tree Pane").
                          Left().Layer(1).Position(1).CloseButton(True).MaximizeButton(True))
                      
        sz = self.GetClientSize()
        self._mgr.AddPane(self.CreateLog(), wx.aui.AuiPaneInfo().
                          Name("test10").Caption("Log Pane").
                          Bottom().Layer(0).Row(0).MinSize(wx.Size(200,sz.y*5/12)).CloseButton(True).MaximizeButton(True))

##         self._mgr.AddPane(SettingsPanel(self, self), wx.aui.AuiPaneInfo().
##                           Name("settings").Caption("Settings").
##                           Right().Layer(0).Position(0).CloseButton(True).MaximizeButton(True))

##                           Name("settings").Caption("Dock Manager Settings").
##                           Dockable(False).Float().Hide().CloseButton(True).MaximizeButton(True))
        
        self.plotpanel = resultsframe.YetAnotherPlot(self, color=[255.0]*3, size=(400, 500))
        self._mgr.AddPane(self.plotpanel, wx.aui.AuiPaneInfo().
                          Name("results").Caption("Results").
                          Bottom().Layer(0).Row(0).Position(1).CloseButton(True).MaximizeButton(True))
        #self.plotpanel.draw()

        # create  center pane
            
        self._mgr.AddPane(self.CreateEditor(), wx.aui.AuiPaneInfo().Name("model_editor").
                          CenterPane().Caption("Editor"))
        # add the toolbars to the manager
                        
        self._mgr.AddPane(tb2, wx.aui.AuiPaneInfo().
                          Name("tb2").Caption("Toolbar 2").
                          ToolbarPane().Top().Row(1).
                          LeftDockable(False).RightDockable(False))
                      
        self._mgr.AddPane(tb3, wx.aui.AuiPaneInfo().
                          Name("tb3").Caption("Toolbar 3").
                          ToolbarPane().Top().Row(1).Position(1).
                          LeftDockable(False).RightDockable(False))
    

        # make some default perspectives
        
        perspective_all = self._mgr.SavePerspective()
        
        all_panes = self._mgr.GetAllPanes()
        
        for ii in xrange(len(all_panes)):
            if not all_panes[ii].IsToolbar():
                all_panes[ii].Hide()
                
        self._mgr.GetPane("test10").Show().Bottom().Layer(0).Row(0).Position(0)
        self._mgr.GetPane("model_editor").Show()
        self._mgr.GetPane("results").Show()
        self._mgr.GetPane("tb3").Hide()

        perspective_default = self._mgr.SavePerspective()

        for ii in xrange(len(all_panes)):
            if not all_panes[ii].IsToolbar():
                all_panes[ii].Hide()

        self._mgr.GetPane("grid_content").Show()
        self._mgr.GetPane("test8").Show().Left().Layer(0).Row(0).Position(0)
        self._mgr.GetPane("test10").Show().Bottom().Layer(0).Row(0).Position(0)
        self._mgr.GetPane("model_editor").Show()

        self._perspectives.append(perspective_default)
        self._perspectives.append(perspective_all)

        self._mgr.GetPane("test8").Hide()
        self._mgr.GetPane("grid_content").Hide()
        self._mgr.GetPane("results").Hide()
        self._mgr.GetPane("tb3").Hide()

        # "commit" all changes made to FrameManager   
        self._mgr.Update()

        self.Bind(wx.EVT_ERASE_BACKGROUND, self.OnEraseBackground)
        self.Bind(wx.EVT_SIZE, self.OnSize)
        self.Bind(wx.EVT_CLOSE, self.OnClose)

        # Show How To Use The Closing Panes Event
        self.Bind(wx.aui.EVT_AUI_PANE_CLOSE, self.OnPaneClose)

        self.Bind(wx.EVT_MENU, self.OnNewMenu, id=ID_File_New)
        self.Bind(wx.EVT_MENU, self.OnOpenMenu, id=ID_File_Open)
        self.Bind(wx.EVT_MENU, self.OnSaveMenu, id=wx.ID_SAVE)
        self.Bind(wx.EVT_MENU, self.OnSaveAsMenu, id=ID_File_Save_As)

        self.Bind(wx.EVT_MENU, self.OnUndo, id=wx.ID_UNDO)
        self.Bind(wx.EVT_MENU, self.OnRedo, id=wx.ID_REDO)
        self.Bind(wx.EVT_MENU, self.OnCutSelection, id=wx.ID_CUT)
        self.Bind(wx.EVT_MENU, self.OnCopySelection, id=wx.ID_COPY)
        self.Bind(wx.EVT_MENU, self.OnPaste, id=wx.ID_PASTE)

        self.Bind(wx.EVT_MENU, self.OnCreatePerspective, id=ID_CreatePerspective)
        self.Bind(wx.EVT_MENU, self.OnCopyPerspective, id=ID_CopyPerspective)

        self.Bind(wx.EVT_MENU, self.OnManagerFlag, id=ID_AllowFloating)
        self.Bind(wx.EVT_MENU, self.OnManagerFlag, id=ID_TransparentHint)
        self.Bind(wx.EVT_MENU, self.OnManagerFlag, id=ID_VenetianBlindsHint)
        self.Bind(wx.EVT_MENU, self.OnManagerFlag, id=ID_RectangleHint)
        self.Bind(wx.EVT_MENU, self.OnManagerFlag, id=ID_NoHint)
        self.Bind(wx.EVT_MENU, self.OnManagerFlag, id=ID_HintFade)
        self.Bind(wx.EVT_MENU, self.OnManagerFlag, id=ID_NoVenetianFade)
        self.Bind(wx.EVT_MENU, self.OnManagerFlag, id=ID_TransparentDrag)
        self.Bind(wx.EVT_MENU, self.OnManagerFlag, id=ID_AllowActivePane)
        
        self.Bind(wx.EVT_MENU, self.OnGradient, id=ID_NoGradient)
        self.Bind(wx.EVT_MENU, self.OnGradient, id=ID_VerticalGradient)
        self.Bind(wx.EVT_MENU, self.OnGradient, id=ID_HorizontalGradient)
        self.Bind(wx.EVT_MENU, self.OnExitMenu, id=wx.ID_EXIT)
        self.Bind(wx.EVT_MENU, self.OnAboutMenu, id=ID_About)

        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_TransparentHint)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_VenetianBlindsHint)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_RectangleHint)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_NoHint)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_HintFade)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_AllowFloating)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_NoVenetianFade)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_TransparentDrag)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_AllowActivePane)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_NoGradient)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_VerticalGradient)
        self.Bind(wx.EVT_UPDATE_UI, self.OnUpdateUI, id=ID_HorizontalGradient)
    
        self.Bind(wx.EVT_MENU_RANGE, self.OnRestorePerspective, id=ID_FirstPerspective,
                  id2=ID_FirstPerspective+1000)

        self.Bind(EVT_UPDATE_GENERATION, self.OnUpdateGeneration)
        self.Bind(EVT_END_COMPUTATION, self.OnEndComputation)
        
        wx.Log_SetActiveTarget(MyLog(self.LogText))


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

    def updateButtons(self):
        canUndo = self.ModelEditor.CanUndo()
        canRedo = self.ModelEditor.CanRedo()
        canSave = self.ModelEditor.IsModified()
        self.tb2.EnableTool(wx.ID_UNDO, canUndo)
        self.tb2.EnableTool(wx.ID_REDO, canRedo)
        self.tb2.EnableTool(wx.ID_SAVE, canSave)
        self.mb.Enable(wx.ID_UNDO, canUndo)
        self.mb.Enable(wx.ID_REDO, canRedo)
        self.mb.Enable(wx.ID_SAVE, canSave)

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
        fileName = os.path.join(os.path.dirname(__file__),'models','glxs_hta.mdl')
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

    def OnCutSelection(self, event):
        self.ModelEditor.Cut()

    def OnCopySelection(self, event):
        self.ModelEditor.Copy()

    def OnPaste(self, event):
        self.ModelEditor.Paste()

    def OnUndo(self, event):
        self.ModelEditor.Undo()

    def OnRedo(self, event):
        self.ModelEditor.Redo()

    def OnPaneClose(self, event):

        caption = event.GetPane().caption

        if caption in ["Tree Pane", "Dock Manager Settings", "Fixed Pane"]:
            msg = "Are You Sure You Want To Close This Pane?"
            dlg = wx.MessageDialog(self, msg, "AUI Question",
                                   wx.YES_NO | wx.NO_DEFAULT | wx.ICON_QUESTION)

            if dlg.ShowModal() in [wx.ID_NO, wx.ID_CANCEL]:
                event.Veto()
            dlg.Destroy()


    def OnClose(self, event):
        self._mgr.UnInit()
        del self._mgr
        self.Destroy()

    def OnExitMenu(self, event):
        self.Close()

    def OnAbout(self, event):
        msg = "wx.aui Demo\n" + \
              "An advanced window management library for wxWidgets\n" + \
              "(c) Copyright 2005-2006, Kirix Corporation"
        dlg = wx.MessageDialog(self, msg, "About wx.aui Demo",
                               wx.OK | wx.ICON_INFORMATION)
        dlg.ShowModal()
        dlg.Destroy()        

    def OnAboutMenu(self, event):
        self.MessageDialog(ABOUT_TEXT, "About S-timator")
        pass

    def GetDockArt(self):
        return self._mgr.GetArtProvider()

    def DoUpdate(self):
        self._mgr.Update()

    def OnEraseBackground(self, event):
        event.Skip()

    def OnSize(self, event):
        event.Skip()

##     def OnSettings(self, event):
##         # show the settings pane, and float it
##         floating_pane = self._mgr.GetPane("settings").Float().Show()
##         if floating_pane.floating_pos == wx.DefaultPosition:
##             floating_pane.FloatingPosition(self.GetStartPosition())
##         self._mgr.Update()


    def OnGradient(self, event):
        gradient = 0
        if event.GetId() == ID_NoGradient:
            gradient = wx.aui.AUI_GRADIENT_NONE
        elif event.GetId() == ID_VerticalGradient:
            gradient = wx.aui.AUI_GRADIENT_VERTICAL
        elif event.GetId() == ID_HorizontalGradient:
            gradient = wx.aui.AUI_GRADIENT_HORIZONTAL
        self._mgr.GetArtProvider().SetMetric(wx.aui.AUI_DOCKART_GRADIENT_TYPE, gradient)
        self._mgr.Update()


    def OnManagerFlag(self, event):
        flag = 0
        eid = event.GetId()

        if eid in [ ID_TransparentHint, ID_VenetianBlindsHint, ID_RectangleHint, ID_NoHint ]:
            flags = self._mgr.GetFlags()
            flags &= ~wx.aui.AUI_MGR_TRANSPARENT_HINT
            flags &= ~wx.aui.AUI_MGR_VENETIAN_BLINDS_HINT
            flags &= ~wx.aui.AUI_MGR_RECTANGLE_HINT
            self._mgr.SetFlags(flags)

        if eid == ID_AllowFloating:
            flag = wx.aui.AUI_MGR_ALLOW_FLOATING
        elif eid == ID_TransparentDrag:
            flag = wx.aui.AUI_MGR_TRANSPARENT_DRAG
        elif eid == ID_HintFade:
            flag = wx.aui.AUI_MGR_HINT_FADE
        elif eid == ID_NoVenetianFade:
            flag = wx.aui.AUI_MGR_NO_VENETIAN_BLINDS_FADE
        elif eid == ID_AllowActivePane:
            flag = wx.aui.AUI_MGR_ALLOW_ACTIVE_PANE
        elif eid == ID_TransparentHint:
            flag = wx.aui.AUI_MGR_TRANSPARENT_HINT
        elif eid == ID_VenetianBlindsHint:
            flag = wx.aui.AUI_MGR_VENETIAN_BLINDS_HINT
        elif eid == ID_RectangleHint:
            flag = wx.aui.AUI_MGR_RECTANGLE_HINT
        
        self._mgr.SetFlags(self._mgr.GetFlags() ^ flag)


    def OnUpdateUI(self, event):
        flags = self._mgr.GetFlags()
        eid = event.GetId()
        
        if eid == ID_NoGradient:
            event.Check(self._mgr.GetArtProvider().GetMetric(wx.aui.AUI_DOCKART_GRADIENT_TYPE) == wx.aui.AUI_GRADIENT_NONE)

        elif eid == ID_VerticalGradient:
            event.Check(self._mgr.GetArtProvider().GetMetric(wx.aui.AUI_DOCKART_GRADIENT_TYPE) == wx.aui.AUI_GRADIENT_VERTICAL)

        elif eid == ID_HorizontalGradient:
            event.Check(self._mgr.GetArtProvider().GetMetric(wx.aui.AUI_DOCKART_GRADIENT_TYPE) == wx.aui.AUI_GRADIENT_HORIZONTAL)

        elif eid == ID_AllowFloating:
            event.Check((flags & wx.aui.AUI_MGR_ALLOW_FLOATING) != 0)

        elif eid == ID_TransparentDrag:
            event.Check((flags & wx.aui.AUI_MGR_TRANSPARENT_DRAG) != 0)

        elif eid == ID_TransparentHint:
            event.Check((flags & wx.aui.AUI_MGR_TRANSPARENT_HINT) != 0)

        elif eid == ID_VenetianBlindsHint:
            event.Check((flags & wx.aui.AUI_MGR_VENETIAN_BLINDS_HINT) != 0)

        elif eid == ID_RectangleHint:
            event.Check((flags & wx.aui.AUI_MGR_RECTANGLE_HINT) != 0)

        elif eid == ID_NoHint:
            event.Check(((wx.aui.AUI_MGR_TRANSPARENT_HINT |
                          wx.aui.AUI_MGR_VENETIAN_BLINDS_HINT |
                          wx.aui.AUI_MGR_RECTANGLE_HINT) & flags) == 0)

        elif eid == ID_HintFade:
            event.Check((flags & wx.aui.AUI_MGR_HINT_FADE) != 0);

        elif eid == ID_NoVenetianFade:
            event.Check((flags & wx.aui.AUI_MGR_NO_VENETIAN_BLINDS_FADE) != 0);

    def OnCreatePerspective(self, event):
        dlg = wx.TextEntryDialog(self, "Enter a name for the new perspective:", "AUI Test")
        dlg.SetValue(("Perspective %d")%(len(self._perspectives)+1))
        if dlg.ShowModal() != wx.ID_OK:
            return
        if len(self._perspectives) == 0:
            self._perspectives_menu.AppendSeparator()
        self._perspectives_menu.Append(ID_FirstPerspective + len(self._perspectives), dlg.GetValue())
        self._perspectives.append(self._mgr.SavePerspective())


    def OnCopyPerspective(self, event):
        s = self._mgr.SavePerspective()
        if wx.TheClipboard.Open():
            wx.TheClipboard.SetData(wx.TextDataObject(s))
            wx.TheClipboard.Close()
        
    def OnRestorePerspective(self, event):
        self._mgr.LoadPerspective(self._perspectives[event.GetId() - ID_FirstPerspective])

    def GetStartPosition(self):
        self.x = self.x + 20
        x = self.x
        pt = self.ClientToScreen(wx.Point(0, 0))
        return wx.Point(pt.x + x, pt.y + x)

    def OnCreateTree(self, event):
        self._mgr.AddPane(self.CreateTreeCtrl(), wx.aui.AuiPaneInfo().
                          Caption("Tree Control").
                          Float().FloatingPosition(self.GetStartPosition()).
                          FloatingSize(wx.Size(150, 300)).CloseButton(True).MaximizeButton(True))
        self._mgr.Update()

    def CreateEditor(self):
        global ID_ME; ID_ME = wx.NewId()
        ed = resultsframe.SDLeditor(self, ID_ME , self)
        self.ModelEditor = ed
        
        ed.SetText(demoText)
        ed.EmptyUndoBuffer()
        # Set up the numbers in the margin for margin #1
        ed.SetMarginType(1, wx.stc.STC_MARGIN_NUMBER)
        # Reasonable value for, say, 4-5 digits using a mono font (40 pix)
        ed.SetMarginWidth(1, 40)
        ed.SetSelBackground(True, 'PLUM')
        ed.SetWrapMode(True)
        ed.SetKeyWords(0, "variables find timecourse rate generations genomesize in reaction title")
        return ed

    def CreateLog(self):
        global ID_LT; ID_LT = wx.NewId()
        ed = resultsframe.SDLeditor(self, ID_LT , self, size = wx.Size(300,800))
        self.LogText = ed
        ed.SetText(u"")
        ed.EmptyUndoBuffer()
        ed.SetIndentationGuides(False)
        if wx.Platform == '__WXMSW__':
            face = 'Courier New'
            pb = 10
        else:
            face = 'Courier'
            pb = 12

        ed.StyleClearAll()
        ed.StyleSetSpec(wx.stc.STC_STYLE_DEFAULT, "size:%d,face:%s" % (pb, face))
        ed.StyleSetSpec(wx.stc.STC_P_WORD, "fore:#00007F,bold")
        ed.SetSelBackground(True, 'PLUM')
        ed.SetWrapMode(True)
        ed.SetKeyWords(0, "TIMECOURSES OPTIMIZATION PARAMETERS")
        return ed


    def CreateTreeCtrl(self):
        tree = wx.TreeCtrl(self, -1, wx.Point(0, 0), wx.Size(160, 250),
                           wx.TR_DEFAULT_STYLE | wx.NO_BORDER)
        root = tree.AddRoot("AUI Project")
        items = []

        imglist = wx.ImageList(16, 16, True, 2)
        imglist.Add(wx.ArtProvider_GetBitmap(wx.ART_FOLDER, wx.ART_OTHER, wx.Size(16,16)))
        imglist.Add(wx.ArtProvider_GetBitmap(wx.ART_NORMAL_FILE, wx.ART_OTHER, wx.Size(16,16)))
        tree.AssignImageList(imglist)

        items.append(tree.AppendItem(root, "Item 1", 0))
        items.append(tree.AppendItem(root, "Item 2", 0))
        items.append(tree.AppendItem(root, "Item 3", 0))
        items.append(tree.AppendItem(root, "Item 4", 0))
        items.append(tree.AppendItem(root, "Item 5", 0))

        for ii in xrange(len(items)):
            id = items[ii]
            tree.AppendItem(id, "Subitem 1", 1)
            tree.AppendItem(id, "Subitem 2", 1)
            tree.AppendItem(id, "Subitem 3", 1)
            tree.AppendItem(id, "Subitem 4", 1)
            tree.AppendItem(id, "Subitem 5", 1)
        tree.Expand(root)
        return tree

    def IndicateError(self, expt):
        #self.MessageDialog("The model description contains errors.\nThe computation was aborted.", "Error")
        if expt.physloc.nstartline == expt.physloc.nendline:
            locmsg = "Error in line %d of model definition" % (expt.physloc.nendline+1)
        else:
            locmsg = "Error in lines %d-%d of model definition" % (expt.physloc.nstartline+1,expt.physloc.nendline+1)
        self.write(locmsg)
        self.write(str(expt))
        self.ModelEditor.SetSelection(expt.physloc.start, expt.physloc.end)

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
        self._mgr.GetPane("test10").Show()
        self._mgr.GetPane("results").Hide()
        self._mgr.Update()

        textlines = [self.ModelEditor.GetLine(i) for i in range(self.ModelEditor.GetLineCount())]
        
        oldout = sys.stdout #parser may need to print messages
        sys.stdout = self
        try:
            self.model = stimator.modelparser.read_model(textlines)
            self.tc = self.model.getData('timecourses')
            self.optSettings = self.model.getData('optSettings')
        except stimator.modelparser.StimatorParserError, expt:
                self.IndicateError(expt)
                sys.stdout = oldout
                return

        ntcread = self.tc.loadTimeCourses (self.GetFileDir(), names = self.tc.defaultnames, verbose=True)
        if ntcread == 0:
           sys.stdout = oldout
           return
        
        sys.stdout = oldout

        #os.chdir(self.GetFileDir())
        
        if self.model.getData('title') == "":
           self.model.setData('title', self.GetFileName())

        self.write("-------------------------------------------------------")
        self.write("Solving %s..."%self.model.getData('title'))
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
            self.PostProcessEnded()
        self.optimizerThread = None

    def PostProcessEnded(self):
        solver = self.optimizerThread.solver        
        # generate report
        reportText = solver.reportResults()
##         reportText = "\n"
##         sections = [solver.optimum[s] for s in ['parameters', 'optimization', 'timecourses']]
##         for section in sections:
##             reportText += "%-20s --------------------------------\n" % section['name'].upper()
##             if section['header']:
##                 reportText += '\t'.join(section['header'])+'\n'
##             reportText += "\n".join([section['format'] % i for i in section['data']])
##             reportText += '\n\n'
        self.write(reportText)
        self.plotpanel.model = self.model
        self.plotpanel.solver = solver
        self.plotpanel.draw()

        self._mgr.GetPane("results").Show()
        self._mgr.Update()


    def generationTick(self, generation, energy):
        evt = UpdateGenerationEvent(generation = generation, energy = energy)
        wx.PostEvent(self, evt)
    
    def finalTick(self, exitCode):
        evt = EndComputationEvent(exitCode = exitCode)
        wx.PostEvent(self, evt)


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

##------------- Optimization thread class

class CalcOptmThread:

    def Start(self, model, optSettings, timecoursecollection, aGenerationTicker, anEndComputationTicker):
        self.solver = stimator.deode.DeODESolver(model,optSettings, timecoursecollection, None, aGenerationTicker, anEndComputationTicker)
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

app = wx.App()
frame = MyFrame(None)
frame.Show()
app.MainLoop()