import  os.path
import  wx
import  wx.grid as  gridlib
import  wx.stc  as  stc

#----------------------------------------------------------------------

debug = 1

#---------------------------------------------------------------------------

def readTCinfo(fullpath):
    retdic = {}
    if fullpath == "": return retdic

    retdic["fullpath"]=fullpath
    retdic["filename"]=os.path.split(fullpath)[1]
    return retdic

#---------------------------------------------------------------------------

class TCGrid(gridlib.Grid): ##, mixins.GridAutoEditMixin):
    def __init__(self, parent, Id, log):
        gridlib.Grid.__init__(self, parent, -1)
        ##mixins.GridAutoEditMixin.__init__(self)
        self.log = log
        self.moveTo = None

        self.Bind(wx.EVT_IDLE, self.OnIdle)

        # test all the events
        self.Bind(gridlib.EVT_GRID_CELL_LEFT_CLICK, self.OnCellLeftClick)
        self.Bind(gridlib.EVT_GRID_CELL_RIGHT_CLICK, self.OnCellRightClick)
        self.Bind(gridlib.EVT_GRID_CELL_LEFT_DCLICK, self.OnCellLeftDClick)
        self.Bind(gridlib.EVT_GRID_CELL_RIGHT_DCLICK, self.OnCellRightDClick)

        self.Bind(gridlib.EVT_GRID_LABEL_LEFT_CLICK, self.OnLabelLeftClick)
        self.Bind(gridlib.EVT_GRID_LABEL_RIGHT_CLICK, self.OnLabelRightClick)
        self.Bind(gridlib.EVT_GRID_LABEL_LEFT_DCLICK, self.OnLabelLeftDClick)
        self.Bind(gridlib.EVT_GRID_LABEL_RIGHT_DCLICK, self.OnLabelRightDClick)

        self.Bind(gridlib.EVT_GRID_ROW_SIZE, self.OnRowSize)
        self.Bind(gridlib.EVT_GRID_COL_SIZE, self.OnColSize)

        self.Bind(gridlib.EVT_GRID_RANGE_SELECT, self.OnRangeSelect)
        self.Bind(gridlib.EVT_GRID_CELL_CHANGE, self.OnCellChange)
        self.Bind(gridlib.EVT_GRID_SELECT_CELL, self.OnSelectCell)

        self.Bind(gridlib.EVT_GRID_EDITOR_SHOWN, self.OnEditorShown)
        self.Bind(gridlib.EVT_GRID_EDITOR_HIDDEN, self.OnEditorHidden)
        self.Bind(gridlib.EVT_GRID_EDITOR_CREATED, self.OnEditorCreated)


    def OnCellLeftClick(self, evt):
        #self.log.write("OnCellLeftClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnCellRightClick(self, evt):
        if not hasattr(self, "popupID1"):
            self.popupID1 = wx.NewId()
            self.popupID2 = wx.NewId()
            self.popupID3 = wx.NewId()

            self.Bind(wx.EVT_MENU, self.log.OnLoadTC, id=self.popupID1)
            self.Bind(wx.EVT_MENU, self.log.OnUnloadTC, id=self.popupID2)
            self.Bind(wx.EVT_MENU, self.log.OnUnloadTCAll, id=self.popupID3)

        menu = wx.Menu()
        menu.Append(self.popupID1, "Add...")
        menu.Append(self.popupID2, "Remove Selected")
        menu.Append(self.popupID3, "Remove all")

        self.PopupMenu(menu, evt.GetPosition())
        menu.Destroy()

    def OnCellLeftDClick(self, evt):
        #self.log.write("OnCellLeftDClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnCellRightDClick(self, evt):
        #self.log.write("OnCellRightDClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelLeftClick(self, evt):
        #self.log.write("OnLabelLeftClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelRightClick(self, evt):
        #self.log.write("OnLabelRightClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelLeftDClick(self, evt):
        #self.log.write("OnLabelLeftDClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelRightDClick(self, evt):
        #self.log.write("OnLabelRightDClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnRowSize(self, evt):
        #self.log.write("OnRowSize: row %d, %s\n" %
        #               (evt.GetRowOrCol(), evt.GetPosition()))
        evt.Skip()

    def OnColSize(self, evt):
        #self.log.write("OnColSize: col %d, %s\n" %
        #               (evt.GetRowOrCol(), evt.GetPosition()))
        evt.Skip()

    def OnRangeSelect(self, evt):
        #if evt.Selecting():
        #    self.log.write("OnRangeSelect: top-left %s, bottom-right %s\n" %
        #                   (evt.GetTopLeftCoords(), evt.GetBottomRightCoords()))
        evt.Skip()


    def OnCellChange(self, evt):
        #self.log.write("OnCellChange: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))

        # Show how to stay in a cell that has bad data.  We can't just
        # call SetGridCursor here since we are nested inside one so it
        # won't have any effect.  Instead, set coordinates to move to in
        # idle time.
        value = self.GetCellValue(evt.GetRow(), evt.GetCol())

        if value == 'no good':
            self.moveTo = evt.GetRow(), evt.GetCol()


    def OnIdle(self, evt):
        if self.moveTo != None:
            self.SetGridCursor(self.moveTo[0], self.moveTo[1])
            self.moveTo = None

        evt.Skip()


    def OnSelectCell(self, evt):
        #self.log.write("OnSelectCell: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))

        # Another way to stay in a cell that has a bad value...
        row = self.GetGridCursorRow()
        col = self.GetGridCursorCol()

        if self.IsCellEditControlEnabled():
            self.HideCellEditControl()
            self.DisableCellEditControl()

        value = self.GetCellValue(row, col)

        if value == 'no good 2':
            return  # cancels the cell selection

        evt.Skip()


    def OnEditorShown(self, evt):
        if evt.GetRow() == 6 and evt.GetCol() == 3 and \
           wx.MessageBox("Are you sure you wish to edit this cell?",
                        "Checking", wx.YES_NO) == wx.NO:
            evt.Veto()
            return

        #self.log.write("OnEditorShown: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()


    def OnEditorHidden(self, evt):
        if evt.GetRow() == 6 and evt.GetCol() == 3 and \
           wx.MessageBox("Are you sure you wish to  finish editing this cell?",
                        "Checking", wx.YES_NO) == wx.NO:
            evt.Veto()
            return

        #self.log.write("OnEditorHidden: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()


    def OnEditorCreated(self, evt):
        pass
        #self.log.write("OnEditorCreated: (%d, %d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetControl()))


#---------------------------------------------------------------------------

class ParamGrid(gridlib.Grid): ##, mixins.GridAutoEditMixin):
    def __init__(self, parent, Id, log):
        gridlib.Grid.__init__(self, parent, -1)
        ##mixins.GridAutoEditMixin.__init__(self)
        self.log = log
        self.moveTo = None

        self.Bind(wx.EVT_IDLE, self.OnIdle)

        # test all the events
        self.Bind(gridlib.EVT_GRID_CELL_LEFT_CLICK, self.OnCellLeftClick)
        self.Bind(gridlib.EVT_GRID_CELL_RIGHT_CLICK, self.OnCellRightClick)
        self.Bind(gridlib.EVT_GRID_CELL_LEFT_DCLICK, self.OnCellLeftDClick)
        self.Bind(gridlib.EVT_GRID_CELL_RIGHT_DCLICK, self.OnCellRightDClick)

        self.Bind(gridlib.EVT_GRID_LABEL_LEFT_CLICK, self.OnLabelLeftClick)
        self.Bind(gridlib.EVT_GRID_LABEL_RIGHT_CLICK, self.OnLabelRightClick)
        self.Bind(gridlib.EVT_GRID_LABEL_LEFT_DCLICK, self.OnLabelLeftDClick)
        self.Bind(gridlib.EVT_GRID_LABEL_RIGHT_DCLICK, self.OnLabelRightDClick)

        self.Bind(gridlib.EVT_GRID_ROW_SIZE, self.OnRowSize)
        self.Bind(gridlib.EVT_GRID_COL_SIZE, self.OnColSize)

        self.Bind(gridlib.EVT_GRID_RANGE_SELECT, self.OnRangeSelect)
        self.Bind(gridlib.EVT_GRID_CELL_CHANGE, self.OnCellChange)
        self.Bind(gridlib.EVT_GRID_SELECT_CELL, self.OnSelectCell)

        self.Bind(gridlib.EVT_GRID_EDITOR_SHOWN, self.OnEditorShown)
        self.Bind(gridlib.EVT_GRID_EDITOR_HIDDEN, self.OnEditorHidden)
        self.Bind(gridlib.EVT_GRID_EDITOR_CREATED, self.OnEditorCreated)


    def OnCellLeftClick(self, evt):
        #self.log.write("OnCellLeftClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnCellRightClick(self, evt):
        if not hasattr(self, "popupID1"):
            self.popupID1 = wx.NewId()
            self.popupID2 = wx.NewId()
            self.popupID3 = wx.NewId()
            self.popupID4 = wx.NewId()
            self.popupID5 = wx.NewId()

            self.Bind(wx.EVT_MENU, self.OnCopyMenu, id=self.popupID1)
            self.Bind(wx.EVT_MENU, self.OnClearMenu, id=self.popupID2)
            self.Bind(wx.EVT_MENU, self.OnSelectAllMenu, id=self.popupID3)
            self.Bind(wx.EVT_MENU, self.log.OnSaveResults, id=self.popupID4)
            self.Bind(wx.EVT_MENU, self.log.OnGenXLS, id=self.popupID5)

        menu = wx.Menu()
        menu.Append(self.popupID1, "Copy")
        menu.Append(self.popupID2, "Clear")
        menu.Append(self.popupID3, "Select all")
        menu.Append(self.popupID4, "Save results as")
        menu.Append(self.popupID5, "Generate xls file")

        self.PopupMenu(menu, evt.GetPosition())
        menu.Destroy()

    def OnCopyMenu(self, event):
        self.log.WriteText("Copy cells not implemented\n")

    def OnClearMenu(self, event):
        self.ClearSelection()
        self.ClearGrid()

    def OnSelectAllMenu(self, event):
        self.SelectAll()


    def OnCellLeftDClick(self, evt):
        #self.log.write("OnCellLeftDClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnCellRightDClick(self, evt):
        #self.log.write("OnCellRightDClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelLeftClick(self, evt):
        #self.log.write("OnLabelLeftClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelRightClick(self, evt):
        #self.log.write("OnLabelRightClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelLeftDClick(self, evt):
        #self.log.write("OnLabelLeftDClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnLabelRightDClick(self, evt):
        #self.log.write("OnLabelRightDClick: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()

    def OnRowSize(self, evt):
        #self.log.write("OnRowSize: row %d, %s\n" %
        #               (evt.GetRowOrCol(), evt.GetPosition()))
        evt.Skip()

    def OnColSize(self, evt):
        #self.log.write("OnColSize: col %d, %s\n" %
        #               (evt.GetRowOrCol(), evt.GetPosition()))
        evt.Skip()

    def OnRangeSelect(self, evt):
        #if evt.Selecting():
        #    self.log.write("OnRangeSelect: top-left %s, bottom-right %s\n" %
        #                   (evt.GetTopLeftCoords(), evt.GetBottomRightCoords()))
        evt.Skip()


    def OnCellChange(self, evt):
        #self.log.write("OnCellChange: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))

        # Show how to stay in a cell that has bad data.  We can't just
        # call SetGridCursor here since we are nested inside one so it
        # won't have any effect.  Instead, set coordinates to move to in
        # idle time.
        value = self.GetCellValue(evt.GetRow(), evt.GetCol())

        if value == 'no good':
            self.moveTo = evt.GetRow(), evt.GetCol()


    def OnIdle(self, evt):
        if self.moveTo != None:
            self.SetGridCursor(self.moveTo[0], self.moveTo[1])
            self.moveTo = None

        evt.Skip()


    def OnSelectCell(self, evt):
        #self.log.write("OnSelectCell: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))

        # Another way to stay in a cell that has a bad value...
        row = self.GetGridCursorRow()
        col = self.GetGridCursorCol()

        if self.IsCellEditControlEnabled():
            self.HideCellEditControl()
            self.DisableCellEditControl()

        value = self.GetCellValue(row, col)

        if value == 'no good 2':
            return  # cancels the cell selection

        evt.Skip()


    def OnEditorShown(self, evt):
        if evt.GetRow() == 6 and evt.GetCol() == 3 and \
           wx.MessageBox("Are you sure you wish to edit this cell?",
                        "Checking", wx.YES_NO) == wx.NO:
            evt.Veto()
            return

        #self.log.write("OnEditorShown: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()


    def OnEditorHidden(self, evt):
        if evt.GetRow() == 6 and evt.GetCol() == 3 and \
           wx.MessageBox("Are you sure you wish to  finish editing this cell?",
                        "Checking", wx.YES_NO) == wx.NO:
            evt.Veto()
            return

        #self.log.write("OnEditorHidden: (%d,%d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetPosition()))
        evt.Skip()


    def OnEditorCreated(self, evt):
        pass
        #self.log.write("OnEditorCreated: (%d, %d) %s\n" %
        #               (evt.GetRow(), evt.GetCol(), evt.GetControl()))



#---------------------------------------------------------------------------

# This shows how to catch the Modified event from the wx.StyledTextCtrl

class SDLeditor(stc.StyledTextCtrl):
    def __init__(self, parent, ID, log):
        stc.StyledTextCtrl.__init__(self, parent, ID)
        self.log = log

        self.Bind(stc.EVT_STC_DO_DROP, self.OnDoDrop)
        self.Bind(stc.EVT_STC_DRAG_OVER, self.OnDragOver)
        self.Bind(stc.EVT_STC_START_DRAG, self.OnStartDrag)
        #self.Bind(stc.EVT_STC_MODIFIED, self.OnModified)

        #self.Bind(wx.EVT_WINDOW_DESTROY, self.OnDestroy)

    #def OnDestroy(self, evt):
        # This is how the clipboard contents can be preserved after
        # the app has exited.
        #wx.TheClipboard.Flush()
        #evt.Skip()


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

#----------------------------------------------------------------------




