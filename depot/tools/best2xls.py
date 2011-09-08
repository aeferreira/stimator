#!/usr/bin/env python
# -*- coding: ISO-8859-1 -*-

#s-timator project
# Converter from best.dat to Excel
# Antonio Ferreira sep 2005

import sys
import os.path
import win32com.client

def genXLS (bestfilename):

    nseries = 2
    maxnseries = 10

    #Open Excel and template.xls
    xl = win32com.client.Dispatch("Excel.Application")
    xl.Visible = 1
    mypath = os.path.split(__file__)[0]
    filename = os.path.join(mypath, "template.xls")
    xlBook = xl.Workbooks.Open(filename)

    #eliminate extra series in "Chart 1"
    shapesname = "Chart 1"
    xl.ActiveSheet.Shapes(shapesname).Select()
    for i in range(maxnseries,nseries,-1):
        xl.ActiveChart.SeriesCollection(i).Delete()

    #Open best.dat and fill spreadsheet, row by row
    savepath = os.path.abspath(bestfilename)
    best = open (savepath, 'r')
    data = best.read()
    best.close()

    lines = data.splitlines()

    irow = 0
    processblock = False
    blockstarts = {}
    blockends = {}
    currblockname = ""

    timecoursenames = []
    workingblock =[]

    for l in lines:
        irow = irow + 1
        l = l.strip()
        if len(l) == 0:
            if not processblock : continue
            processblock = False
            blockends[currblockname] = irow - 1
            if currblockname == "" : continue
            xl.ActiveSheet.Range(xl.ActiveSheet.Cells(blockstarts[currblockname]+1, 1),
                                 xl.ActiveSheet.Cells(irow-1, len(workingblock[0] ))).Value = workingblock
            del(workingblock[:])
            currblockname == ""
        elif processblock:
            tokens = l.split("\t")
            workingblock.append(tokens)
            #xl.ActiveSheet.Range(xl.ActiveSheet.Cells(irow, 1), xl.ActiveSheet.Cells(irow, len(tokens))).Value = tokens
        elif l.startswith("Score"):
            tokens = l.split(":")
            tokens[1] = float(tokens[1])
            xl.ActiveSheet.Range(xl.ActiveSheet.Cells(irow, 1), xl.ActiveSheet.Cells(irow, len(tokens))).Value = tokens
        elif l.startswith("Parameters:"):
            xl.ActiveSheet.Cells(irow,1).Value = l
            processblock = True
            currblockname = "Parsblock"
            blockstarts[currblockname] = irow
        elif l.startswith("Time courses used:"):
            xl.ActiveSheet.Cells(irow,1).Value = l
            processblock = True
            currblockname = l
            blockstarts[currblockname] = irow
        elif l.startswith("Time course "):
            timecoursenames.append(l)
            processblock = True
            xl.ActiveSheet.Cells(irow,1).Value = l
            currblockname = l
            blockstarts[currblockname] = irow
        else: pass

    #Clone template.xls chart and adjust series pointers
    for k in range(len(timecoursenames)-1):
           xl.ActiveSheet.ChartObjects(1).Duplicate()

    k = 1
    for tc in timecoursenames:
          l,c = divmod(k-1,3)
          xl.ActiveSheet.ChartObjects(k).Left = 200 + 290 * c
          xl.ActiveSheet.ChartObjects(k).Top  = 30 + 221 * l
          cht = xl.ActiveSheet.ChartObjects(k).Chart
          cht.ChartTitle.Characters.Text = tc
          for i in range(1,nseries+1,1):
             cht.SeriesCollection(i).Name = "=Sheet1!R%dC%d" %(blockstarts[tc]+1, i+1 )
             cht.SeriesCollection(i).XValues = "=Sheet1!R%dC1:R%dC1" % ( blockstarts[tc]+2,blockends[tc])
             cht.SeriesCollection(i).Values = "=Sheet1!R%dC%d:R%dC%d" % ( blockstarts[tc]+2,i+1,blockends[tc],i+1)
          k+=1

    #Clean original data from template.xls
    xl.ActiveSheet.Range(xl.ActiveSheet.Cells(2, 15), xl.ActiveSheet.Cells(12, 25)).Value = None
    xl.ActiveSheet.Cells(1,1).Select()

    savepath, savename = os.path.split(savepath)
    savename, saveext = os.path.splitext(savename)
    filename = savename + ".xls"
    filename = os.path.join(savepath,filename)
    fn = xl.GetSaveAsFilename(InitialFilename = filename)
    if fn:
       xl.ActiveWorkbook.SaveAs(fn)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print 'No file specified.'
        raw_input("press any key...")
        sys.exit()

    genXLS(sys.argv[1])



