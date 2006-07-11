@ECHO OFF

MKDIR build


IF EXIST build\AGEDOLib.lib DEL build\AGEDOLib.lib


:BuildLib

SET CXXFLAGS=-I.. -DNDEBUG -O2 -Oi -Ov -Oc -w-8027

CD build

ECHO    *** Compiling GALib...

bcc32 %CXXFLAGS% -P -c ..\ga\garandom.C ..\ga\gaerror.C ..\ga\GAParameter.C ..\ga\GAStatistics.C ..\ga\GABaseGA.C ..\ga\GASStateGA.C ..\ga\GASimpleGA.C ..\ga\GAIncGA.C ..\ga\GADemeGA.C ..\ga\GADCrowdingGA.C ..\ga\GASelector.C ..\ga\GAScaling.C ..\ga\GAPopulation.C ..\ga\GAGenome.C ..\ga\GABinStr.C ..\ga\gabincvt.C ..\ga\GAAllele.C ..\ga\GAStringGenome.C ..\ga\GA1DBinStrGenome.C ..\ga\GA2DBinStrGenome.C ..\ga\GA3DBinStrGenome.C ..\ga\GABin2DecGenome.C ..\ga\GA1DArrayGenome.C ..\ga\GA2DArrayGenome.C ..\ga\GA3DArrayGenome.C ..\ga\GATreeBASE.C ..\ga\GATree.C ..\ga\GATreeGenome.C ..\ga\GAListBASE.C ..\ga\GAList.C ..\ga\GAListGenome.C

ECHO.

ECHO    *** Compiling GADifferentialEvolution...

bcc32 %CXXFLAGS% -c ..\DE\GADifferentialEvolution.cpp

ECHO.

ECHO    *** Compiling LSODA...

bcc32 %CXXFLAGS% -c ..\LSODA\LivermoreSolver.cpp

ECHO.

ECHO    *** Creating AGEDOLib Library...

tlib AGEDOLib.lib + garandom.obj + gaerror.obj + GAParameter.obj + GAStatistics.obj + GABaseGA.obj + GASStateGA.obj + GASimpleGA.obj + GAIncGA.obj + GADemeGA.obj + GADCrowdingGA.obj + GASelector.obj + GAScaling.obj + GAPopulation.obj + GAGenome.obj + GABinStr.obj + gabincvt.obj + GAAllele.obj + GAStringGenome.obj + GA1DBinStrGenome.obj + GA2DBinStrGenome.obj + GA3DBinStrGenome.obj + GABin2DecGenome.obj + GA1DArrayGenome.obj + GA2DArrayGenome.obj + GA3DArrayGenome.obj + GATreeBASE.obj + GATree.obj + GATreeGenome.obj + GAListBASE.obj + GAList.obj + GAListGenome.obj + GADifferentialEvolution.obj + LivermoreSolver.obj

