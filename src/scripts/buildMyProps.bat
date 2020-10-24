setlocal EnableDelayedExpansion

SET CWD=%cd%

cd build

SET INCLUDE_COOLPROP=../externals/CoolProp/include

SET INCLUDE_FMT_LIB=../externals/CoolProp/externals/fmtlib

SET TARGET_LIB=MyProps.dll

SET MODELICA_INCLUDE=..\Modelica\Steps\Resources\Include

SET MODELICA_LIBRARY=..\Modelica\Steps\Resources\Library

SET MODELICA_WS[0]="C:\Users\uqxsui\AppData\Local\Temp\OpenModelica\OMEdit\Steps.Cycle.OffDesignRCBCycle_v2"
SET MODELICA_WS[1]="C:\Users\uqxsui\AppData\Local\Temp\OpenModelica\OMEdit\Steps.Test.TestComponentSeries"
SET MODELICA_WS[2]="C:\Users\uqxsui\AppData\Local\Temp\OpenModelica\OMEdit\Steps.Test.TestPCHECImpl"

ECHO OFF

REM build

REM generate dynamic lib

REM libCoolProp.a (static library) should be exsits current directory

ECHO ON

g++ -g -ggdb -c ../MyPropsLib.cpp ../PCHE.cpp -m64 -DBUILDING_DLL -DWITH_SHARED_LIB_WRAPPER -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I%INCLUDE_COOLPROP% -I%INCLUDE_FMT_LIB%

g++ -shared -m64 -L. -o %TARGET_LIB% MyPropsLib.o PCHE.o -lCoolProp

cd %CWD%

copy build\%TARGET_LIB% ..\..\build\

copy build\%TARGET_LIB% ..\..\lib\

copy MyPropsLib.h %MODELICA_INCLUDE%

copy PCHE.h %MODELICA_INCLUDE%

copy build\%TARGET_LIB% %MODELICA_LIBRARY%

REM copy build\%TARGET_LIB% %TARGETDIR%

for /L %%i in (0,1,2) do (
    call echo !MODELICA_WS[%%i]!
    call copy build\%TARGET_LIB% !MODELICA_WS[%%i]!
)