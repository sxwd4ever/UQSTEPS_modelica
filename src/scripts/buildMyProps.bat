SET CWD=%cd%

cd build

SET INCLUDE_COOLPROP=../externals/CoolProp/include

SET INCLUDE_FMT_LIB=../externals/CoolProp/externals/fmtlib

SET TARGET_LIB=MyProps.dll

SET MODELICA_INCLUDE=..\Modelica\Steps\Resources\Include

SET MODELICA_LIBRARY=..\Modelica\Steps\Resources\Library

SET MODELICA_TEST=C:\Users\uqxsui\AppData\Local\Temp\OpenModelica\OMEdit\Steps.Utilities.CoolProp

REM build

REM generate dynamic lib

REM libCoolProp.a (static library) should be exsits current directory

g++ -g -ggdb -c ../MyPropsLib.cpp ../PCHE.cpp -m64 -DBUILDING_DLL -DWITH_SHARED_LIB_WRAPPER -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I%INCLUDE_COOLPROP% -I%INCLUDE_FMT_LIB%

g++ -shared -m64 -L. -o %TARGET_LIB% MyPropsLib.o PCHE.o -lCoolProp

cd %CWD%

copy build\%TARGET_LIB% ..\..\build\

copy build\%TARGET_LIB% ..\..\lib\

copy MyPropsLib.h %MODELICA_INCLUDE%

copy PCHE.h %MODELICA_INCLUDE%

copy build\%TARGET_LIB% %MODELICA_LIBRARY%

copy build\%TARGET_LIB% "%MODELICA_TEST%"
