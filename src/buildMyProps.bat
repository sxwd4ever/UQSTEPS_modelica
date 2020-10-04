SET CWD=%cd%

cd build

SET INCLUDE_COOLPROP=../externals/CoolProp/include

SET INCLUDE_FMT_LIB=../externals/CoolProp/externals/fmtlib

SET TARGET_LIB=MyProps.dll

REM build

REM generate dynamic lib

REM libCoolProp.a (static library) should be exsits current directory

g++ -c ../MyPropsLib.cpp ../PCHE.cpp ../Utils.cpp -m64 -DBUILDING_DLL -DWITH_SHARED_LIB_WRAPPER -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I%INCLUDE_COOLPROP% -I%INCLUDE_FMT_LIB%

g++ -shared  -m64 -L. -o %TARGET_LIB% MyPropsLib.o PCHE.o Utils.o -lCoolProp

cd %CWD%

copy build\%TARGET_LIB% ..\build\

move build\%TARGET_LIB% ..\lib\

del build\%TARGET_LIB%