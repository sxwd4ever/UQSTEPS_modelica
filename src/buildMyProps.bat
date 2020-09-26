SET CWD=%cd%

cd build

SET INCLUDE_COOLPROP=../externals/CoolProp/include

SET INCLUDE_FMT_LIB=../externals/CoolProp/externals/fmtlib

SET TARGET_LIB=MyProps.dll

REM build

REM generate dynamic lib

gcc -shared -fPIC -std=c++11 -L. -L../../lib -I%INCLUDE_COOLPROP% -I%INCLUDE_FMT_LIB% -c ../MyProps.cpp -o %TARGET_LIB% -lCoolProp

cd %CWD%

copy build\%TARGET_LIB% ..\build\

move build\%TARGET_LIB% ..\lib\

del build\%TARGET_LIB%