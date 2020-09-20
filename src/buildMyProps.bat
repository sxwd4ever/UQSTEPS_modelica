SET CWD=%cd%

cd build

SET INCLUDE_COOLPROP=../externals/CoolProp/include

SET INCLUDE_FMT_LIB=../externals/CoolProp/externals/fmtlib

REM build

g++ -c -std=c++11 -L. -L../../lib -I%INCLUDE_COOLPROP% -I%INCLUDE_FMT_LIB% ../myProps.cpp -lCoolProp

REM generate static lib

ar rvs myProps.a myProps.o

cd %CWD%

copy build\myProps.a ..\build\

move build\myProps.a ..\lib\

del build\myProps.o