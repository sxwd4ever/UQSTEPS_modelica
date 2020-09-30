echo OFF

SET CWD=%cd%

SET APP_NAME=PCHE

SET CP_BUILD_PATH=build

SET INCLUDE_DIR=-I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib

echo ON

cd %CP_BUILD_PATH%

g++ -g -c ../%APP_NAME%.cpp -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER %INCLUDE_DIR% 

g++ -o %APP_NAME%.exe %APP_NAME%.o -L. -lCoolProp
