echo OFF

SET CWD=%cd%

SET APP_NAME=TestPCHE

SET CP_BUILD_PATH=build

SET INCLUDE_DIR=-I.. -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib

echo ON

cd %CP_BUILD_PATH%

g++ -g -ggdb -c ../%APP_NAME%.cpp -DBUILDING_DLL -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER %INCLUDE_DIR% 

g++ -o %APP_NAME%.exe %APP_NAME%.o MyProps.dll -L. -lCoolProp
