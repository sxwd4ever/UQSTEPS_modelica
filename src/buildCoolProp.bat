SET CWD=%cd%

SET CP_BUILD_PATH=externals\CoolProp\build

cd %CP_BUILD_PATH%

cmake .. -G "MinGW Makefiles" -DFORCE_BITNESS_64=ON -DCOOLPROP_STATIC_LIBRARY=ON

cmake --build .

cd %CWD%

copy %CP_BUILD_PATH%\libCoolProp.a ..\build

move %CP_BUILD_PATH%\libCoolProp.a ..\lib