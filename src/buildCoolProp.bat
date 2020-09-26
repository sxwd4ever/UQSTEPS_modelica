SET CWD=%cd%

SET CP_BUILD_PATH=externals\CoolProp\build

cd %CP_BUILD_PATH%

REM true: static library; false: shared library

SET BUILD_STATIC_LIB=false

echo %BUILD_STATIC_LIB%

IF /I "%BUILD_STATIC_LIB%" == "true" (

SET CP_TARGET=libCoolProp.a

REM cmake .. -DFORCE_BITNESS_64=ON -DCOOLPROP_STATIC_LIBRARY=ON -DCOOLPROP_OBJECT_LIBRARY=OFF -DCOOLPROP_SHARED_LIBRARY=OFF -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release

cmake .. -G "MinGW Makefiles" -DFORCE_BITNESS_64=ON -DCOOLPROP_STATIC_LIBRARY=ON 

) ELSE (

SET CP_TARGET=libCoolProp.dll 

cmake .. -DFORCE_BITNESS_64=ON -DCOOLPROP_STATIC_LIBRARY=OFF -DCOOLPROP_OBJECT_LIBRARY=OFF -DCOOLPROP_SHARED_LIBRARY=ON -DCOOLPROP_EXTERNC_LIBRARY=ON -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release

) 

cmake --build .

cd %CWD%

copy %CP_BUILD_PATH%\%CP_TARGET% ..\build

move %CP_BUILD_PATH%\%CP_TARGET% ..\lib