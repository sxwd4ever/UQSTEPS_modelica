#!/bin/bash

cmake .. -DFORCE_BITNESS_64=ON -DCOOLPROP_STATIC_LIBRARY=OFF -DCOOLPROP_OBJECT_LIBRARY=OFF -DCOOLPROP_SHARED_LIBRARY=ON -DCMAKE_BUILD_TYPE=Release

cmake --build .

# test with main.cpp 
g++ -std=c++11 -Wall -O2 -ldl -o main -DCOOLPROP_LIB -I../include -L. main.cpp libCoolProp.so 


g++ -c ../MyPropsLib.cpp -m64  -DBUILDING_DLL -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib -lCoolProp


    g++ -m64 -DWITH_SHARED_LIB_DIRECT -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib -L. -o main.exe main.o libCoolProp.so 