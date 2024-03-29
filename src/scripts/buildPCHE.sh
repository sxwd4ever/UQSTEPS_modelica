#!/bin/bash

APP_NAME=PCHE

CP_BUILD_PATH=build

INCLUDE_DIR='-I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib'

cd $CP_BUILD_PATH

g++ -g -c ../$APP_NAME.cpp -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER $INCLUDE_DIR

g++ -o $APP_NAME.exe $APP_NAME.o -L. -lCoolProp
