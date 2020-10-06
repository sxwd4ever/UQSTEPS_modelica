# See also http://stackoverflow.com/a/18697099
cmake_minimum_required (VERSION 2.8.11)
project (main)
set(COOLPROP_STATIC_LIBRARY true)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_COMPILER C:/MinGW/mingw64/mingw64/bin/g++.exe)
set(CMAKE_C_COMPILER C:/MinGW/mingw64/mingw64/bin/gcc.exe)
message(${CMAKE_CXX_COMPILER})
message(${CMAKE_C_COMPILER})

SET(BITNESS "64")

# set(CMAKE_CXX_FLAGS 0)

# add_subdirectory ("${CMAKE_SOURCE_DIR}/../lib" CoolProp)

# message("${CMAKE_SOURCE_DIR}")
set(MY_INCLUDE, "${CMAKE_SOURCE_DIR}/externals/CoolProp/include")
message("${MY_INCLUDE}")

include_directories("${CMAKE_SOURCE_DIR}/externals/CoolProp/include")
include_directories("${CMAKE_SOURCE_DIR}")
include_directories("${CMAKE_SOURCE_DIR}/externals/CoolProp/externals/fmtlib")

add_executable (main "${CMAKE_SOURCE_DIR}/main.cpp")
target_link_libraries(main ${CMAKE_SOURCE_DIR}/../lib/libCoolProp.a)
target_link_libraries(main ${CMAKE_SOURCE_DIR}/../lib/myProps.a)