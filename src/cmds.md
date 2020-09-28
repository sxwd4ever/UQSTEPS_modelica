### In Progress

### Static Library - Works alone, not with Modelica

To complaint with CoolName naming convention. MyProps.* work with static library, and MyPropsLib.* (with a -lib suffix) work with shared library. 

- generate static library of CoolProp (libCoolProp.a) by setting 'BUILD_STATIC_LIBRARY=true' in buildCoolProp.bat 

- WITHOUT wrapper - call CoolProps directly

    - update the Main.cpp respectively (include 'MyProps.h in it).

    - call CoolProp directly in Main.cpp

``` bash

g++ -c ../main.cpp -m64 -DWITH_STATIC_LIB -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib

# works

g++ -o main.exe ../main.cpp -m64 -DWITH_STATIC_LIB -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../ -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib -L../../lib -lCoolProp

# also works - flag '-std=c++11' removed

g++ -o main.exe ../main.cpp -m64 -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../ -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib -L../../lib -lCoolProp

```

- WITH Myprops.a as a wrapper

    - generate MyProps.a

``` bash
# compile MyProps

g++ -o MyProps.o -c ../MyProps.cpp -m64 -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../ -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib -L../../lib -lCoolProp

# package it as lib

ar rvs MyProps.a myProps.o

# test it with main.cpp, which is updated accordingly

g++ ../main.cpp -m64 -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../ -I../externals/CoolProp/include -I ../externals/CoolProp/externals/fmtlib -L. MyProps.a -lCoolProp

```

#### Shared Library

- generate shared library of CoolProp (libCoolProp.dll) by setting 'BUILD_STATIC_LIBRARY=false' in buildCoolProp.bat

- For a Demo dll
    
    follow the instructions in [dll with MinGW](http://www.mingw.org/wiki/sampleDLL)

    - prepare dll lib files. 

    example_dll.h

    ```C++
    #ifndef EXAMPLE_DLL_H
    #define EXAMPLE_DLL_H

    #ifdef __cplusplus
    extern "C" {
    #endif

    #ifdef BUILDING_EXAMPLE_DLL
    #define EXAMPLE_DLL __declspec(dllexport)
    #else
    #define EXAMPLE_DLL __declspec(dllimport)
    #endif

    void __stdcall EXAMPLE_DLL hello(const char *s);

    int EXAMPLE_DLL Double(int x);

    #ifdef __cplusplus
    }
    #endif

    // NOTE: this function is not declared extern "C"
    void EXAMPLE_DLL CppFunc(void);

    // NOTE: this class must not be declared extern "C"
    class EXAMPLE_DLL MyClass
    {
    public:
            MyClass() {};
            virtual ~MyClass() {};
            void func(void);
    };

    #endif  // EXAMPLE_DLL_H
    ```

    example_dll.cpp

    ``` C++
    #include <stdio.h>
    #include "example_dll.h"

    __stdcall void hello(const char *s)
    {
            printf("Hello %s\n", s);
    }
    int Double(int x)
    {
            return 2 * x;
    }
    void CppFunc(void)
    {
            puts("CppFunc");
    }
    void MyClass::func(void)
    {
            puts("MyClass.func()");
    }
    ```
    - prepare direct caller main.cpp

    ``` C++
    #include <stdio.h>
    #include "example_dll.h"

    int main(void)
    {
            hello("World");
            printf("%d\n", Double(333));
            CppFunc();

            MyClass a;
            a.func();

            return 0;
    }
    ```
    - compile and run

    ``` bash
    # current working dir is src/build

    # build dll

    g++ -c -DBUILDING_EXAMPLE_DLL ../example_dll.cpp
    g++ -shared -o example_dll.dll example_dll.o -Wl,--out-implib,libexample_dll.a

    # build caller
    g++ -c ../main.cpp
    g++ -o main.exe main.o -L. -lexample_dll

    # see the result
    example_exe
    ```

    - call dll through a wrapper

    ``` bash
    # build the wrapper - dll example_dll.dll generated already

    g++ -c ../MyPropsLib.cpp -m64  -DBUILDING_DLL -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib -lCoolProp

    # this flag '-Wl,--out-implib,libMyPropsLib.a' is not required
    # record it in case need it. 

    g++ -shared -o MyPropsLib.dll MyPropsLib.o -Wl,--out-implib,libMyPropsLib.a -L. -lexample_dll

    # build the (updated) caller 
    
    g++ -c ../main.cpp -m64 -DWITH_SHARED_LIB_WRAPPER -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib 

    g++ -o main.exe main.o -L. -lMyPropsLib  
    ```

    - possible errors

        Case I

        ``` bash
        # declare the MyCall in Header
        ../main.cpp:22:12: error: 'MyCall' was not declared in this scope
        ```

        Case II

        ``` bash
        # recompile MyPropsLib the wrapper
        # the function is declared in header but not compiled in .o
        main.o:main.cpp:(.text+0x10): undefined reference to `__imp_MyCall'
    
        ```   
- For CoolProp

    - **DEAD END** Direct way - need a wrapper

        test if CoolProp the shared lib can work on windows
    
    ``` bash
    g++ -c ../main.cpp -m64 -DWITH_SHARED_LIB_DIRECT -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib

    g++ -m64 -DWITH_SHARED_LIB_DIRECT -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib -L. -o main.exe main.o -lCoolProp
    ```

        I suspect the bug has sth to do with following in MyPropsLib.h

        ``` C
        #ifdef BUILDING_DLL
        #define EXPORT_MY_CODE __declspec(dllexport)
        #else
        #define EXPORT_MY_CODE __declspec(dllimport)
        #endif
        ```

    - through a wrapper dll - USE this one - My wrapper and CoolProp's dll
    
        The combination of my wraper and CoolProp could be:

        1. dll and static lib for my wrapper.(MyProps.a and MyPropsLib.dll)
        2. dll and static lib for CoolProp. (libCoolProp.dll and libCoolProp.a)
        3. call functions in CoolProp.h or CoolPropLib.h

        these will create 8 combinations. And after several trials, I found:

        1. Modelica (windows) works with dll, not static lib. (undefined reference error);        
        2. Compliation for wrapper of MyProps success, while calling the functions declared in CoolProp.h lead to (memory?) crash;
        3. Compilation for wrapper of MyProps success and call functions in CoolPropLib.h works. 

        So, both a dynamic library of CoolProp (libCoolProp.dll) and a dynamic library of my wrapper, which is compiled with static library libCoolProp.a exsits in runtime directory. My Wrapper called the functions declared in CoolPropLib.h

        since Modelica works with coolProp's dll, and myPropsLib works with coolProps.

    ``` bash
    # compile wrapper with libCoolProp.a the static library
    g++ -c ../MyPropsLib.cpp -m64 -DBUILDING_DLL -DWITH_SHARED_LIB_WRAPPER -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib

    g++ -shared  -m64 -L. -o MyPropsLib.dll MyPropsLib.o -lCoolProp
   
    # compile main.cpp the application with MyProps.dll
    g++ -c ../main.cpp -m64 -DWITH_SHARED_LIB_WRAPPER -DSHARED_PTR_TR1_NAMESPACE -DSHARED_PTR_TR1_MEMORY_HEADER -I../externals/CoolProp/include -I../externals/CoolProp/externals/fmtlib


    g++ -m64 -L. -o main.exe main.o -lMyPropsLib -lCoolProp
    ```    

