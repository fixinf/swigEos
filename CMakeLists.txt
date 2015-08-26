cmake_minimum_required(VERSION 3.1)
project(eosWrap)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Ofast -std=c++11")

set(SOURCE_FILES
    ap.cpp
    ap.h
    aux.cpp
    aux.h
    constants.h
    DriverBase.cpp
    DriverBase.h
    EnthalpyTOV.cpp
    EnthalpyTOV.h
    eos.cpp
    eos.h
    InterpolatedScalings.cpp
    InterpolatedScalings.h
    KVDriver.cpp
    KVDriver.h
    KVOR.cpp
    KVOR.h
    KVORmod.cpp
    KVORmod.h
    KVORmod2.cpp
    KVORmod2.h
#    KVORmodnos.cpp
#    KVORmodnos.h
#    linalg.cpp
#    linalg.h
    numpy.i
    setconst.cpp
    setconst.h
    solve.cpp
    solve.h
    stdafx.h
    TOV.cpp
    TOV.h
    Walecka.cpp
    Walecka.h
)

include_directories("/opt/levmar-2.6")
include_directories("/usr/include/python3.4m")
include_directories("/usr/lib/python3.4/site-packages/numpy/core/include/")
set(LIBRARY_OUTPUT_PATH ${eosWrap_SOURCE_DIR})

add_library(wrapper OBJECT ${eosWrap_SOURCE_DIR}/wrap.cpp)

add_library(eosWrap SHARED ${SOURCE_FILES} $<TARGET_OBJECTS:wrapper> SCDelta.cpp SCDelta.h)
#target_link_libraries(wrapper python3)
#add_dependencies(wrapper wrapper_src)
add_dependencies(eosWrap wrapper)
target_link_libraries(eosWrap gsl gslcblas levmar lapack)

#add_custom_target(
#wrapper_src
#COMMAND swig -DSWIG_PYTHON_LEGACY_BOOL -c++ -python -o ${eosWrap_SOURCE_DIR}/wrap.cpp -outdir ${eosWrap_BINARY_DIR} ${eosWrap_SOURCE_DIR}/wrap.i
#)

add_custom_command(
OUTPUT ${eosWrap_SOURCE_DIR}/wrap.cpp
COMMAND swig -DSWIG_PYTHON_LEGACY_BOOL -c++ -python -o ${eosWrap_SOURCE_DIR}/wrap.cpp -outdir ${eosWrap_BINARY_DIR} ${eosWrap_SOURCE_DIR}/wrap.i
COMMAND echo "Command run"
DEPENDS ${eosWrap_SOURCE_DIR}/wrap.i
)

add_executable(DirectUse DirectUse.cpp KVOR.h)

target_link_libraries(DirectUse gsl gslcblas levmar lapack eosWrap python3)