cmake_minimum_required(VERSION 3.0)

project(gmxtraj)
set(CMAKE_VERBOSE_MAKEFILE OFF)

set(SOURCES 
    src/xdrfile.c
    src/xdrfile_trr.c
    src/xdrfile_xtc.c
    src/gmx_traj.cpp
    )
set(HEADERS
    include/gmx_traj.hpp
    include/xdrfile.h
    include/xdrfile_trr.h
    include/xdrfile_xtc.h
    )

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/libgmxtraj"
           CACHE PATH "default install path" FORCE)
endif()

# Include header files
include_directories(include)

# Optimization flag 
if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "-O2 -Wall  -std=c++11")
    set(CMAKE_C_FLAGS "-O2 -Wno-incompatible-function-pointer-types")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-O2")
    set(CMAKE_C_FLAGS "-O2 -Wno-incompatible-pointer-types")
endif()

# Create shared library
add_library(${PROJECT_NAME} SHARED ${SOURCES})

configure_file(libgmxtraj.pc.in libgmxtraj.pc @ONLY)

# Install library
install(TARGETS ${PROJECT_NAME} DESTINATION lib/)

# Install library headers
install(FILES ${HEADERS} DESTINATION include/)

install(FILES ${CMAKE_BINARY_DIR}/libgmxtraj.pc
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)

install(FILES ${CMAKE_BINARY_DIR}/libgmxtraj.pc
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/pkgconfig)

# example 
install(FILES ${CMAKE_SOURCE_DIR}/example/main.cpp 
    DESTINATION ${CMAKE_INSTALL_PREFIX}/example/)
install(FILES ${CMAKE_SOURCE_DIR}/example/nvt.gro 
    DESTINATION ${CMAKE_INSTALL_PREFIX}/example/)
install(FILES ${CMAKE_SOURCE_DIR}/example/nvt.xtc 
    DESTINATION ${CMAKE_INSTALL_PREFIX}/example/)
configure_file(${CMAKE_SOURCE_DIR}/example/CMakeLists.txt.in 
    ${CMAKE_INSTALL_PREFIX}/example/CMakeLists.txt @ONLY)
