cmake_minimum_required(VERSION 3.18)
project(tunneling LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)          

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(TUNNELING_NATIVE "Оптимизировать под текущий процессор (-march=native)" ON)

add_executable(tunneling CN.cpp)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(tunneling PRIVATE -Wall -Wextra -Wpedantic)
    if(TUNNELING_NATIVE)
        target_compile_options(tunneling PRIVATE $<$<CONFIG:Release>:-march=native>)
    endif()
endif()

find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(LAPACKE QUIET IMPORTED_TARGET lapacke)
endif()

if(LAPACKE_FOUND)
    message(STATUS "LAPACKE найдена через pkg-config: ${LAPACKE_VERSION}")
    target_link_libraries(tunneling PRIVATE PkgConfig::LAPACKE)
else()
    find_package(LAPACK REQUIRED)
    find_path(LAPACKE_INCLUDE_DIR lapacke.h
              PATH_SUFFIXES lapacke openblas
              REQUIRED)
    find_library(LAPACKE_LIBRARY NAMES lapacke openblas REQUIRED)
    message(STATUS "LAPACKE: ${LAPACKE_LIBRARY}")

    target_include_directories(tunneling PRIVATE ${LAPACKE_INCLUDE_DIR})
    target_link_libraries(tunneling PRIVATE ${LAPACKE_LIBRARY} LAPACK::LAPACK)
endif()
