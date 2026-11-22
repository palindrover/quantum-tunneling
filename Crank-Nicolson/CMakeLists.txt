cmake_minimum_required(VERSION 3.18)
project(tunneling LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Тип сборки по умолчанию: Release (только для одноконфигурационных генераторов)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 Debug Release RelWithDebInfo MinSizeRel)
endif()

option(TUNNELING_NATIVE "Оптимизировать под текущий процессор (-march=native)" ON)

# ---------------------------------------------------------------------------
# Зависимости
# ---------------------------------------------------------------------------
find_package(OpenMP REQUIRED)

# LAPACKE: сначала через pkg-config, иначе ищем вручную.
# В обоих случаях получаем одно имя цели TUNNELING_LAPACKE,
# чтобы потом подключать её к любым целям (бенчмарки, тесты, CUDA-версия).
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
    pkg_check_modules(LAPACKE QUIET IMPORTED_TARGET lapacke)
endif()

if(LAPACKE_FOUND)
    message(STATUS "LAPACKE найдена через pkg-config: ${LAPACKE_VERSION}")
    set(TUNNELING_LAPACKE PkgConfig::LAPACKE)
else()
    find_package(LAPACK REQUIRED)
    find_path(LAPACKE_INCLUDE_DIR lapacke.h
              PATH_SUFFIXES lapacke openblas
              REQUIRED)
    find_library(LAPACKE_LIBRARY NAMES lapacke openblas REQUIRED)
    message(STATUS "LAPACKE: ${LAPACKE_LIBRARY}")

    add_library(tunneling_lapacke INTERFACE)
    target_include_directories(tunneling_lapacke INTERFACE ${LAPACKE_INCLUDE_DIR})
    target_link_libraries(tunneling_lapacke INTERFACE ${LAPACKE_LIBRARY} LAPACK::LAPACK)
    set(TUNNELING_LAPACKE tunneling_lapacke)
endif()

# ---------------------------------------------------------------------------
# Флаги компилятора
# ---------------------------------------------------------------------------
include(CheckCXXCompilerFlag)
# Есть в GCC; в Clang только в новых версиях — проверяем, а не надеемся.
check_cxx_compiler_flag(-fcx-limited-range TUNNELING_HAVE_CX_LIMITED_RANGE)

# ---------------------------------------------------------------------------
# Цель
# ---------------------------------------------------------------------------
add_executable(tunneling cn.cpp)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(tunneling PRIVATE -Wall -Wextra -Wpedantic)

    if(TUNNELING_HAVE_CX_LIMITED_RANGE)
        # Комплексное умножение без вызова __muldc3 -> цикл векторизуется
        target_compile_options(tunneling PRIVATE -fcx-limited-range)
    else()
        message(WARNING "-fcx-limited-range не поддерживается: "
                        "циклы с std::complex могут не векторизоваться")
    endif()

    if(TUNNELING_NATIVE)
        # Во всех оптимизированных сборках, кроме Debug
        target_compile_options(tunneling PRIVATE
            $<$<NOT:$<CONFIG:Debug>>:-march=native>)
    endif()
endif()

# OpenMP::OpenMP_CXX сам добавляет -fopenmp и при компиляции, и при линковке
target_link_libraries(tunneling PRIVATE OpenMP::OpenMP_CXX ${TUNNELING_LAPACKE})
