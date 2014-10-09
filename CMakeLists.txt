# CMake build script for LASM.
#
# Copyright (C) 2014.
#
# Authors:
#
#   - Li Dong <dongli@lasg.iap.ac.cn>

cmake_minimum_required (VERSION 2.8)

project (lasm CXX)

option (OPENMP "Turn OpenMP compiler flag ON or OFF" OFF)
option (SHARED "Turn building shared libraries ON of OFF" OFF)

if (OPENMP)
    message ("@@ LASM uses OpenMP compiler flag.")
    find_package (OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    else ()
        message (WARNING "CXX compiler does not support OpenMP!")
    endif ()
else ()
    message ("@@ LASM does not use OpenMP compiler flag.")
endif ()
if (SHARED)
    set (shared_or_static SHARED)
else ()
    set (shared_or_static STATIC)
endif ()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated-register")

# Handle the external libraries.
# MLPACK
if (DEFINED ENV{MLPACK_ROOT})
    message ("@@ Use user provided library MLPACK.")
    message ("@@ MLPACK_ROOT = $ENV{MLPACK_ROOT}")
    set (MLPACK_INCLUDE_DIRS "$ENV{MLPACK_ROOT}/include")
    set (MLPACK_LIBRARY_DIRS "$ENV{MLPACK_ROOT}/lib")
    find_package(LibXml2 2.6.0 REQUIRED)
    include_directories (${MLPACK_INCLUDE_DIRS} ${LIBXML2_INCLUDE_DIR})
    link_directories (${MLPACK_LIBRARY_DIRS})
else ()
    message (FATAL_ERROR
        "CMake couldn't find library MLPACK! "
        "If it have been installed and you know where it is, "
        "set MLPACK_ROOT (e.g. in .bashrc) to it."
    )
endif ()
# GEOMTK
add_subdirectory ("${PROJECT_SOURCE_DIR}/external/geomtk")
include_directories (${GEOMTK_INCLUDE_DIRS})
# BAROTROPIC-MODEL
add_subdirectory ("${PROJECT_SOURCE_DIR}/external/barotropic-model")
include_directories (${BAROTROPIC_MODEL_INCLUDE_DIRS})
# Set the source directory structure.
set (source_directories
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/src/Utils"
    "${PROJECT_SOURCE_DIR}/src/Parcel"
    "${PROJECT_SOURCE_DIR}/src/Advection"
    "${PROJECT_SOURCE_DIR}/src/Advection/Tracer"
    "${PROJECT_SOURCE_DIR}/src/Advection/TestCase"
)
# Collect sources and headers.
foreach (dir ${source_directories})
    include_directories ("${dir}")
    # header files
    file (GLOB tmp1 "${dir}/*.h")
    list (APPEND headers ${tmp1})
    # source files
    aux_source_directory ("${dir}" tmp2)
    list (APPEND sources ${tmp2})
endforeach ()
# Define library targets.
add_library (lasm ${shared_or_static} ${headers} ${sources})
target_link_libraries (lasm geomtk barotropic-model mlpack)
add_dependencies (lasm geomtk barotropic-model)
# Define executable targets.
foreach (dir ${source_directories})
    if (EXISTS "${dir}/demo")
        include_directories ("${dir}/demo")
        aux_source_directory ("${dir}/demo" demo_paths)
        foreach (demo_path ${demo_paths})
            get_filename_component (demo ${demo_path} NAME_WE)
            add_executable (${demo} ${demo_path})
            target_link_libraries (${demo} geomtk lasm barotropic-model)
            add_dependencies (${demo} geomtk lasm barotropic-model)
        endforeach ()
    endif ()
endforeach ()
# Define install rules.
if (${CMAKE_GENERATOR} MATCHES "Unix Makefiles")
    install (TARGETS lasm ARCHIVE DESTINATION lib)
    foreach (header ${headers})
        if (${header} MATCHES "lasm.h")
            install (FILES ${header}
                DESTINATION "include"
            )
        else ()
            install (FILES ${header}
                DESTINATION "include/lasm"
            )
        endif ()
    endforeach ()
endif ()
