# ------------------------------------------------------------------------------
# CMake build script for LASM.
#
# Copyright (C) 2014.
#
# Authors:
#
#   - Li Dong <dongli@lasg.iap.ac.cn>
# ------------------------------------------------------------------------------

cmake_minimum_required (VERSION 2.8)

# ------------------------------------------------------------------------------
project (lasm CXX)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-c++11-extensions")

# ------------------------------------------------------------------------------
# external libraries
# mlpack
if (DEFINED ENV{MLPACK_ROOT})
    message ("@@ Use user provided library MLPACK.")
    message ("@@ MLPACK_ROOT = $ENV{MLPACK_ROOT}")
    set (MLPACK_INCLUDE_DIR "$ENV{MLPACK_ROOT}/include")
    set (MLPACK_LIBRARY_DIR "$ENV{MLPACK_ROOT}/lib")
    find_package(LibXml2 2.6.0 REQUIRED)
    include_directories (${MLPACK_INCLUDE_DIR} ${LIBXML2_INCLUDE_DIR})
    link_directories (${MLPACK_LIBRARY_DIR})
else ()
    message (FATAL_ERROR
        "CMake couldn't find library MLPACK! "
        "If it have been installed and you know where it is, "
        "set MLPACK_ROOT (e.g. in .bashrc) to it."
    )
endif ()

# geomtk
add_subdirectory ("${PROJECT_SOURCE_DIR}/external/geomtk")
include_directories (${GEOMTK_INCLUDE_DIRS})

# barotropic-model
add_subdirectory ("${PROJECT_SOURCE_DIR}/external/barotropic-model")
include_directories (${BAROTROPIC_MODEL_INCLUDE_DIRS})

# ------------------------------------------------------------------------------
# source directory structure
set (source_directories
    "${PROJECT_SOURCE_DIR}/src"
    "${PROJECT_SOURCE_DIR}/src/Utils"
    "${PROJECT_SOURCE_DIR}/src/Parcel"
    "${PROJECT_SOURCE_DIR}/src/Advection"
    "${PROJECT_SOURCE_DIR}/src/Advection/Tracer"
    "${PROJECT_SOURCE_DIR}/src/Advection/TestCase"
)

# ------------------------------------------------------------------------------
# collect sources and headers
foreach (dir ${source_directories})
    include_directories ("${dir}")
    # header files
    file (GLOB tmp1 "${dir}/*.h")
    list (APPEND headers ${tmp1})
    # source files
    aux_source_directory ("${dir}" tmp2)
    list (APPEND sources ${tmp2})
endforeach ()

# ------------------------------------------------------------------------------
# library targets
add_library (lasm STATIC ${headers} ${sources})
target_link_libraries (lasm geomtk barotropic-model mlpack)
add_dependencies (lasm geomtk barotropic-model)

# ------------------------------------------------------------------------------
# executable targets
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

# ------------------------------------------------------------------------------
# install rules
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
