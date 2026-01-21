#-----------------------------------------------------------------------#
#  file: CMakeLists.txt                                                 #
#                                                                       #
#  version: 1.0   date 9/3/2022                                         #
#                                                                       #
#  Copyright (C) 2022                                                   #
#                                                                       #
#      Enrico Bertolazzi, Francesco Biral and Paolo Bosetti             #
#      Dipartimento di Ingegneria Industriale                           #
#      Università degli Studi di Trento                                 #
#      Via Sommarive 9, I-38123, Trento, Italy                          #
#      email: enrico.bertolazzi@unitn.it                                #
#             francesco.biral@unitn.it                                  #
#             paolo.bosetti@unitn.it                                    #
#-----------------------------------------------------------------------#

cmake_minimum_required( VERSION 3.14 )

project( Splines )

include(./cmake_utils/CMakeLists-common.txt)

# Fix the version warning by providing a default version
if(NOT EB_PROJECT_VERSION)
    set(EB_PROJECT_VERSION "1.1.5")
endif()

project( ${PROJECT_NAME} VERSION ${EB_PROJECT_VERSION} )

if ( NOT DEFINED CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release )
endif()

set(Matlab_ROOT_DIR "/Applications/MATLAB_R2025b.app")
find_package( Matlab REQUIRED COMPONENTS MX_LIBRARY MEX_LIBRARY MAIN_PROGRAM )

# Aggiunge le directory di include di MATLAB
include_directories(${Matlab_INCLUDE_DIRS})

message( STATUS "Matlab_ROOT_DIR = ${Matlab_ROOT_DIR}" )
message( STATUS "PROJECT_NAME    = ${PROJECT_NAME}" )
message( STATUS "Matlab_MEX_LIBRARY = ${Matlab_MEX_LIBRARY}" )
message( STATUS "Matlab_MX_LIBRARY = ${Matlab_MX_LIBRARY}" )

if ( (UNIX OR LINUX) AND NOT APPLE )
  find_package( Threads REQUIRED )
  set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined" )
  set( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -static-libstdc++ -static-libgcc -Wl,--no-undefined" )
endif()

include(./cmake_utils/CMakeLists-cflags.txt)
include(./cmake_utils/CMakeLists-utilities.txt)

# Evaluating the suffix to append to target built by this CMake file
utils_artifacts_suffix(ARTIFACTS_STATIC_SUFFIX TRUE)
if( UTILS_BUILD_SHARED )
  utils_artifacts_suffix(ARTIFACTS_DYNAMIC_SUFFIX FALSE)
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )
message( STATUS "BASE SUFFIXES = ${ARTIFACTS_STATIC_SUFFIX} ${ARTIFACTS_DYNAMIC_SUFFIX}" )

#   ____                                _     _ _
#  |  _ \ _ __ ___ _ __ ___  __ _ _   _(_)___(_) |_ ___
#  | |_) | '__/ _ \ '__/ _ \/ _` | | | | / __| | __/ _ \
#  |  __/| | |  __/ | |  __/ (_| | |_| | \__ \ | ||  __/
#  |_|   |_|  \___|_|  \___|\__, |\__,_|_|___/_|\__\___|
#                              |_|
# MUST be executed before the build of the header list

#   ___
#  / __| ___ _  _ _ _ __ ___ ___
#  \__ \/ _ \ || | '_/ _/ -_|_-<
#  |___/\___/\_,_|_| \__\___/__/
#

set( SOURCES )
file(
  GLOB_RECURSE S
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/src_lib/*.cc
)
foreach(F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF ${CMAKE_CURRENT_SOURCE_DIR} "${F}" )
  list( APPEND SOURCES ${RF} )
endforeach(F ${S})

#   _  _             _
#  | || |___ __ _ __| |___ _ _ ___
#  | __ / -_) _` / _` / -_) '_(_-<
#  |_||_\___\__,_\__,_\___|_| /__/
#
set( PUBLIC_HEADERS )
file(
  GLOB_RECURSE S
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.h*
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.h*
)
foreach (F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF "${CMAKE_CURRENT_SOURCE_DIR}" "${F}" )
  list( APPEND PUBLIC_HEADERS ${RF} )
endforeach(F ${S})

list(FILTER PUBLIC_HEADERS EXCLUDE REGEX "\.DS")

set( PRECOMPILED_HEADERS )
file(
  GLOB_RECURSE S
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*.hh
  ${CMAKE_CURRENT_SOURCE_DIR}/src/*/*.hh
)
foreach (F ${S})
  # build relative path respect to source directory
  file( RELATIVE_PATH RF "${CMAKE_CURRENT_SOURCE_DIR}" "${F}" )
  list( APPEND PRECOMPILED_HEADERS ${RF} )
endforeach(F ${S})

list(FILTER PRECOMPILED_HEADERS EXCLUDE REGEX "\.DS")

# extra include directory for compilation
include_directories( src src_lib src/Utils )

message( STATUS "SOURCES = ${SOURCES}" )

#   _____                  _
#  |_   _|_ _ _ _ __ _ ___| |_ ___
#    | |/ _` | '_/ _` / -_)  _(_-<
#    |_|\__,_|_| \__, \___|\__/__/
#                |___/
#

utils_final_messages()

add_library( Splines STATIC ${SOURCES} )

# Set C++17 for the main library
set_target_properties(Splines PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)

set(
  MEX_LIST
  BaseHermiteMexWrapper
  Spline1DMexWrapper
  Spline1DblendMexWrapper
  Spline2DMexWrapper
  Spline2DblendMexWrapper
  SplineSetMexWrapper
  SplineVecMexWrapper
)

message( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")
if ( CMAKE_BUILD_TYPE MATCHES Debug )
  message( STATUS "DEBUG MODE\n\n")
  add_compile_definitions(DEBUG)
endif()

if ( (UNIX OR LINUX) AND NOT APPLE )
  set( EXTRA_LIBS Threads::Threads ${CMAKE_DL_LIBS} )
else()
  set( EXTRA_LIBS ${CMAKE_DL_LIBS} )
endif()

set( EXTRA_LIBS ${EXTRA_LIBS} Splines )

add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# Set C++ standard to 17 for C++ files only
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

# ----------------------------------------------------------------------
# Calcola automaticamente il suffisso MEX per la piattaforma corrente
# ----------------------------------------------------------------------
if(APPLE)
    # macOS: determina se è Apple Silicon (ARM64) o Intel
    execute_process(
        COMMAND uname -m
        OUTPUT_VARIABLE ARCHITECTURE
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(ARCHITECTURE STREQUAL "arm64" OR ARCHITECTURE STREQUAL "aarch64")
        set(MEX_SUFFIX ".mexmaca64")
        message(STATUS "Detected Apple Silicon (ARM64), using MEX suffix: ${MEX_SUFFIX}")
    else()
        set(MEX_SUFFIX ".mexmaci64")
        message(STATUS "Detected macOS Intel, using MEX suffix: ${MEX_SUFFIX}")
    endif()
elseif(WIN32)
    set(MEX_SUFFIX ".mexw64")
    message(STATUS "Detected Windows, using MEX suffix: ${MEX_SUFFIX}")
elseif(UNIX AND NOT APPLE)
    set(MEX_SUFFIX ".mexa64")
    message(STATUS "Detected Linux, using MEX suffix: ${MEX_SUFFIX}")
else()
    set(MEX_SUFFIX ".mex")
    message(STATUS "Unknown platform, using default MEX suffix: ${MEX_SUFFIX}")
endif()

# ----------------------------------------------------------------------
# Funzione per creare target MEX con nome di output corretto
# ----------------------------------------------------------------------
function(add_mex_target MEX_NAME)
    message(STATUS "Building MEX: ${MEX_NAME}")
    
    # Crea il target MEX
    add_library(mex_${MEX_NAME} MODULE 
        src_mex/mex_${MEX_NAME}.cc
        ${Matlab_ROOT_DIR}/extern/version/c_mexapi_version.c
    )
    
    # Imposta le proprietà con OUTPUT_NAME che include il suffisso
    set_target_properties(mex_${MEX_NAME} PROPERTIES
        PREFIX ""
        OUTPUT_NAME "${MEX_NAME}${MEX_SUFFIX}"  # Nome completo con suffisso
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF
        BUNDLE TRUE
    )
    
    # Aggiungi le directory di include
    target_include_directories(mex_${MEX_NAME} 
        PRIVATE 
        ${Matlab_INCLUDE_DIRS}
        src 
        src_lib 
        src/Utils
    )
    
    # Aggiungi le definizioni di compilazione per MEX
    target_compile_definitions(mex_${MEX_NAME} PRIVATE
        MATLAB_MEX_FILE
        MATLAB_DEFAULT_RELEASE=R2018a
        _USE_MATH_DEFINES
    )
    
    # Link delle librerie
    target_link_libraries(mex_${MEX_NAME}
        ${EXTRA_LIBS}
        ${Matlab_MEX_LIBRARY}
        ${Matlab_MX_LIBRARY}
    )
    
    # Opzioni di linker specifiche per macOS
    if(APPLE)
        target_link_options(mex_${MEX_NAME} PRIVATE
            "-Wl,-undefined,dynamic_lookup"
            "-L${Matlab_ROOT_DIR}/bin/maca64"
            "-Wl,-rpath,${Matlab_ROOT_DIR}/bin/maca64"
        )
    endif()
    
    # Copia il file MEX nella directory bin
    add_custom_command(
        TARGET mex_${MEX_NAME}
        POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "$<TARGET_FILE:mex_${MEX_NAME}>"
        "${CMAKE_CURRENT_SOURCE_DIR}/bin"
        COMMENT "Copying ${MEX_NAME}${MEX_SUFFIX} to bin directory"
    )
endfunction()

# ----------------------------------------------------------------------
# Crea tutti i target MEX
# ----------------------------------------------------------------------
foreach(MEX ${MEX_LIST})
    add_mex_target(${MEX})
endforeach()

# Aggiungere un target per la pulizia
add_custom_target(clean_all
    COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
    COMMENT "Cleaning all build files"
)

# Aggiungere informazioni di debug
message(STATUS "Matlab libraries found:")
message(STATUS "  MEX Library: ${Matlab_MEX_LIBRARY}")
message(STATUS "  MX Library: ${Matlab_MX_LIBRARY}")
message(STATUS "  Include dirs: ${Matlab_INCLUDE_DIRS}")
message(STATUS "  Platform: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "  MEX suffix: ${MEX_SUFFIX}")