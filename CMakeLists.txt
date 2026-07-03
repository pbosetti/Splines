############################################################################
#                                                                          #
#  Splines -- modern, idiomatic CMake build.                              #
#                                                                          #
#  Self-contained: no dependency on the external `cmake_utils` grab-bag   #
#  of .cmake fragments (it isn't even a proper CMake subproject -- it has #
#  no CMakeLists.txt of its own -- and its Clang flags hardcode           #
#  `-stdlib=libc++`, which breaks on Linux unless libc++-dev happens to   #
#  be installed). Structure mirrors GenericContainer's CMakeLists.txt,    #
#  which this project itself depends on.                                 #
#                                                                          #
############################################################################

cmake_minimum_required( VERSION 3.25 )

project(
  Splines
  VERSION 1.1.5
  DESCRIPTION "A C++ library for univariate and bivariate spline interpolation"
  HOMEPAGE_URL "https://ebertolazzi.github.io/Clothoids/"
  LANGUAGES CXX C
)

# Default to a Release build when the generator uses a single configuration
# and the user hasn't picked one (multi-config generators like Ninja
# Multi-Config or Visual Studio ignore CMAKE_BUILD_TYPE entirely).
if ( NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES )
  set( CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE )
  set_property( CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel" )
endif()

# Standalone vs. consumed-as-a-dependency: tests only make sense standalone,
# and installing shouldn't be forced on a project that FetchContent's us.
option( BUILD_SHARED_LIBS "Build Splines as a shared library" OFF )
option( SPLINES_INSTALL "Generate the install target" ${PROJECT_IS_TOP_LEVEL} )

include( GNUInstallDirs )
include( FetchContent )

if ( APPLE )
  execute_process( COMMAND uname -m COMMAND tr -d '\n' OUTPUT_VARIABLE ARCHITECTURE )
  set( CMAKE_OSX_ARCHITECTURES ${ARCHITECTURE} )
endif()

#   ___                          _                 _
#  |   \ ___ _ __  ___ _ _  __| |___ _ _  __ ___(_)___ ___
#  | |) / -_) '_ \/ -_) ' \/ _` / -_) ' \/ _/ -_)_(_-</ -_)
#  |___/\___| .__/\___|_||_\__,_\___|_||_\__\___(_)__/\___|
#            |_|
#
FetchContent_Declare(
  nlohmann_json
  URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz
  DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable( nlohmann_json )

FetchContent_Declare(
  GenericContainer
  GIT_REPOSITORY https://github.com/pbosetti/GenericContainer
  GIT_TAG        main
  GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable( GenericContainer )
set( GC GenericContainer::Yaml )

# GenericContainer's own CMakeLists.txt doesn't force a standard library for
# Clang, so on Linux it compiles against libstdc++ by default -- but
# UtilsLite/quarticRootsFlocke (below) force libc++ for any Clang via their
# own vendored cmake_utils, and Splines itself matches them (see the Splines
# target further down). Patch GenericContainer's own targets the same way,
# purely from here, so every Clang-compiled object in the final link agrees
# on one standard library.
if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  foreach( _gc_target GenericContainer GenericContainerYaml )
    if ( TARGET ${_gc_target} )
      target_compile_options( ${_gc_target} PUBLIC -stdlib=libc++ )
      target_link_options( ${_gc_target} PUBLIC -stdlib=libc++ )
    endif()
  endforeach()
endif()

# UtilsLite and quarticRootsFlocke (both below) share the same vendored
# cmake_utils option to gate their own internal test suites -- turn it off
# so `ctest` in this project only runs Splines' own test00..test18, not the
# ~70 tests that come with these two dependencies.
set( UTILS_ENABLE_TESTS OFF )

FetchContent_Declare(
  UtilsLite
  GIT_REPOSITORY https://github.com/pbosetti/UtilsLite.git
  GIT_TAG        main
)
FetchContent_MakeAvailable( UtilsLite )
set( UTILS utils::UtilsLite_Static )
# Deliberately not redirecting Eigen3::Eigen at UtilsLite's vendored copy
# (src/Utils/3rd/Eigen) -- UtilsLite's own CMakeLists.txt doesn't reference
# that target at all, so GenericContainer's Eigen3::Eigen (fetched above,
# GLOBAL IMPORTED) is used everywhere, with no duplicate Eigen in the tree.

FetchContent_Declare(
  quarticRootsFlocke
  GIT_REPOSITORY https://github.com/ebertolazzi/quarticRootsFlocke.git
  GIT_TAG        b13bd707f713e19d4d1e7c22ff8c00c60b34e90b
)
FetchContent_MakeAvailable( quarticRootsFlocke )
set( ROOTS utils_quarticRootsFlocke_Static )

message( STATUS "----------------------------------------------------------------" )
message( STATUS "UTILS = ${UTILS}" )
message( STATUS "ROOTS = ${ROOTS}" )
message( STATUS "GC    = ${GC} Eigen3::Eigen" )
message( STATUS "----------------------------------------------------------------" )

#   _     _ _
#  | |   (_) |__  _ __ __ _ _ __ _   _
#  | |   | | '_ \| '__/ _` | '__| | | |
#  | |___| | |_) | | | (_| | |  | |_| |
#  |_____|_|_.__/|_|  \__,_|_|   \__, |
#                                |___/
#
add_library(
  Splines
  src/SplineAkima.cc
  src/SplineBiCubic.cc
  src/SplineBiCubicBase.cc
  src/SplineBilinear.cc
  src/SplineBiQuintic.cc
  src/SplineBiQuinticBase.cc
  src/SplineBuild.cc
  src/SplineConstant.cc
  src/SplineCubic.cc
  src/SplineCubicBase.cc
  src/SplineHermite.cc
  src/SplineLinear.cc
  src/SplineParametrization.cc
  src/SplinePchip.cc
  src/SplineQuintic.cc
  src/SplineQuinticBase.cc
  src/Splines.cc
  src/Splines1D.cc
  src/SplinesCinterface.cc
  src/SplineSet.cc
  src/SplineSurf.cc
  src/SplineUtils.cc
  src/SplineVanLeer.cc
  src/SplineVec.cc
)
add_library( Splines::Splines ALIAS Splines )

target_compile_features( Splines PUBLIC cxx_std_20 )
target_compile_definitions( Splines PUBLIC SPLINES_NO_COMPATIBILITY )
target_include_directories(
  Splines PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
  "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)
target_link_libraries( Splines PUBLIC ${ROOTS} ${GC} ${UTILS} nlohmann_json::nlohmann_json )

set_target_properties(
  Splines PROPERTIES
  VERSION     ${PROJECT_VERSION}
  SOVERSION   ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME Splines
)

if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang|GNU" )
  target_compile_options( Splines PRIVATE -Wall -Wextra )
endif()

# quarticRootsFlocke and UtilsLite build via their own vendored cmake_utils,
# which forces -stdlib=libc++ for any Clang (harmless/redundant on Apple,
# where Clang already defaults to libc++). Match it here so Splines' own
# object files and the final test executables agree with those static
# libraries' ABI -- without it, linking mixes libc++ (std::__1::...) and
# libstdc++ symbols and fails. On Linux this needs libc++-dev/libc++abi-dev
# installed alongside Clang (see README).
if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  target_compile_options( Splines PUBLIC -stdlib=libc++ )
  target_link_options( Splines PUBLIC -stdlib=libc++ )
endif()

# Splines' own public headers pull in fmt (Utils_fmt.hh, used throughout for
# error messages), and fmt static_asserts on MSVC unless the source is
# compiled as UTF-8 explicitly (MSVC's default source/execution charset
# isn't UTF-8). PUBLIC because anything that includes Splines' headers --
# including our own test executables -- hits the same fmt headers.
# /bigobj avoids "too many sections" from the heavy SplineSet_eval*.hxx
# template instantiations.
if ( MSVC )
  target_compile_options( Splines PUBLIC /utf-8 /bigobj )
endif()

message( STATUS "Compiler used: ${CMAKE_CXX_COMPILER_ID}" )

#   _____         _
#  |_   _|__  ___| |_ ___
#    | |/ _ \/ __| __/ __|
#    | |  __/\__ \ |_\__ \
#    |_|\___||___/\__|___/
#
# Tests only make sense when Splines is the top-level project -- a parent
# project that FetchContent's us shouldn't get our test executables too.
include( CTest )
if ( BUILD_TESTING AND PROJECT_IS_TOP_LEVEL )
  set( EXECUTABLE_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/bin" )
  set(
    SPLINES_TESTS
    test00 test01 test02 test03 test04 test05 test06 test07 test08 test09
    test10 test11 test12 test13 test14 test15 test16 test17 test18
  )
  add_custom_target( "${PROJECT_NAME}_all_tests" ALL )
  foreach( S ${SPLINES_TESTS} )
    add_executable( ${S} "${CMAKE_CURRENT_SOURCE_DIR}/src_tests/${S}.cc" )
    target_link_libraries( ${S} PRIVATE Splines::Splines ${CMAKE_DL_LIBS} )
    add_test( NAME "${S}" COMMAND ./bin/${S} WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
    add_dependencies( "${PROJECT_NAME}_all_tests" ${S} )
  endforeach()
endif()

#   ___         _        _ _
#  |_ _|_ _  __| |_ __ _| | |
#   | || ' \(_-<  _/ _` | | |
#  |___|_||_/__/\__\__,_|_|_|
#
# No `install(EXPORT ...)` / `SplinesConfig.cmake` here: Splines PUBLIC-links
# GenericContainer/UtilsLite/quarticRootsFlocke/nlohmann_json targets that are
# themselves plain (non-IMPORTED) targets built in-tree via FetchContent, and
# none of them export themselves when nested as someone else's subdependency
# (GenericContainer skips its own install() entirely unless it's the
# top-level project, by design -- see its own README). CMake requires every
# target reachable through an exported target's link interface to belong to
# some export set, which none of these do, so `install(EXPORT)` hard-fails
# configure. `FetchContent` (see README) is the supported consumption path;
# this section still gives a real, working `install()` of the compiled
# library and headers, plus a CPack `package` target.
if ( SPLINES_INSTALL )
  install(
    TARGETS Splines
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )
  install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/include/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  # `cpack` / the "package" target: source and binary archives.
  set( CPACK_PACKAGE_NAME "${PROJECT_NAME}" )
  set( CPACK_PACKAGE_VERSION "${PROJECT_VERSION}" )
  set( CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}" )
  set( CPACK_GENERATOR "TGZ;ZIP" )
  set( CPACK_SOURCE_GENERATOR "TGZ;ZIP" )
  set( CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}/;/\\\\.git/;/build.*/" )
  include( CPack )
endif()
