################################################################################
# Project: RooUnfold
################################################################################

if (${RunningAtlasCI})
  message("Setup for CI in ATLAS analysis release")
  include (CMakeLists_forAtlasCI.txt)
endif()

cmake_minimum_required(VERSION 2.8)

find_package(AnalysisBase QUIET)
find_package(AthAnalysis QUIET)

if (${AnalysisBase_FOUND} OR ${AthAnalysis_FOUND})
  message ("Configuring for build within analysis release")
  set (foundAnalysisRelease  1)
else() 
  set (foundAnalysisRelease  0)
endif()

# The project version number.
project(RooUnfold)
set(VERSION_MAJOR   2   CACHE STRING "Project major version number.")
set(VERSION_MINOR   0   CACHE STRING "Project minor version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR)

find_package( ROOT COMPONENTS Tree Unfold Matrix Hist RIO MathCore Physics RooFitCore RooFit HistFactory Graf Postscript Gpad XMLParser)

file(GLOB RooUnfoldLinkDef src/*_LinkDef.h)
set(RooUnfoldDagostini ${CMAKE_CURRENT_SOURCE_DIR}/src/bayes.for)
file(GLOB RooUnfoldSources src/*.cxx)
file(GLOB RooUnfoldHeaders src/*.h src/*.tpp)
list(REMOVE_ITEM RooUnfoldHeaders ${RooUnfoldLinkDef})
file(GLOB RooUnfoldExecSources examples/*.cxx)

if(EXISTS ${RooUnfoldDagostini})
  enable_language (Fortran)
  set(EXTRA_FLAGS "-DHAVE_DAGOSTINI")
  get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)
  if (Fortran_COMPILER_NAME MATCHES "gfortran")
    set (CMAKE_Fortran_FLAGS "-funroll-all-loops -fno-f2c -O3 -std=f95 -fimplicit-none")
  elseif (Fortran_COMPILER_NAME MATCHES "g77")
    set (CMAKE_Fortran_FLAGS "-funroll-all-loops -fno-f2c -O3 -std=f95 -fimplicit-none -m32")
  endif()
else()
  message("disabling d'Agostini implementation")
  file(GLOB RooUnfoldDagostini src/RooUnfoldDagostini.*)
  list(REMOVE_ITEM RooUnfoldHeaders ${RooUnfoldDagostini})
  list(REMOVE_ITEM RooUnfoldSources ${RooUnfoldDagostini})
  set(RooUnfoldDagostini "")  
endif()

execute_process( COMMAND ln -sf ${RooUnfoldHeaders} -t ${CMAKE_CURRENT_BINARY_DIR} )

file(GLOB Tests "test/*.sh")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_FLAGS} -g -lRooFitCore")

if(${foundAnalysisRelease})
  atlas_subdir( RooUnfold )

  atlas_add_root_dictionary( RooUnfold _dictSource
    ROOT_HEADERS
    ${RooUnfoldHeaders}
    ${RooUnfoldLinkDef})

  # Ensure that the 'RooUnfold' folder expected by atlas_add_library exists at source level
  # by linking the appropriate files from src/ 
  set( _RooUnfold_header_dir ${CMAKE_CURRENT_SOURCE_DIR}/RooUnfold )
  file(MAKE_DIRECTORY ${_RooUnfold_header_dir})
  execute_process( COMMAND ln -sf ${RooUnfoldHeaders} -t ${_RooUnfold_header_dir} )

  
  atlas_add_library( RooUnfold
    ${RooUnfoldHeaders} ${RooUnfoldSources} ${_dictSource}
    PUBLIC_HEADERS RooUnfold
    PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
    PRIVATE_LINK_LIBRARIES ${ROOT_LIBRARIES})

  foreach(ExecSource ${RooUnfoldExecSources})
    get_filename_component(ExecName ${ExecSource} NAME_WE)    
    atlas_add_executable( ${ExecName} ${ExecSource}
      INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} RooUnfold ${CMAKE_CURRENT_SOURCE_DIR}/examples ${CMAKE_CURRENT_SOURCE_DIR}/src
      LINK_LIBRARIES ${ROOT_LIBRARIES} RooUnfold)      
  endforeach()
else()

  execute_process( COMMAND ln -sf ${RooUnfoldHeaders} -t ${CMAKE_CURRENT_BINARY_DIR} )
  set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
  file(WRITE ${SETUP} "#!/bin/bash\n")
  file(APPEND ${SETUP} "# this is an auto-generated setup script\n" )


  # register all the files and directories
  include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/src")
  include_directories ("${ROOT_INCLUDE_DIRS}")

  file(APPEND ${SETUP} "export PATH=\${PATH}:${CMAKE_CURRENT_BINARY_DIR}\n")  
  file(APPEND ${SETUP} "export PYTHONPATH=\${PYTHONPATH}:${CMAKE_CURRENT_BINARY_DIR}\n")
  file(APPEND ${SETUP} "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${CMAKE_CURRENT_BINARY_DIR}\n")  
  
  # generate the dictionary source code
  include(${ROOT_USE_FILE})
  ROOT_GENERATE_DICTIONARY(G__RooUnfold ${RooUnfoldHeaders} LINKDEF ${RooUnfoldLinkDef} OPTIONS ${EXTRA_FLAGS} -Wno-inconsistent-missing-override)

  # register the shared object to include both sources and dictionaries
  add_library( RooUnfold SHARED ${RooUnfoldSources} ${RooUnfoldDagostini} G__RooUnfold.cxx)
 
  # link everything together at the end
  target_link_libraries( RooUnfold ${ROOT_LIBRARIES} )
  
  # Add all targets to the build-tree export set
  export(TARGETS RooUnfold FILE "${PROJECT_BINARY_DIR}/RooUnfoldTargets.cmake")
  
  # Export the package for use from the build-tree
  # (this registers the build-tree with a global CMake-registry)
  export(PACKAGE RooUnfold)
  
  set(CONF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}" "${PROJECT_BINARY_DIR}")
  set(CONF_LIBRARY_DIRS "${PROJECT_BINARY_DIR}")
  set(CONF_LIBRARIES    RooUnfold)
  configure_file(RooUnfoldConfig.cmake.in
    "${PROJECT_BINARY_DIR}/RooUnfoldConfig.cmake" @ONLY)
  
  # Install the RooUnfoldConfig.cmake and RooUnfoldConfigVersion.cmake
  install(FILES
    "${PROJECT_BINARY_DIR}/${CMAKE_FILES_DIRECTORY}/RooUnfoldConfig.cmake"
    "${PROJECT_BINARY_DIR}/RooUnfoldConfigVersion.cmake"
    DESTINATION "${PROJECT_SOURCE_DIR}" COMPONENT dev)

  include(CTest)
  enable_testing()

  foreach(ExecSource ${RooUnfoldExecSources})
    get_filename_component(ExecName ${ExecSource} NAME_WE)    
    add_executable( ${ExecName} ${ExecSource} )
    target_link_libraries ( ${ExecName} RooUnfold ${ROOT_LIBRARIES} )
  endforeach()
endif()

file(GLOB Tests "test/*.sh")

foreach(TestScript ${Tests})
  get_filename_component(TestName ${TestScript} NAME)
  add_test(
    NAME ${TestName}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND bash ${TestScript}
    )
endforeach()

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.cmake)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

