set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultigridSolversApplication ****")

######### add include directory here
include_directories( ${CMAKE_SOURCE_DIR}/kratos )
#######################################################################


######### add compilation units here
set( KRATOS_MULTIGRID_SOLVERS_APPLICATION_SOURCES
    # main application compilation unit
    ${CMAKE_CURRENT_SOURCE_DIR}/multigrid_solvers_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_level_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_amg_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_gmg_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_linear_solvers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/multigrid_solvers_python_application.cpp
)
#######################################################################


## define library Kratos which defines the basic python interface
add_library(KratosMultigridSolversApplication SHARED ${KRATOS_MULTIGRID_SOLVERS_APPLICATION_SOURCES})
target_link_libraries(KratosMultigridSolversApplication KratosCore)
set_target_properties(KratosMultigridSolversApplication PROPERTIES PREFIX "")
install(TARGETS KratosMultigridSolversApplication DESTINATION libs )
#######################################################################


## changing the .dll suffix to .pyd for Windows compatibility
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMultigridSolversApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#######################################################################


## install application's python script if required
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)
#######################################################################


## Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MultigridSolversApplication.py" DESTINATION KratosMultiphysics )
#######################################################################

