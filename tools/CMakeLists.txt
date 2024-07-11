cmake_minimum_required(VERSION 3.29)

project(ModernCFD VERSION 0.1.0)

set ( PRJ_COMPILE_FEATURES )
set ( PRJ_COMPILE_DEFINITIONS )
set ( PRJ_LIBRARIES )
set ( PRJ_INCLUDE_DIRS )

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets)

list ( APPEND PRJ_LIBRARIES Qt6::Widgets )
list ( APPEND PRJ_COMPILE_FEATURES cxx_std_23 )

if ( MSVC )
    set_property( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME} )
endif()

set( PROJECT_SOURCES
    main.cpp
    mainwindow.cpp
    mainwindow.h
    mainwindow.ui
    CfdThread.h CfdThread.cpp
    Terminal.h Terminal.cpp
    terminal.ui
)

set( PRJ_SOURCES )
foreach(source ${PROJECT_SOURCES})
  list( APPEND PRJ_SOURCES "codes/${source}" )
endforeach()

foreach(source ${PRJ_SOURCES})
  #message( STATUS "source=${source}" )
endforeach()

add_executable( ${PROJECT_NAME}
)

target_sources( ${PROJECT_NAME}
  PRIVATE
	${PRJ_SOURCES}
)

target_include_directories ( ${PROJECT_NAME}
  PRIVATE
    ${PRJ_INCLUDE_DIRS}
)

target_link_libraries( ${PROJECT_NAME} 
  PRIVATE 
    ${PRJ_LIBRARIES}
)

target_compile_features ( ${PROJECT_NAME} 
  PRIVATE 
    ${PRJ_COMPILE_FEATURES}
)

target_compile_definitions ( ${PROJECT_NAME}
  PRIVATE
    ${PRJ_COMPILE_DEFINITIONS} 
)

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
