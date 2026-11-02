project ( OneFLOW_QT LANGUAGES CXX )

set ( CMAKE_INCLUDE_CURRENT_DIR ON )

set ( CMAKE_AUTOUIC ON )
set ( CMAKE_AUTOMOC ON )
set ( CMAKE_AUTORCC ON )

set ( CMAKE_PREFIX_PATH $ENV{QTDIR} )
find_package ( Qt6Core REQUIRED )
find_package ( Qt6Widgets REQUIRED )

set ( PROJECT_SOURCES
        main.cpp
        mainwindow.cpp
        mainwindow.h
)

add_executable ( ${PROJECT_NAME} WIN32 ${PROJECT_SOURCES} )

target_compile_features ( ${PROJECT_NAME} 
	PUBLIC 
		cxx_std_11
)

target_link_libraries ( ${PROJECT_NAME} 
	PRIVATE
		Qt6::Widgets
)