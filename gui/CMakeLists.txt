project ( OneFLOW_GUI )

set ( use_qt false )
#set ( use_qt true )
if ( use_qt )
    add_subdirectory ( qt )
endif ()

set ( use_opengl false )
#set ( use_opengl true )

if ( use_opengl )
    add_subdirectory ( opengl )
endif ()

set ( use_vtk false )
#set ( use_vtk true )

if ( use_vtk )
    add_subdirectory ( vtk )
endif ()