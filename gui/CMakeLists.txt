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