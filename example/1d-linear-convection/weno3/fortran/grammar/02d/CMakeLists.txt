cmake_minimum_required(VERSION 3.20)

project ( testprj )

set ( PRJ_COMPILE_FEATURES )
set ( PRJ_COMPILE_DEFINITIONS )

enable_language(Fortran)

list ( APPEND PRJ_COMPILE_FEATURES cxx_std_20 )

add_subdirectory(src)
add_subdirectory(tests)



