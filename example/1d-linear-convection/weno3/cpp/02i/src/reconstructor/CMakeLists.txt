# src/reconstructor/CMakeLists.txt

set(RECONSTRUCTOR_SRC_DIR ${CMAKE_SOURCE_DIR}/src/reconstructor)
set(RECONSTRUCTOR_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/include/reconstructor)

set(RECONSTRUCTOR_HEADERS
    ${RECONSTRUCTOR_INCLUDE_DIR}/Reconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/ReconstructorFactory.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/advanced/EnoReconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/advanced/MPWenoReconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/advanced/WenoReconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/basic/FirstOrderReconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/basic/MUSCLReconstructor.h
    ${RECONSTRUCTOR_INCLUDE_DIR}/basic/SecondOrderReconstructor.h
)

set(RECONSTRUCTOR_SOURCES
    ${RECONSTRUCTOR_SRC_DIR}/ReconstructorFactory.cpp
    ${RECONSTRUCTOR_SRC_DIR}/advanced/EnoReconstructor.cpp
    ${RECONSTRUCTOR_SRC_DIR}/advanced/MPWenoReconstructor.cpp
    ${RECONSTRUCTOR_SRC_DIR}/advanced/WenoReconstructor.cpp
    ${RECONSTRUCTOR_SRC_DIR}/basic/FirstOrderReconstructor.cpp
    ${RECONSTRUCTOR_SRC_DIR}/basic/MUSCLReconstructor.cpp
    ${RECONSTRUCTOR_SRC_DIR}/basic/SecondOrderReconstructor.cpp
)


# 创建通量库（包含公共部分）
add_library(cfd_reconstructor STATIC)
target_sources(cfd_reconstructor PRIVATE ${RECONSTRUCTOR_HEADERS} ${RECONSTRUCTOR_SOURCES})

target_include_directories(cfd_reconstructor PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

source_group("include/reconstructor" FILES
    ${CMAKE_SOURCE_DIR}/include/reconstructor/Reconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/ReconstructorFactory.h
)
source_group("include/reconstructor/advanced" FILES
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/EnoReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/MPWenoReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/advanced/WenoReconstructor.h
)
source_group("include/reconstructor/basic" FILES
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/FirstOrderReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/MUSCLReconstructor.h
    ${CMAKE_SOURCE_DIR}/include/reconstructor/basic/SecondOrderReconstructor.h
)

source_group("src/reconstructor" FILES
    reconstructor/ReconstructorFactory.cpp
)
source_group("src/reconstructor/advanced" FILES
    reconstructor/advanced/EnoReconstructor.cpp
    reconstructor/advanced/MPWenoReconstructor.cpp
    reconstructor/advanced/WenoReconstructor.cpp
)
source_group("src/reconstructor/basic" FILES
    reconstructor/basic/FirstOrderReconstructor.cpp
    reconstructor/basic/MUSCLReconstructor.cpp
    reconstructor/basic/SecondOrderReconstructor.cpp
)