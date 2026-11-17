cmake_minimum_required(VERSION 3.30)
project(nanoflann LANGUAGES CXX)

add_library(nanoflann INTERFACE)
target_include_directories(nanoflann INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_library(nanoflann::nanoflann ALIAS nanoflann)