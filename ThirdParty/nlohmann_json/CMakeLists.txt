cmake_minimum_required(VERSION 3.30)
project(nlohmann_json LANGUAGES CXX)

add_library(nlohmann_json INTERFACE)
target_include_directories(nlohmann_json INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_library(nlohmann_json::nlohmann_json ALIAS nlohmann_json)