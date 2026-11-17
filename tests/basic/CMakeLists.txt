set(CMAKE_FOLDER ${CMAKE_FOLDER}/basic)

add_executable(HXVectorTest HXVectorTest.cpp)

target_include_directories(HXVectorTest
    PRIVATE
        ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(HXVectorTest
    PRIVATE
        GTest::gtest_main
)

gtest_discover_tests(HXVectorTest DISCOVERY_MODE PRE_TEST)

add_executable(MarrayTest MarrayTest.cpp)

target_include_directories(MarrayTest
    PRIVATE
        ${CMAKE_SOURCE_DIR}/codes/array/include
        ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(MarrayTest
    PRIVATE
        GTest::gtest_main
)

gtest_discover_tests(MarrayTest DISCOVERY_MODE PRE_TEST)
