set(CMAKE_FOLDER ${CMAKE_FOLDER}/partition)

# The structured partitioner is a pure function over cell coordinates, so the
# test links the one translation unit instead of the whole runtime.
add_executable(StructuredPartitionTest
    StructuredPartitionTest.cpp
    ${CMAKE_SOURCE_DIR}/codes/partition/src/StructuredPartition.cpp
)

target_include_directories(StructuredPartitionTest
    PRIVATE
        ${CMAKE_SOURCE_DIR}/codes/partition/include
        ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(StructuredPartitionTest
    PRIVATE
        GTest::gtest_main
)

gtest_discover_tests(StructuredPartitionTest DISCOVERY_MODE PRE_TEST)
