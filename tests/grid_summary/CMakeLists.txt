add_executable(GridSummaryTest GridSummaryTest.cpp)

target_include_directories(GridSummaryTest PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/geometry/include
)

target_link_libraries(GridSummaryTest PRIVATE GTest::gtest_main)

gtest_discover_tests(GridSummaryTest DISCOVERY_MODE PRE_TEST)
