set(CMAKE_FOLDER ${CMAKE_FOLDER}/hello)

add_executable(hello_test
  hello_test.cpp
)

target_link_libraries(hello_test
  PRIVATE
    GTest::gtest_main
)

gtest_discover_tests(hello_test
  DISCOVERY_MODE PRE_TEST          # 关键！
)
