set(CMAKE_FOLDER ${CMAKE_FOLDER}/io)

add_executable(TextFileParser_test
  TextFileParser_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/basic/src/Word.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/CommentLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/TextFileParser.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/PrjCmdLineArgs.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/FileUtils.cpp
  
)

target_include_directories(TextFileParser_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/io/include
	${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(TextFileParser_test PRIVATE GTest::gtest_main)
target_compile_features(TextFileParser_test PRIVATE cxx_std_20)

gtest_discover_tests(TextFileParser_test DISCOVERY_MODE PRE_TEST)