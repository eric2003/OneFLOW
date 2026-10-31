set(CMAKE_FOLDER ${CMAKE_FOLDER}/io)

add_executable(TextFileParser_test
  TextFileParser_test.cpp
  control_manifest_reader_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlManifestReader.cpp
  ${CMAKE_SOURCE_DIR}/codes/basic/src/Word.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/CommentLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/TextFileParser.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/SimulationCommandLine.cpp
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
add_executable(ConfigLoader_test
  ConfigLoader_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ConfigLoader.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/LegacyParameterSyntax.cpp
  ${CMAKE_SOURCE_DIR}/codes/basic/src/Word.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/CommentLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/TextFileParser.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/SimulationCommandLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/PrjCmdLineArgs.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/FileUtils.cpp
)

target_include_directories(ConfigLoader_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/io/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
)

target_link_libraries(ConfigLoader_test PRIVATE GTest::gtest_main)
target_compile_features(ConfigLoader_test PRIVATE cxx_std_20)

gtest_discover_tests(ConfigLoader_test DISCOVERY_MODE PRE_TEST)

# ---------------------------------------------------------------------------
# Legacy control-parameter array-size resolver (real DataBase lookup)
# ---------------------------------------------------------------------------
set(CONTROL_RESOLVER_DB_DIR ${CMAKE_SOURCE_DIR}/codes/database)
set(CONTROL_RESOLVER_BASIC_DIR ${CMAKE_SOURCE_DIR}/codes/basic)
set(CONTROL_RESOLVER_PARALLEL_DIR ${CMAKE_SOURCE_DIR}/codes/parallel)
set(CONTROL_RESOLVER_ZONE_DIR ${CMAKE_SOURCE_DIR}/codes/zone)

add_executable(control_configuration_test
  control_parameter_resolver_test.cpp
  control_configuration_ingestion_test.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuContextIO.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuContextBootstrap.cpp
  ${CMAKE_SOURCE_DIR}/codes/main/src/SimuContextCaseEnvironment.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlManifestReader.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlScriptReader.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlFileReader.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlDatabaseIO.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/PIOFileAccess.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/LogFile.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlScriptParser.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ConfigLoader.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/LegacyParameterSyntax.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/TextFileParser.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/CommentLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/SimulationCommandLine.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/PrjCmdLineArgs.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/FileUtils.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ControlParameterResolver.cpp
  ${CMAKE_SOURCE_DIR}/codes/io/src/ConfigDatabaseAdapter.cpp
  ${CONTROL_RESOLVER_DB_DIR}/src/DataBase.cpp
  ${CONTROL_RESOLVER_DB_DIR}/src/DataPara.cpp
  ${CONTROL_RESOLVER_DB_DIR}/src/DataBook.cpp
  ${CONTROL_RESOLVER_DB_DIR}/src/DataField.cpp
  ${CONTROL_RESOLVER_DB_DIR}/src/DataBaseType.cpp
  ${CONTROL_RESOLVER_DB_DIR}/src/DataPage.cpp
  ${CONTROL_RESOLVER_DB_DIR}/src/DataBaseIO.cpp
  ${CONTROL_RESOLVER_BASIC_DIR}/src/Fatal.cpp
  ${CONTROL_RESOLVER_BASIC_DIR}/src/Word.cpp
  ${CONTROL_RESOLVER_PARALLEL_DIR}/src/BasicParallel.cpp
  ${CONTROL_RESOLVER_PARALLEL_DIR}/src/Parallel.cpp
  ${CONTROL_RESOLVER_ZONE_DIR}/src/ZoneState.cpp
)

target_include_directories(control_configuration_test
  PRIVATE
    ${CMAKE_SOURCE_DIR}/codes/io/include
    ${CMAKE_SOURCE_DIR}/codes/main/include
    ${CMAKE_SOURCE_DIR}/codes/accel/include
    ${CONTROL_RESOLVER_DB_DIR}/include
    ${CONTROL_RESOLVER_BASIC_DIR}/include
    ${CMAKE_SOURCE_DIR}/codes/array/include
    ${CONTROL_RESOLVER_PARALLEL_DIR}/include
    ${CONTROL_RESOLVER_ZONE_DIR}/include
    ${CMAKE_SOURCE_DIR}/codes/global/include
)

target_link_libraries(control_configuration_test
  PRIVATE
    GTest::gtest_main
)

target_compile_features(control_configuration_test PRIVATE cxx_std_20)

gtest_discover_tests(control_configuration_test
  DISCOVERY_MODE PRE_TEST
)
