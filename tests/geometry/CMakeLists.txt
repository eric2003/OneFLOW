set(CMAKE_FOLDER ${CMAKE_FOLDER}/geometry)

# ==============================================================================
# 1. 第三方库配置 (METIS & CGNS)
# ==============================================================================
option ( METIS_ENABLE "ON for using METIS partitioning library" ON )
if ( METIS_ENABLE )
    list ( APPEND PRJ_COMPILE_DEF ENABLE_METIS )
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
        if ( NOT METIS_ROOT )
            set ( METIS_ROOT "C:/dev/METIS/METIS-${LIB_COMPILER_TYPE}-STATIC" )
        endif()
        if ( NOT METIS_INCLUDE_DIRS )
            set ( METIS_INCLUDE_DIRS "${METIS_ROOT}/include" CACHE PATH "path to METIS headers" )
        endif()
        if ( NOT METIS_LIBRARIES )
            set ( METIS_LIBRARIES "${METIS_ROOT}/lib/metis.lib" CACHE PATH "path to METIS library" )
        endif()
    else()
        set ( METIS_INCLUDE_DIRS $ENV{METIS_HOME_INC} )
        set ( METIS_LIBRARIES $ENV{METIS_HOME_LIB} )
    endif()
endif()

option ( CGNS_ENABLE "ON for using CGNS library" ON )
if ( CGNS_ENABLE )
    if ( ${CMAKE_SYSTEM_NAME} MATCHES "Windows" )
        if ( NOT CGNS_ROOT )
            set ( CGNS_ROOT "C:/dev/cgns/4.4.0" )
        endif()
        set ( CGNS_LIBNAME "cgnsdll.lib" )
        if ( NOT CGNS_INCLUDE_DIRS )
            set ( CGNS_INCLUDE_DIRS "${CGNS_ROOT}/include" CACHE PATH "path to CGNS headers" )
        endif()
        if ( NOT CGNS_LIBRARIES )
            set ( CGNS_LIBRARIES "${CGNS_ROOT}/lib/${CGNS_LIBNAME}" CACHE PATH "path to CGNS library" )
        endif()
    else()
        set ( CGNS_INCLUDE_DIRS $ENV{CGNS_HOME_INC} CACHE PATH "path to CGNS headers" )
        set ( CGNS_LIBRARIES $ENV{CGNS_HOME_LIB} CACHE PATH "path to CGNS library" )
    endif()
    message ( STATUS "The CGNS_INCLUDE_DIRS is ${CGNS_INCLUDE_DIRS}" )
    message ( STATUS "The CGNS_LIBRARIES is ${CGNS_LIBRARIES}" ) 
endif()

# ==============================================================================
# 2. 源码与头文件目录定义
# ==============================================================================
set ( INCLUDE_DIRS
    ${CMAKE_SOURCE_DIR}/codes/geometry/include
    ${CMAKE_SOURCE_DIR}/codes/global/include
    ${CMAKE_SOURCE_DIR}/codes/basic/include
    ${CMAKE_SOURCE_DIR}/codes/io/include
    ${CMAKE_SOURCE_DIR}/codes/task/include
    ${CMAKE_SOURCE_DIR}/codes/cgns/include
    ${CMAKE_SOURCE_DIR}/codes/database/include
    ${CMAKE_SOURCE_DIR}/codes/array/include
    ${CMAKE_SOURCE_DIR}/codes/special/include
    ${CMAKE_SOURCE_DIR}/codes/adt/include
    ${CMAKE_SOURCE_DIR}/codes/zone/include
    ${CMAKE_SOURCE_DIR}/codes/partition/include
    ${CMAKE_SOURCE_DIR}/codes/multigrid/include
    ${CMAKE_SOURCE_DIR}/codes/parallel/include
    ${CMAKE_SOURCE_DIR}/codes/scalar/include
    ${CMAKE_SOURCE_DIR}/codes/main/include
    ${METIS_INCLUDE_DIRS}
    ${CGNS_INCLUDE_DIRS}
)

set ( SRC_DIRS
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/GridGeneration.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/GridWorkflow.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/GridConversion.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/GridTypes.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/GridLayout.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/GridMediator.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/DomainInp.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/Su2Grid.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/Plot3D.cpp
    ${CMAKE_SOURCE_DIR}/codes/partition/src/Partition.cpp
    ${CMAKE_SOURCE_DIR}/codes/basic/src/Word.cpp
    ${CMAKE_SOURCE_DIR}/codes/io/src/Prj.cpp
    ${CMAKE_SOURCE_DIR}/codes/io/src/PrjCmdLineArgs.cpp
    ${CMAKE_SOURCE_DIR}/codes/io/src/CommentLine.cpp
    ${CMAKE_SOURCE_DIR}/codes/io/src/TextFileParser.cpp
    ${CMAKE_SOURCE_DIR}/codes/io/src/FileO.cpp
    ${CMAKE_SOURCE_DIR}/codes/io/src/FileUtils.cpp
    ${CMAKE_SOURCE_DIR}/codes/io/src/LogFile.cpp
    ${CMAKE_SOURCE_DIR}/codes/io/src/PIO.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/ClassicGridWorkflow.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/GridCreate.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/Cavity.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/Rae2822.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/Cylinder.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/CgnsTest.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/GridMachine.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/GridLayoutParser.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/Transfinite.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/CurveLine.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/CurveMachine.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/PointMachine.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/LineMachine.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/DomainMachine.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/BlockMachine.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/SegmentCtrl.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/LineInfo.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/CircleInfo.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/BlockFaceSolver.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/CurveMesh.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/CurveInfo.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/LineMeshImp.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/SimpleDomain.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/MLine.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/SDomain.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/MDomain.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/BlkMesh.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/Block3D.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/Block2D.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/BlockElem.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/LineMesh.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/CircleLineMesh.cpp
    ${CMAKE_SOURCE_DIR}/codes/special/src/CalcCoor.cpp
    ${CMAKE_SOURCE_DIR}/codes/task/src/GridState.cpp
    ${CMAKE_SOURCE_DIR}/codes/task/src/ActionState.cpp
    ${CMAKE_SOURCE_DIR}/codes/task/src/TaskState.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsFactory.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsZbase.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsZbaseUtils.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsGlobal.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsSection.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsZsection.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsZbcBoco.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsBcBoco.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsZone.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsCoor.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsFile.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsBase.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsBaseUtils.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsVariable.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/HXCgns.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsZoneUtils.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsZbc.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsFamilyBc.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsZbc1to1.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsZbcConn.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsBcLink.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsBc1to1.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsBcConn.cpp
    ${CMAKE_SOURCE_DIR}/codes/cgns/src/CgnsPeriod.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/CalcGrid.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/GridElem.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/StrGrid.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/BcRecord.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/PointLocator.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/ElementHome.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/NodeMesh.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/UnsGrid.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/FaceTopo.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/CellTopo.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/InterFace.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/Boundary.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/PointManager.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/UnitElement.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/FaceSolver.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/FaceSearch.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/IFaceLink.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/ElemFeature.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/Grid.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/FaceMesh.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/CellMesh.cpp
    ${CMAKE_SOURCE_DIR}/codes/geometry/src/SlipFace.cpp
    ${CMAKE_SOURCE_DIR}/codes/global/src/Ctrl.cpp
    ${CMAKE_SOURCE_DIR}/codes/global/src/Dimension.cpp
    ${CMAKE_SOURCE_DIR}/codes/global/src/Tolerance.cpp
    ${CMAKE_SOURCE_DIR}/codes/database/src/DataPara.cpp
    ${CMAKE_SOURCE_DIR}/codes/database/src/DataBase.cpp
    ${CMAKE_SOURCE_DIR}/codes/database/src/DataPage.cpp
    ${CMAKE_SOURCE_DIR}/codes/database/src/DataBook.cpp
    ${CMAKE_SOURCE_DIR}/codes/database/src/DataStorage.cpp
    ${CMAKE_SOURCE_DIR}/codes/database/src/DataBaseType.cpp
    ${CMAKE_SOURCE_DIR}/codes/database/src/DataField.cpp
    ${CMAKE_SOURCE_DIR}/codes/database/src/DataBaseIO.cpp
    ${CMAKE_SOURCE_DIR}/codes/multigrid/src/BgGrid.cpp
    ${CMAKE_SOURCE_DIR}/codes/zone/src/ZoneState.cpp
    ${CMAKE_SOURCE_DIR}/codes/zone/src/Zone.cpp
    ${CMAKE_SOURCE_DIR}/codes/zone/src/GridGroup.cpp
    ${CMAKE_SOURCE_DIR}/codes/main/src/SimuBase.cpp
    ${CMAKE_SOURCE_DIR}/codes/parallel/src/BasicParallel.cpp
    ${CMAKE_SOURCE_DIR}/codes/parallel/src/Parallel.cpp
    ${CMAKE_SOURCE_DIR}/codes/scalar/src/ScalarZone.cpp
    ${CMAKE_SOURCE_DIR}/codes/scalar/src/ScalarGrid.cpp
    ${CMAKE_SOURCE_DIR}/codes/scalar/src/ScalarCgns.cpp
    ${CMAKE_SOURCE_DIR}/codes/scalar/src/ScalarIFace.cpp
)

# ==============================================================================
# 3. 核心优化：构建共享静态库 (避免每个测试重复编译底层源码)
# ==============================================================================
add_library(OneFLOW_GeometryCore STATIC ${SRC_DIRS})

target_include_directories(OneFLOW_GeometryCore PUBLIC ${INCLUDE_DIRS})
target_compile_features(OneFLOW_GeometryCore PUBLIC cxx_std_20)

if ( METIS_ENABLE )
    target_compile_definitions(OneFLOW_GeometryCore PUBLIC ENABLE_METIS)
    target_link_libraries(OneFLOW_GeometryCore PUBLIC ${METIS_LIBRARIES})
endif()

if ( CGNS_ENABLE )
    target_compile_definitions(OneFLOW_GeometryCore PUBLIC ENABLE_CGNS)
    target_link_libraries(OneFLOW_GeometryCore PUBLIC ${CGNS_LIBRARIES})
endif()

if ( WIN32 )
    target_compile_definitions(OneFLOW_GeometryCore PUBLIC USE_DLL)
endif()

function(add_oneflow_gtest TEST_NAME TEST_SRC)
    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME} PRIVATE OneFLOW_GeometryCore GTest::gtest_main)
    gtest_discover_tests(${TEST_NAME} DISCOVERY_MODE PRE_TEST)
endfunction()

add_oneflow_gtest(GridGenerationTest GridGenerationTest.cpp)
add_oneflow_gtest(CgnsFactoryLifecycleTest CgnsFactoryLifecycleTest.cpp)
add_oneflow_gtest(FaceSolverLifecycleTest FaceSolverLifecycleTest.cpp)
add_oneflow_gtest(ElemFeatureLifecycleTest ElemFeatureLifecycleTest.cpp)
add_oneflow_gtest(PointLocatorLifecycleTest PointLocatorLifecycleTest.cpp)
add_oneflow_gtest(PointMachineLifecycleTest PointMachineLifecycleTest.cpp)
add_oneflow_gtest(CellMeshLifecycleTest CellMeshLifecycleTest.cpp)
add_oneflow_gtest(BcManagerLifecycleTest BcManagerLifecycleTest.cpp)
add_oneflow_gtest(GridElemLifecycleTest GridElemLifecycleTest.cpp)

add_oneflow_gtest(GridLayoutParserTest GridLayoutParserTest.cpp)
add_oneflow_gtest(GridLayoutGenerationTest GridLayoutGenerationTest.cpp)
add_oneflow_gtest(DomainMachineTest DomainMachineTest.cpp)
target_compile_definitions(GridLayoutParserTest PRIVATE ONEFLOW_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
target_compile_definitions(GridLayoutGenerationTest PRIVATE ONEFLOW_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
