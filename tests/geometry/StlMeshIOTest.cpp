#include "StlMeshIO.h"
#include "GridConversion.h"
#include "GmshMeshReader.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>

#include <chrono>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <set>
#include <string>
#include <vector>

namespace {

class TemporaryDirectory
{
public:
    TemporaryDirectory()
    {
        const auto stamp = std::chrono::steady_clock::now()
            .time_since_epoch().count();
        path_ = std::filesystem::temp_directory_path()
            / ("oneflow_stl_io_" + std::to_string(stamp));
        std::filesystem::create_directories(path_);
    }

    ~TemporaryDirectory()
    {
        if ( preserve_ )
        {
            std::cerr << "Preserving failed test artifacts at: "
                      << path_.string() << std::endl;
            return;
        }
        std::error_code error;
        std::filesystem::remove_all(path_, error);
    }

    void Preserve() const { preserve_ = true; }

    std::filesystem::path File(const std::string& name) const
    {
        return path_ / name;
    }

private:
    std::filesystem::path path_;
    mutable bool preserve_{ false };
};

const std::string kAsciiStl =
    "solid triangle\r\n"
    " facet normal 0 0 1\r\n"
    "  outer loop\r\n"
    "   vertex 0 0 0\r\n"
    "   vertex 1 0 0\r\n"
    "   vertex 0 1 0\r\n"
    "  endloop\r\n"
    " endfacet\r\n"
    "endsolid triangle\r\n";

void AppendU16(std::vector<unsigned char>& bytes, std::uint16_t value)
{
    bytes.push_back(static_cast<unsigned char>(value & 0xff));
    bytes.push_back(static_cast<unsigned char>((value >> 8) & 0xff));
}

void AppendU32(std::vector<unsigned char>& bytes, std::uint32_t value)
{
    for (int shift = 0; shift < 32; shift += 8)
        bytes.push_back(static_cast<unsigned char>((value >> shift) & 0xff));
}

void AppendFloat(std::vector<unsigned char>& bytes, float value)
{
    std::uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    AppendU32(bytes, bits);
}

std::vector<unsigned char> MakeBinaryStl(
    const std::string& header, std::uint16_t attribute = 0)
{
    std::vector<unsigned char> bytes(80, 0);
    for (std::size_t i = 0; i < header.size() && i < 80; ++i)
        bytes[i] = static_cast<unsigned char>(header[i]);
    AppendU32(bytes, 1);
    for (float value : {0.0f, 0.0f, 1.0f,
                        0.0f, 0.0f, 0.0f,
                        1.0f, 0.0f, 0.0f,
                        0.0f, 1.0f, 0.0f})
        AppendFloat(bytes, value);
    AppendU16(bytes, attribute);
    return bytes;
}

void WriteBytes(const std::filesystem::path& path,
                const std::vector<unsigned char>& bytes)
{
    std::ofstream output(path, std::ios::binary);
    output.write(reinterpret_cast<const char*>(bytes.data()),
                 static_cast<std::streamsize>(bytes.size()));
}

void WriteText(const std::filesystem::path& path, const std::string& content)
{
    std::ofstream output(path, std::ios::binary);
    output.write(content.data(), static_cast<std::streamsize>(content.size()));
}

std::vector<unsigned char> ReadBytes(const std::filesystem::path& path)
{
    std::ifstream input(path, std::ios::binary);
    return std::vector<unsigned char>(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>());
}

TEST(StlMeshIOTest, ReadsAsciiFacetAndPreservesVertexOrder)
{
    TemporaryDirectory temp;
    const auto path = temp.File("triangle_ascii.stl");
    WriteText(path, kAsciiStl);
    const auto mesh = ONEFLOW::StlMeshReader::Read(path.string());

    EXPECT_EQ(mesh.header, "triangle");
    ASSERT_EQ(mesh.triangles.size(), 1u);
    EXPECT_DOUBLE_EQ(mesh.triangles[0].normal[2], 1.0);
    EXPECT_DOUBLE_EQ(mesh.triangles[0].vertices[1][0], 1.0);
    EXPECT_DOUBLE_EQ(mesh.triangles[0].vertices[2][1], 1.0);
}

TEST(StlMeshIOTest, ReadsBinaryAndPreservesSolidHeaderAndAttribute)
{
    TemporaryDirectory temp;
    const auto path = temp.File("triangle_binary.stl");
    WriteBytes(path, MakeBinaryStl("solid binary header", 0x1234));
    const auto mesh = ONEFLOW::StlMeshReader::Read(path.string());

    ASSERT_EQ(mesh.triangles.size(), 1u);
    EXPECT_EQ(mesh.header.substr(0, 19), "solid binary header");
    EXPECT_EQ(mesh.triangles[0].attribute, 0x1234);
    EXPECT_DOUBLE_EQ(mesh.triangles[0].vertices[2][1], 1.0);
}

TEST(StlMeshIOTest, RejectsTruncatedBinaryFacet)
{
    TemporaryDirectory temp;
    auto bytes = MakeBinaryStl("truncated");
    bytes.pop_back();
    const auto path = temp.File("truncated.stl");
    WriteBytes(path, bytes);
    EXPECT_THROW(
        ONEFLOW::StlMeshReader::Read(path.string(), ONEFLOW::StlFormat::Binary),
        std::runtime_error);
}

TEST(StlMeshIOTest, RejectsNonFiniteAsciiCoordinates)
{
    TemporaryDirectory temp;
    std::string invalid = kAsciiStl;
    const auto position = invalid.find("vertex 1 0 0");
    ASSERT_NE(position, std::string::npos);
    invalid.replace(position, std::string("vertex 1 0 0").size(),
                    "vertex nan 0 0");
    const auto path = temp.File("non_finite.stl");
    WriteText(path, invalid);
    EXPECT_THROW(ONEFLOW::StlMeshReader::Read(path.string()), std::runtime_error);
}

TEST(StlMeshIOTest, WritesBinaryAndRoundTripsAttributeAndFacetData)
{
    TemporaryDirectory temp;
    ONEFLOW::StlMesh mesh;
    mesh.header = "OneFLOW binary output";
    ONEFLOW::StlTriangle triangle;
    triangle.normal = {0.0, 0.0, 1.0};
    triangle.vertices = {{{0.0, 0.0, 0.0},
                          {1.0, 0.0, 0.0},
                          {0.0, 1.0, 0.0}}};
    triangle.attribute = 0xabcd;
    mesh.triangles.push_back(triangle);

    const auto path = temp.File("round_trip_binary.stl");
    ONEFLOW::StlMeshWriter::WriteBinary(path.string(), mesh);
    const auto bytes = ReadBytes(path);
    ASSERT_EQ(bytes.size(), 134u);
    EXPECT_EQ(bytes[80], 1u);
    EXPECT_EQ(bytes[132], 0xcdu);
    EXPECT_EQ(bytes[133], 0xabu);

    const auto readBack = ONEFLOW::StlMeshReader::Read(path.string());
    ASSERT_EQ(readBack.triangles.size(), 1u);
    EXPECT_EQ(readBack.triangles[0].attribute, 0xabcdu);
    EXPECT_DOUBLE_EQ(readBack.triangles[0].vertices[1][0], 1.0);
}

TEST(StlMeshIOTest, WritesAsciiAndRoundTripsFacetGeometry)
{
    TemporaryDirectory temp;
    ONEFLOW::StlMesh mesh;
    mesh.header = "round trip";
    ONEFLOW::StlTriangle triangle;
    triangle.normal = {0.0, 0.0, 1.0};
    triangle.vertices = {{{0.0, 0.0, 0.0},
                          {1.25, 0.0, 0.0},
                          {0.0, 2.5, 0.0}}};
    mesh.triangles.push_back(triangle);

    const auto path = temp.File("round_trip_ascii.stl");
    ONEFLOW::StlMeshWriter::WriteASCII(path.string(), mesh);
    const auto readBack = ONEFLOW::StlMeshReader::Read(path.string());
    ASSERT_EQ(readBack.triangles.size(), 1u);
    EXPECT_EQ(readBack.header, "round trip");
    EXPECT_DOUBLE_EQ(readBack.triangles[0].vertices[1][0], 1.25);
    EXPECT_DOUBLE_EQ(readBack.triangles[0].vertices[2][1], 2.5);
}

TEST(StlMeshIOTest, RejectsAsciiWriteWhenBinaryAttributesWouldBeLost)
{
    TemporaryDirectory temp;
    ONEFLOW::StlMesh mesh;
    ONEFLOW::StlTriangle triangle;
    triangle.attribute = 1;
    mesh.triangles.push_back(triangle);
    EXPECT_THROW(
        ONEFLOW::StlMeshWriter::WriteASCII(
            temp.File("lossy.stl").string(), mesh),
        std::invalid_argument);
}

TEST(StlMeshIOTest, RejectsBinaryWriterValuesOutsideFloat32Range)
{
    TemporaryDirectory temp;
    ONEFLOW::StlMesh mesh;
    ONEFLOW::StlTriangle triangle;
    triangle.vertices[0][0] = 1.0e100;
    mesh.triangles.push_back(triangle);
    EXPECT_THROW(
        ONEFLOW::StlMeshWriter::WriteBinary(
            temp.File("overflow.stl").string(), mesh),
        std::invalid_argument);
}


ONEFLOW::StlMesh MakeRectangleStl(double width)
{
    ONEFLOW::StlMesh mesh;
    ONEFLOW::StlTriangle first;
    first.normal = {0.0, 0.0, 1.0};
    first.vertices = {{{0.0, 0.0, 0.0},
                       {width, 0.0, 0.0},
                       {width, 1.0, 0.0}}};
    ONEFLOW::StlTriangle second;
    second.normal = {0.0, 0.0, 1.0};
    second.vertices = {{{0.0, 0.0, 0.0},
                        {width, 1.0, 0.0},
                        {0.0, 1.0, 0.0}}};
    mesh.triangles = {first, second};
    return mesh;
}


ONEFLOW::StlMesh MakeCircleStl(int segments, double radius)
{
    ONEFLOW::StlMesh mesh;
    for (int i = 0; i < segments; ++i)
    {
        const double firstAngle = 2.0 * 3.14159265358979323846 * i / segments;
        const double secondAngle = 2.0 * 3.14159265358979323846 * (i + 1) / segments;
        ONEFLOW::StlTriangle triangle;
        triangle.normal = {0.0, 0.0, 1.0};
        triangle.vertices = {{
            {0.0, 0.0, 0.0},
            {radius * std::cos(firstAngle), radius * std::sin(firstAngle), 0.0},
            {radius * std::cos(secondAngle), radius * std::sin(secondAngle), 0.0}}};
        mesh.triangles.push_back(triangle);
    }
    return mesh;
}

void ValidateNativeStlExternalFlowCase(const TemporaryDirectory& temporary)
{
    const auto inputDirectory = temporary.File("input");
    const auto scriptDirectory = temporary.File("script");
    std::filesystem::create_directories(inputDirectory);
    std::filesystem::create_directories(scriptDirectory);
    ONEFLOW::StlMeshWriter::WriteASCII(
        (inputDirectory / "model.stl").string(), MakeCircleStl(32, 0.5));

    const auto layoutPath = scriptDirectory / "layout.json";
    std::ofstream layout(layoutPath);
    ASSERT_TRUE(layout.good());
    layout
        << "{\n"
        << "  \"surface_name\": \"fluid\",\n"
        << "  \"body_boundary_name\": \"wall\",\n"
        << "  \"mesh_size\": 0.15,\n"
        << "  \"outer_mesh_size\": 0.4,\n"
        << "  \"domain\": {\n"
        << "    \"xmin\": -3, \"xmax\": 6, \"ymin\": -3, \"ymax\": 3,\n"
        << "    \"boundaries\": {\n"
        << "      \"left\": \"inflow\", \"right\": \"outflow\",\n"
        << "      \"bottom\": \"farfield\", \"top\": \"farfield\"\n"
        << "    }\n"
        << "  }\n"
        << "}\n";
    layout.close();

    ONEFLOW::GridConfig config;
    config.sourceType = ONEFLOW::GridFileType::STL;
    config.targetType = ONEFLOW::GridFileType::Gmsh;
    config.sourceCaseDir = temporary.File("").string();
    config.sourceFile = "input/model.stl";
    config.layoutFile = "script/layout.json";
    config.targetFile = "results/circle_external_flow.msh";
    config.gmshExecutable = "gmsh";
    config.stlWeldTolerance = 1.0e-10;

    try
    {
        ONEFLOW::ConvertGrid(config, temporary.File("").string());
    }
    catch (const std::exception& error)
    {
        temporary.Preserve();
        std::cerr << "Native STL external-flow Case failed: "
                  << error.what() << std::endl;
        throw;
    }

    const auto meshPath = temporary.File("results/circle_external_flow.msh");
    const auto reportPath =
        temporary.File("results/circle_external_flow.stl-to-gmsh-report.json");
    ASSERT_TRUE(std::filesystem::is_regular_file(meshPath));
    ASSERT_TRUE(std::filesystem::is_regular_file(reportPath));

    std::ifstream meshInput(meshPath, std::ios::binary);
    ASSERT_TRUE(meshInput.good());
    const std::string meshText{
        std::istreambuf_iterator<char>(meshInput),
        std::istreambuf_iterator<char>()};
    const auto mesh = ONEFLOW::ReadGmsh41AsciiMesh(meshText);
    EXPECT_EQ(mesh.version, "4.1");

    std::size_t lineCount = 0;
    std::size_t triangleCount = 0;
    std::size_t quadrilateralCount = 0;
    for (const auto& element : mesh.elements)
    {
        if (element.type == 1) ++lineCount;
        else if (element.type == 2) ++triangleCount;
        else if (element.type == 3) ++quadrilateralCount;
    }
    EXPECT_GE(mesh.nodes.size(), 20u);
    EXPECT_GE(lineCount, 36u);
    EXPECT_GT(triangleCount, 0u);
    EXPECT_EQ(quadrilateralCount, 0u);

    const std::set<std::pair<int, std::string>> expectedGroups{
        {1, "wall"}, {1, "inflow"}, {1, "outflow"},
        {1, "farfield"}, {2, "fluid"}};
    std::set<std::pair<int, std::string>> actualGroups;
    for (const auto& physical : mesh.physicalNames)
        actualGroups.emplace(physical.first.first, physical.second);
    EXPECT_EQ(actualGroups, expectedGroups);
}

void ValidateNativeStlToGmshCase(
    const TemporaryDirectory& temporary,
    double width,
    int cellsX,
    int cellsY,
    std::size_t expectedNodes,
    std::size_t expectedLines,
    std::size_t expectedQuads)
{
    const auto inputDirectory = temporary.File("input");
    const auto scriptDirectory = temporary.File("script");
    std::filesystem::create_directories(inputDirectory);
    std::filesystem::create_directories(scriptDirectory);

    const auto stlPath = inputDirectory / "model.stl";
    ONEFLOW::StlMeshWriter::WriteASCII(stlPath.string(), MakeRectangleStl(width));

    const auto layoutPath = scriptDirectory / "layout.json";
    std::ofstream layout(layoutPath);
    ASSERT_TRUE(layout.good());
    layout
        << "{\n"
        << "  \"surface_name\": \"fluid_domain\",\n"
        << "  \"mesh_size\": 1.0,\n"
        << "  \"transfinite_cells\": {\"x\": " << cellsX
        << ", \"y\": " << cellsY << "},\n"
        << "  \"boundaries\": [\n"
        << "    {\"name\": \"inflow\", \"endpoints\": [[0, 0], [0, 1]]},\n"
        << "    {\"name\": \"outflow\", \"endpoints\": [[" << width
        << ", 0], [" << width << ", 1]]},\n"
        << "    {\"name\": \"symmetry\", \"endpoints\": [[0, 0], ["
        << width << ", 0]]},\n"
        << "    {\"name\": \"farfield\", \"endpoints\": [[0, 1], ["
        << width << ", 1]]}\n"
        << "  ]\n"
        << "}\n";
    layout.close();

    ONEFLOW::GridConfig config;
    config.sourceType = ONEFLOW::GridFileType::STL;
    config.targetType = ONEFLOW::GridFileType::Gmsh;
    config.sourceCaseDir = temporary.File("").string();
    config.sourceFile = "input/model.stl";
    config.layoutFile = "script/layout.json";
    config.targetFile = "results/generated.msh";
    config.gmshExecutable = "gmsh";
    config.stlWeldTolerance = 1.0e-10;

    try
    {
        ONEFLOW::ConvertGrid(config, temporary.File("").string());
    }
    catch ( const std::exception & error )
    {
        temporary.Preserve();
        std::cerr << "STL-to-Gmsh case failed: " << error.what() << std::endl;
        throw;
    }

    const auto meshPath = temporary.File("results/generated.msh");
    const auto reportPath =
        temporary.File("results/generated.stl-to-gmsh-report.json");
    ASSERT_TRUE(std::filesystem::is_regular_file(meshPath));
    ASSERT_TRUE(std::filesystem::is_regular_file(reportPath));

    std::ifstream meshInput(meshPath, std::ios::binary);
    ASSERT_TRUE(meshInput.good());
    const std::string meshText{
        std::istreambuf_iterator<char>(meshInput),
        std::istreambuf_iterator<char>()};
    const auto mesh = ONEFLOW::ReadGmsh41AsciiMesh(meshText);
    EXPECT_EQ(mesh.version, "4.1");
    EXPECT_EQ(mesh.nodes.size(), expectedNodes);

    std::size_t lineCount = 0;
    std::size_t quadCount = 0;
    std::size_t triangleCount = 0;
    for (const auto& element : mesh.elements)
    {
        if (element.type == 1) ++lineCount;
        else if (element.type == 2) ++triangleCount;
        else if (element.type == 3) ++quadCount;
    }
    EXPECT_EQ(lineCount, expectedLines);
    EXPECT_EQ(quadCount, expectedQuads);
    EXPECT_EQ(triangleCount, 0u);

    const std::set<std::pair<int, std::string>> expectedGroups{
        {1, "inflow"}, {1, "outflow"}, {1, "symmetry"},
        {1, "farfield"}, {2, "fluid_domain"}};
    std::set<std::pair<int, std::string>> actualGroups;
    for (const auto& physical : mesh.physicalNames)
        actualGroups.emplace(physical.first.first, physical.second);
    EXPECT_EQ(actualGroups, expectedGroups);
}

#if ONEFLOW_ENABLE_GMSH_INTEGRATION_TESTS
TEST(StlMeshIOTest, NativeStlToGmshGeneratesOneSquareCell)
{
    TemporaryDirectory temporary;
    ValidateNativeStlToGmshCase(temporary, 1.0, 1, 1, 4, 4, 1);
}

TEST(StlMeshIOTest, NativeStlToGmshGeneratesTwoSquareCells)
{
    TemporaryDirectory temporary;
    ValidateNativeStlToGmshCase(temporary, 2.0, 2, 1, 6, 6, 2);
}

TEST(StlMeshIOTest, NativeStlToGmshGeneratesFourCellsFromSquareGeometry)
{
    TemporaryDirectory temporary;
    ValidateNativeStlToGmshCase(temporary, 1.0, 2, 2, 9, 8, 4);
}

TEST(StlMeshIOTest, NativeStlToGmshGeneratesCircleExternalFlowDomain)
{
    TemporaryDirectory temporary;
    ValidateNativeStlExternalFlowCase(temporary);
}
#endif

TEST(StlMeshIOTest, ExtractsSquareBoundaryAndRemovesStlDiagonal)
{
    const auto boundary = ONEFLOW::ExtractStlPlanarBoundary(MakeRectangleStl(1.0));
    ASSERT_EQ(boundary.vertices.size(), 4u);
    EXPECT_EQ(boundary.vertices[0], (std::array<double, 3>{0.0, 0.0, 0.0}));
    EXPECT_EQ(boundary.vertices[1], (std::array<double, 3>{1.0, 0.0, 0.0}));
    EXPECT_EQ(boundary.vertices[2], (std::array<double, 3>{1.0, 1.0, 0.0}));
    EXPECT_EQ(boundary.vertices[3], (std::array<double, 3>{0.0, 1.0, 0.0}));
}

TEST(StlMeshIOTest, ExtractsTwoSquareCellsAsOneOuterRectangle)
{
    const auto boundary = ONEFLOW::ExtractStlPlanarBoundary(MakeRectangleStl(2.0));
    ASSERT_EQ(boundary.vertices.size(), 4u);
    double twiceArea = 0.0;
    for (std::size_t i = 0; i < boundary.vertices.size(); ++i) {
        const auto& a = boundary.vertices[i];
        const auto& b = boundary.vertices[(i + 1) % boundary.vertices.size()];
        twiceArea += a[0] * b[1] - b[0] * a[1];
    }
    EXPECT_DOUBLE_EQ(std::abs(twiceArea) * 0.5, 2.0);
    EXPECT_EQ(
        std::count_if(
            boundary.vertices.begin(), boundary.vertices.end(),
            [](const auto& point) { return point[0] == 1.0; }),
        0);
}

TEST(StlMeshIOTest, RejectsNonPlanarSurfaceForTwoDimensionalBoundaryExtraction)
{
    auto mesh = MakeRectangleStl(1.0);
    mesh.triangles[1].vertices[2][2] = 0.1;
    EXPECT_THROW(
        ONEFLOW::ExtractStlPlanarBoundary(mesh, 1.0e-3),
        std::invalid_argument);
}

TEST(StlMeshIOTest, RejectsNonManifoldTriangleEdge)
{
    ONEFLOW::StlMesh mesh;
    for (const auto& vertices : {
        std::array<std::array<double, 3>, 3>{{
            {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}}},
        std::array<std::array<double, 3>, 3>{{
            {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {1.0, -1.0, 0.0}}},
        std::array<std::array<double, 3>, 3>{{
            {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, -1.0, 0.0}}}}) {
        ONEFLOW::StlTriangle triangle;
        triangle.vertices = vertices;
        mesh.triangles.push_back(triangle);
    }
    EXPECT_THROW(
        ONEFLOW::ExtractStlPlanarBoundary(mesh),
        std::invalid_argument);
}

} // namespace
