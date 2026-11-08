#include "StlMeshIO.h"

#include <gtest/gtest.h>

#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
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
        std::error_code error;
        std::filesystem::remove_all(path_, error);
    }

    std::filesystem::path File(const std::string& name) const
    {
        return path_ / name;
    }

private:
    std::filesystem::path path_;
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

} // namespace
