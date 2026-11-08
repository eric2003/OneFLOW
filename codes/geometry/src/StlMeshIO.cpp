#include "StlMeshIO.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace ONEFLOW {
namespace {

using ByteBuffer = std::vector<unsigned char>;

std::string Lower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(),
        [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return value;
}

std::string Trim(const std::string& value)
{
    const std::size_t first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const std::size_t last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

ByteBuffer ReadBytes(const std::string& fileName)
{
    std::ifstream input(fileName, std::ios::binary);
    if (!input) throw std::runtime_error("Cannot open STL input: " + fileName);
    ByteBuffer bytes;
    bytes.assign(
        std::istreambuf_iterator<char>(input),
        std::istreambuf_iterator<char>());
    if (input.bad()) throw std::runtime_error("Failed reading STL input: " + fileName);
    return bytes;
}

std::uint32_t ReadU32(const ByteBuffer& bytes, std::size_t offset)
{
    if (offset > bytes.size() || bytes.size() - offset < 4)
        throw std::runtime_error("Truncated binary STL integer");
    return static_cast<std::uint32_t>(bytes[offset])
        | (static_cast<std::uint32_t>(bytes[offset + 1]) << 8)
        | (static_cast<std::uint32_t>(bytes[offset + 2]) << 16)
        | (static_cast<std::uint32_t>(bytes[offset + 3]) << 24);
}

std::uint16_t ReadU16(const ByteBuffer& bytes, std::size_t offset)
{
    if (offset > bytes.size() || bytes.size() - offset < 2)
        throw std::runtime_error("Truncated binary STL attribute");
    return static_cast<std::uint16_t>(
        static_cast<std::uint16_t>(bytes[offset])
        | (static_cast<std::uint16_t>(bytes[offset + 1]) << 8));
}

float ReadFloat32(const ByteBuffer& bytes, std::size_t offset)
{
    const std::uint32_t bits = ReadU32(bytes, offset);
    float value = 0.0f;
    std::memcpy(&value, &bits, sizeof(value));
    if (!std::isfinite(value))
        throw std::runtime_error("Binary STL contains a non-finite float");
    return value;
}

void AppendU16(ByteBuffer& bytes, std::uint16_t value)
{
    bytes.push_back(static_cast<unsigned char>(value & 0xff));
    bytes.push_back(static_cast<unsigned char>((value >> 8) & 0xff));
}

void AppendU32(ByteBuffer& bytes, std::uint32_t value)
{
    bytes.push_back(static_cast<unsigned char>(value & 0xff));
    bytes.push_back(static_cast<unsigned char>((value >> 8) & 0xff));
    bytes.push_back(static_cast<unsigned char>((value >> 16) & 0xff));
    bytes.push_back(static_cast<unsigned char>((value >> 24) & 0xff));
}

void AppendFloat32(ByteBuffer& bytes, double source)
{
    if (!std::isfinite(source))
        throw std::invalid_argument("Cannot write non-finite STL coordinate or normal");
    const float value = static_cast<float>(source);
    if (!std::isfinite(value))
        throw std::invalid_argument("STL value is outside the float32 range");
    std::uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof(bits));
    AppendU32(bytes, bits);
}

std::string NextToken(std::istringstream& input, const std::string& context)
{
    std::string token;
    if (!(input >> token))
        throw std::runtime_error("Unexpected end of ASCII STL while reading " + context);
    return token;
}

void ExpectToken(
    std::istringstream& input,
    const std::string& expected,
    const std::string& context)
{
    const std::string actual = NextToken(input, context);
    if (Lower(actual) != expected)
        throw std::runtime_error(
            "Expected '" + expected + "' while reading ASCII STL " + context
            + ", got '" + actual + "'");
}

double ReadAsciiReal(std::istringstream& input, const std::string& context)
{
    const std::string token = NextToken(input, context);
    std::size_t consumed = 0;
    double value = 0.0;
    try {
        value = std::stod(token, &consumed);
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid ASCII STL number for " + context + ": " + token);
    }
    if (consumed != token.size() || !std::isfinite(value))
        throw std::runtime_error("Invalid ASCII STL number for " + context + ": " + token);
    return value;
}

StlMesh ReadASCII(const ByteBuffer& bytes)
{
    const std::string content(bytes.begin(), bytes.end());
    std::istringstream input(content);
    ExpectToken(input, "solid", "header");

    StlMesh mesh;
    std::string restOfHeader;
    std::getline(input, restOfHeader);
    mesh.header = Trim(restOfHeader);

    bool ended = false;
    std::string token;
    while (input >> token) {
        const std::string keyword = Lower(token);
        if (keyword == "endsolid") {
            std::string endName;
            std::getline(input, endName);
            input >> std::ws;
            if (input.peek() != std::char_traits<char>::eof())
                throw std::runtime_error("Unexpected data after ASCII STL endsolid");
            ended = true;
            break;
        }
        if (keyword != "facet")
            throw std::runtime_error("Expected ASCII STL facet or endsolid, got '" + token + "'");

        ExpectToken(input, "normal", "facet normal");
        StlTriangle triangle;
        for (double& value : triangle.normal)
            value = ReadAsciiReal(input, "facet normal");
        ExpectToken(input, "outer", "outer loop");
        ExpectToken(input, "loop", "outer loop");
        for (auto& vertex : triangle.vertices) {
            ExpectToken(input, "vertex", "triangle vertex");
            for (double& value : vertex)
                value = ReadAsciiReal(input, "triangle vertex");
        }
        ExpectToken(input, "endloop", "facet endloop");
        ExpectToken(input, "endfacet", "facet endfacet");
        mesh.triangles.push_back(triangle);
    }
    if (!ended) throw std::runtime_error("ASCII STL is missing endsolid");
    return mesh;
}

StlMesh ReadBinary(const ByteBuffer& bytes)
{
    if (bytes.size() < 84)
        throw std::runtime_error("Truncated binary STL header");
    const std::uint32_t count = ReadU32(bytes, 80);
    const std::uint64_t expected =
        84ULL + 50ULL * static_cast<std::uint64_t>(count);
    if (expected != static_cast<std::uint64_t>(bytes.size())) {
        throw std::runtime_error(
            "Binary STL size mismatch: expected " + std::to_string(expected)
            + " bytes, got " + std::to_string(bytes.size()));
    }

    StlMesh mesh;
    mesh.header.assign(
        reinterpret_cast<const char*>(bytes.data()), static_cast<std::size_t>(80));
    mesh.triangles.reserve(static_cast<std::size_t>(count));
    std::size_t offset = 84;
    for (std::uint32_t i = 0; i < count; ++i) {
        StlTriangle triangle;
        for (double& value : triangle.normal) {
            value = ReadFloat32(bytes, offset);
            offset += 4;
        }
        for (auto& vertex : triangle.vertices) {
            for (double& value : vertex) {
                value = ReadFloat32(bytes, offset);
                offset += 4;
            }
        }
        triangle.attribute = ReadU16(bytes, offset);
        offset += 2;
        mesh.triangles.push_back(triangle);
    }
    return mesh;
}

void ValidateTriangle(const StlTriangle& triangle)
{
    for (double value : triangle.normal)
        if (!std::isfinite(value))
            throw std::invalid_argument("STL mesh contains a non-finite normal");
    for (const auto& vertex : triangle.vertices)
        for (double value : vertex)
            if (!std::isfinite(value))
                throw std::invalid_argument("STL mesh contains a non-finite vertex");
}

std::string AsciiSolidName(const std::string& header)
{
    std::string name;
    for (unsigned char c : header) {
        if (c < 32 || c > 126) {
            if (!name.empty()) break;
            continue;
        }
        name.push_back(static_cast<char>(c));
    }
    name = Trim(name);
    return name.empty() ? "OneFLOW" : name;
}

void CheckOutput(std::ofstream& output, const std::string& fileName)
{
    output.flush();
    if (!output) throw std::runtime_error("Failed writing STL output: " + fileName);
}

} // namespace

StlMesh StlMeshReader::Read(const std::string& fileName, StlFormat format)
{
    const ByteBuffer bytes = ReadBytes(fileName);
    if (format == StlFormat::ASCII) return ReadASCII(bytes);
    if (format == StlFormat::Binary) return ReadBinary(bytes);

    // Binary STL has no magic number. Its declared count and exact conventional
    // length are the strongest structural discriminator; the header may say solid.
    if (bytes.size() >= 84) {
        const std::uint32_t count = ReadU32(bytes, 80);
        const std::uint64_t expected =
            84ULL + 50ULL * static_cast<std::uint64_t>(count);
        if (expected == static_cast<std::uint64_t>(bytes.size()))
            return ReadBinary(bytes);
    }
    return ReadASCII(bytes);
}

void StlMeshWriter::WriteASCII(
    const std::string& fileName, const StlMesh& mesh)
{
    for (const StlTriangle& triangle : mesh.triangles) {
        ValidateTriangle(triangle);
        if (triangle.attribute != 0)
            throw std::invalid_argument(
                "ASCII STL cannot preserve nonzero binary facet attributes");
    }

    std::ofstream output(fileName, std::ios::binary);
    if (!output) throw std::runtime_error("Cannot open STL output: " + fileName);
    const std::string name = AsciiSolidName(mesh.header);
    output << "solid " << name << '\n';
    output << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const StlTriangle& triangle : mesh.triangles) {
        output << "  facet normal " << triangle.normal[0] << ' '
               << triangle.normal[1] << ' ' << triangle.normal[2] << '\n';
        output << "    outer loop\n";
        for (const auto& vertex : triangle.vertices) {
            output << "      vertex " << vertex[0] << ' ' << vertex[1]
                   << ' ' << vertex[2] << '\n';
        }
        output << "    endloop\n  endfacet\n";
    }
    output << "endsolid " << name << '\n';
    CheckOutput(output, fileName);
}

void StlMeshWriter::WriteBinary(
    const std::string& fileName, const StlMesh& mesh)
{
    if (mesh.triangles.size() > std::numeric_limits<std::uint32_t>::max())
        throw std::invalid_argument("STL facet count exceeds the binary uint32 limit");
    if (mesh.header.size() > 80)
        throw std::invalid_argument("Binary STL header must not exceed 80 bytes");
    if (mesh.triangles.size() >
        (std::numeric_limits<std::size_t>::max() - 84) / 50)
        throw std::invalid_argument("STL binary output size overflows size_t");

    ByteBuffer bytes;
    bytes.reserve(84 + mesh.triangles.size() * 50);
    for (std::size_t i = 0; i < 80; ++i) {
        bytes.push_back(i < mesh.header.size()
            ? static_cast<unsigned char>(mesh.header[i]) : 0);
    }
    AppendU32(bytes, static_cast<std::uint32_t>(mesh.triangles.size()));
    for (const StlTriangle& triangle : mesh.triangles) {
        ValidateTriangle(triangle);
        for (double value : triangle.normal) AppendFloat32(bytes, value);
        for (const auto& vertex : triangle.vertices)
            for (double value : vertex) AppendFloat32(bytes, value);
        AppendU16(bytes, triangle.attribute);
    }

    std::ofstream output(fileName, std::ios::binary);
    if (!output) throw std::runtime_error("Cannot open STL output: " + fileName);
    output.write(
        reinterpret_cast<const char*>(bytes.data()),
        static_cast<std::streamsize>(bytes.size()));
    CheckOutput(output, fileName);
}

} // namespace ONEFLOW
