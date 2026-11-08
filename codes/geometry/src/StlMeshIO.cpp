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

StlPlanarBoundary ExtractStlPlanarBoundary(
    const StlMesh& mesh, double weldTolerance)
{
    if (!std::isfinite(weldTolerance) || weldTolerance <= 0.0)
        throw std::invalid_argument("STL weld tolerance must be positive and finite");
    const double toleranceSquared = weldTolerance * weldTolerance;
    if (!std::isfinite(toleranceSquared))
        throw std::invalid_argument("STL weld tolerance is too large");

    std::vector<std::array<double, 3>> rawVertices;
    rawVertices.reserve(mesh.triangles.size() * 3);
    std::array<double, 3> minimum{
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()};
    double minimumZ = std::numeric_limits<double>::max();
    double maximumZ = -std::numeric_limits<double>::max();

    for (const StlTriangle& triangle : mesh.triangles) {
        for (const auto& vertex : triangle.vertices) {
            for (double value : vertex) {
                if (!std::isfinite(value))
                    throw std::invalid_argument("STL contains a non-finite vertex");
            }
            rawVertices.push_back(vertex);
            for (std::size_t axis = 0; axis < 3; ++axis)
                minimum[axis] = std::min(minimum[axis], vertex[axis]);
            minimumZ = std::min(minimumZ, vertex[2]);
            maximumZ = std::max(maximumZ, vertex[2]);
        }
    }
    if (rawVertices.size() < 3)
        throw std::invalid_argument("STL contains no triangles");
    if (maximumZ - minimumZ > weldTolerance)
        throw std::invalid_argument("2-D STL boundary must lie in a constant-Z XY plane");

    using BucketKey = std::array<long long, 3>;
    std::map<BucketKey, std::vector<std::size_t>> buckets;
    std::vector<std::array<double, 3>> points;
    std::vector<std::array<std::size_t, 3>> trianglePoints;
    trianglePoints.reserve(mesh.triangles.size());

    for (std::size_t start = 0; start < rawVertices.size(); start += 3) {
        std::array<std::size_t, 3> ids{};
        for (std::size_t corner = 0; corner < 3; ++corner) {
            const auto& point = rawVertices[start + corner];
            BucketKey key{};
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const double relative = (point[axis] - minimum[axis]) / weldTolerance;
                if (!std::isfinite(relative) ||
                    static_cast<long double>(relative) <
                        static_cast<long double>(std::numeric_limits<long long>::min()) + 2.0L ||
                    static_cast<long double>(relative) >
                        static_cast<long double>(std::numeric_limits<long long>::max()) - 2.0L)
                    throw std::invalid_argument(
                        "STL coordinate range is too large for the weld tolerance");
                key[axis] = static_cast<long long>(std::floor(relative));
            }

            std::size_t best = points.size();
            double bestDistance = weldTolerance;
            for (int dx = -1; dx <= 1; ++dx) {
                for (int dy = -1; dy <= 1; ++dy) {
                    for (int dz = -1; dz <= 1; ++dz) {
                        const BucketKey neighbor{
                            key[0] + dx, key[1] + dy, key[2] + dz};
                        const auto bucket = buckets.find(neighbor);
                        if (bucket == buckets.end()) continue;
                        for (std::size_t candidate : bucket->second) {
                            const auto& other = points[candidate];
                            const double distance = std::hypot(
                                point[0] - other[0],
                                point[1] - other[1],
                                point[2] - other[2]);
                            if (distance < bestDistance ||
                                (best != points.size() && distance == bestDistance &&
                                 candidate < best)) {
                                best = candidate;
                                bestDistance = distance;
                            }
                        }
                    }
                }
            }
            if (best == points.size()) {
                best = points.size();
                points.push_back(point);
                buckets[key].push_back(best);
            }
            ids[corner] = best;
        }
        if (ids[0] == ids[1] || ids[1] == ids[2] || ids[2] == ids[0])
            throw std::invalid_argument("Degenerate STL triangle after vertex welding");
        const auto& a = points[ids[0]];
        const auto& b = points[ids[1]];
        const auto& c = points[ids[2]];
        const double twiceArea =
            (b[0] - a[0]) * (c[1] - a[1]) -
            (b[1] - a[1]) * (c[0] - a[0]);
        if (!std::isfinite(twiceArea) || std::abs(twiceArea) <= toleranceSquared)
            throw std::invalid_argument("STL contains a zero-area XY triangle");
        trianglePoints.push_back(ids);
    }

    using Edge = std::pair<std::size_t, std::size_t>;
    std::map<Edge, std::size_t> edgeIncidence;
    for (const auto& triangle : trianglePoints) {
        for (const auto edge : {
            Edge{triangle[0], triangle[1]},
            Edge{triangle[1], triangle[2]},
            Edge{triangle[2], triangle[0]}}) {
            const Edge key = std::minmax(edge.first, edge.second);
            const std::size_t count = ++edgeIncidence[key];
            if (count > 2)
                throw std::invalid_argument("STL contains a non-manifold edge");
        }
    }

    std::map<std::size_t, std::vector<std::size_t>> adjacency;
    for (const auto& entry : edgeIncidence) {
        if (entry.second != 1) continue;
        adjacency[entry.first.first].push_back(entry.first.second);
        adjacency[entry.first.second].push_back(entry.first.first);
    }
    if (adjacency.size() < 3)
        throw std::invalid_argument("STL has no closed outer boundary");
    for (auto& entry : adjacency) {
        auto& neighbors = entry.second;
        if (neighbors.size() != 2)
            throw std::invalid_argument("STL boundary is open, branched, or has touching loops");
        std::sort(neighbors.begin(), neighbors.end());
    }

    const std::size_t first = adjacency.begin()->first;
    std::vector<std::size_t> order{first};
    std::size_t previous = points.size();
    std::size_t current = first;
    while (true) {
        const auto& neighbors = adjacency.at(current);
        const std::size_t next = neighbors[0] == previous
            ? neighbors[1] : neighbors[0];
        if (next == first) break;
        if (std::find(order.begin(), order.end(), next) != order.end())
            throw std::invalid_argument("STL boundary traversal repeats a vertex");
        order.push_back(next);
        previous = current;
        current = next;
        if (order.size() > adjacency.size())
            throw std::invalid_argument("STL boundary traversal did not close");
    }
    if (order.size() != adjacency.size())
        throw std::invalid_argument("STL contains multiple boundary loops or components");

    double twiceArea = 0.0;
    for (std::size_t i = 0; i < order.size(); ++i) {
        const auto& a = points[order[i]];
        const auto& b = points[order[(i + 1) % order.size()]];
        twiceArea += a[0] * b[1] - b[0] * a[1];
    }
    if (!std::isfinite(twiceArea) || std::abs(twiceArea) <= toleranceSquared)
        throw std::invalid_argument("STL boundary encloses zero XY area");
    if (twiceArea < 0.0)
        std::reverse(order.begin() + 1, order.end());

    StlPlanarBoundary boundary;
    boundary.vertices.reserve(order.size());
    for (std::size_t id : order)
        boundary.vertices.push_back(points[id]);
    return boundary;
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
