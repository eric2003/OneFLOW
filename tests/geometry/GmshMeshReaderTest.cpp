#include <gtest/gtest.h>
#include "GmshMeshReader.h"
#include <stdexcept>
#include <string>
#include <vector>
namespace {
const std::string mesh=R"MSH($MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
4 4 10 40
0 1 0 1
10
0 0 0
0 2 0 1
20
1 0 0
0 3 0 1
30
1 1 0
0 4 0 1
40
0 1 0
$EndNodes
$Elements
5 6 1 6
1 1 1 1
1 10 20
1 2 1 1
2 20 30
1 3 1 1
3 30 40
1 4 1 1
4 40 10
2 1 2 2
5 10 20 30
6 10 30 40
$EndElements
)MSH";
}
TEST(GmshMeshReaderTest, ReadsCoordinatesEntityAndOrderedConnectivity) {
 auto m=ONEFLOW::ReadGmsh41AsciiMesh(mesh);ASSERT_EQ(m.nodes.size(),4u);ASSERT_EQ(m.elements.size(),6u);
 EXPECT_EQ(m.nodes[0].tag,10);EXPECT_DOUBLE_EQ(m.nodes[2].coordinates[1],1.0);
 EXPECT_EQ(m.elements[4].tag,5);EXPECT_EQ(m.elements[4].type,2);EXPECT_EQ(m.elements[4].entityDimension,2);
 EXPECT_EQ(m.elements[4].nodeTags,(std::vector<long long>{10,20,30}));
}
TEST(GmshMeshReaderTest, RejectsDuplicateNodeTags) {
 auto bad=mesh;auto p=bad.find("20\n1 0 0");ASSERT_NE(p,std::string::npos);bad.replace(p,8,"10\n1 0 0");
 EXPECT_THROW(ONEFLOW::ReadGmsh41AsciiMesh(bad),std::runtime_error);
}
TEST(GmshMeshReaderTest, RejectsUnknownNodeReferences) {
 auto bad=mesh;auto p=bad.find("5 10 20 30");ASSERT_NE(p,std::string::npos);bad.replace(p,10,"5 10 20 999");
 EXPECT_THROW(ONEFLOW::ReadGmsh41AsciiMesh(bad),std::runtime_error);
}
TEST(GmshMeshReaderTest, RejectsUnsupportedElementTypes) {
 auto bad=mesh;auto p=bad.find("2 1 2 2");ASSERT_NE(p,std::string::npos);bad.replace(p,7,"2 1 99 2");
 EXPECT_THROW(ONEFLOW::ReadGmsh41AsciiMesh(bad),std::runtime_error);
}
