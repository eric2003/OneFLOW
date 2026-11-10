#include <gtest/gtest.h>
#include "GridConversion.h"
#include "GridTypes.h"
#include <chrono>
#include <filesystem>
#include <string>

#ifdef ENABLE_CGNS
namespace
{
class TemporaryCaseDirectory
{
public:
    TemporaryCaseDirectory()
    {
        const auto stamp = std::chrono::steady_clock::now().time_since_epoch().count();
        path = std::filesystem::temp_directory_path() /
            ("oneflow-gmsh-two-quads-" + std::to_string(stamp));
        std::filesystem::create_directories(path);
    }
    ~TemporaryCaseDirectory()
    {
        std::error_code error;
        std::filesystem::remove_all(path, error);
    }
    const std::filesystem::path & GetPath() const { return path; }
private:
    std::filesystem::path path;
};
}

TEST(GmshTwoQuadsOneFlowTest, ConvertsTwoAdjacentQuadsThroughCgnsIntoNativeGrid)
{
    namespace fs = std::filesystem;
    const fs::path source = fs::path(ONEFLOW_SOURCE_DIR) /
        "grid" / "conversion" / "gmsh_to_cgns" / "examples" / "convert_2d" / "two_quads" / "grid" / "two_quads.msh";
    ASSERT_TRUE(fs::is_regular_file(source)) << source.string();

    TemporaryCaseDirectory temporaryCase;
    const fs::path caseDirectory = temporaryCase.GetPath();
    const fs::path input = caseDirectory / "two_quads.msh";
    const fs::path results = caseDirectory / "results";
    fs::create_directories(results);
    fs::copy_file(source, input);

    ONEFLOW::GridConfig config;
    config.sourceType = ONEFLOW::GridFileType::Gmsh;
    config.targetType = ONEFLOW::GridFileType::OneFLOW;
    config.sourceFile = input.filename().string();
    // CalcGrid resolves relative output names under the configured mesh-output
    // root, not under caseDirectory. Use an absolute path so this isolated
    // conversion test writes only inside its temporary workspace.
    config.targetFile = (results / "two_quads.oneflow.grid").string();
    config.topology = ONEFLOW::GridTopology::Unstructured;

    ASSERT_NO_THROW(ONEFLOW::ConvertGrid(config, caseDirectory.string()));

    const fs::path interchange = results / "two_quads.gmsh-to-oneflow.cgns";
    const fs::path nativeGrid = results / "two_quads.oneflow.grid";
    ASSERT_TRUE(fs::is_regular_file(interchange)) << interchange.string();
    ASSERT_GT(fs::file_size(interchange), 0u);
    ASSERT_TRUE(fs::is_regular_file(nativeGrid)) << nativeGrid.string();
    EXPECT_GT(fs::file_size(nativeGrid), 0u);
}
#else
TEST(GmshTwoQuadsOneFlowTest, CGNSDisabled)
{
    GTEST_SKIP() << "CGNS is not enabled in this build.";
}
#endif
