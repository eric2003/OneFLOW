// C++17 Gmsh/OpenCASCADE geometry driver for the OECD/NEA-Vattenfall T-junction.
// Surface serialization and checksum manifests use the shared Python MSH2 exporter.
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
std::string quote(const std::string& value) {
#ifdef _WIN32
  std::string out = "\"";
  for (char ch : value) {
    if (ch == '"') out += "\\\"";
    else out += ch;
  }
  out += "\"";
  return out;
#else
  std::string out = "'";
  for (char ch : value) {
    if (ch == '\'') out += "'\\''";
    else out += ch;
  }
  out += "'";
  return out;
#endif
}
std::string value(int& i, int argc, char** argv) {
  if (i + 1 >= argc) throw std::invalid_argument("missing value after " + std::string(argv[i]));
  return argv[++i];
}
}
int main(int argc, char** argv) {
  try {
    std::map<std::string, std::string> args;
    for (int i = 1; i < argc; ++i) args[argv[i]] = value(i, argc, argv);
    for (const auto* key : {"--output", "--gmsh", "--python", "--extractor"})
      if (!args.count(key)) throw std::invalid_argument(std::string("required option missing: ") + key);
    const std::filesystem::path output = std::filesystem::absolute(args["--output"]);
    const std::filesystem::path geo = output.replace_extension(".geo");
    const std::filesystem::path msh = output.parent_path() / (output.stem().string() + ".surface.msh");
    std::filesystem::create_directories(output.parent_path());
    const double meshSize = args.count("--mesh-size") ? std::stod(args["--mesh-size"]) : 0.03;
    if (!(meshSize > 0.0)) throw std::invalid_argument("mesh size must be positive");
    constexpr double mainRadius = 0.070, branchRadius = 0.050;
    constexpr double mainInlet = 1.070, mainOutlet = 1.070, branchLength = 0.470;
    constexpr double eps = 1.0e-7;
    std::ofstream script(geo);
    if (!script) throw std::runtime_error("cannot create Gmsh geometry script");
    script.precision(12);
    script << "SetFactory(\"OpenCASCADE\");\n"
           << "Geometry.OCCBooleanPreserveNumbering = 0;\n"
           << "Mesh.MshFileVersion = 2.2;\n"
           << "Mesh.MeshSizeMin = " << meshSize << ";\n"
           << "Mesh.MeshSizeMax = " << meshSize << ";\n"
           << "Mesh.MeshSizeFromCurvature = 16;\n"
           << "Mesh.MeshSizeExtendFromBoundary = 1;\n"
           << "Cylinder(1) = {" << -mainInlet << ", 0, 0, " << mainInlet + mainOutlet << ", 0, 0, " << mainRadius << "};\n"
           << "Cylinder(2) = {0, 0, 0, 0, 0, " << branchLength << ", " << branchRadius << "};\n"
           << "fluid[] = BooleanUnion{ Volume{1}; Delete; }{ Volume{2}; Delete; };\n"
           << "surfaces[] = Boundary{ Volume{fluid[]}; };\n"
           << "inlet[] = Surface In BoundingBox{" << -mainInlet-eps << ", " << -mainRadius-eps << ", " << -mainRadius-eps << ", " << -mainInlet+eps << ", " << mainRadius+eps << ", " << mainRadius+eps << "};\n"
           << "outlet[] = Surface In BoundingBox{" << mainOutlet-eps << ", " << -mainRadius-eps << ", " << -mainRadius-eps << ", " << mainOutlet+eps << ", " << mainRadius+eps << ", " << mainRadius+eps << "};\n"
           << "branch_inlet[] = Surface In BoundingBox{" << -branchRadius-eps << ", " << -branchRadius-eps << ", " << branchLength-eps << ", " << branchRadius+eps << ", " << branchRadius+eps << ", " << branchLength+eps << "};\n"
           << "wall[] = surfaces[];\nwall[] -= inlet[];\nwall[] -= outlet[];\nwall[] -= branch_inlet[];\n"
           << "Physical Surface(\"wall\", 1) = {wall[]};\n"
           << "Physical Surface(\"inflow\", 2) = {inlet[], branch_inlet[]};\n"
           << "Physical Surface(\"outflow\", 3) = {outlet[]};\n";
    script.close();
    const std::string gmshCommand = quote(args["--gmsh"]) + " " + quote(geo.string())
      + " -2 -format msh2 -o " + quote(msh.string()) + " -v 2";
    if (std::system(gmshCommand.c_str()) != 0) throw std::runtime_error("Gmsh surface meshing failed");
    const std::string pythonCommand = quote(args["--python"]) + " " + quote(args["--extractor"])
      + " --extract-msh " + quote(msh.string()) + " --output " + quote(output.string())
      + " --generator-kind cpp";
    if (std::system(pythonCommand.c_str()) != 0) throw std::runtime_error("shared MSH2-to-STL exporter failed");
    std::cout << "PASS: native C++ geometry driver generated " << output << '\n';
  } catch (const std::exception& error) {
    std::cerr << "T-junction C++ geometry failed: " << error.what() << '\n';
    return 2;
  }
  return 0;
}
