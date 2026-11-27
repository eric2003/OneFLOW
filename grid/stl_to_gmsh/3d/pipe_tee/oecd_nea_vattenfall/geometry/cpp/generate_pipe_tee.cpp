// C++17 Gmsh/OpenCASCADE geometry driver for the OECD/NEA-Vattenfall T-junction.
// Surface serialization and checksum manifests use the shared Python MSH2 exporter.
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _WIN32
#include <process.h>
#endif

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

int run_program(const std::string& executable, const std::vector<std::string>& arguments) {
#ifdef _WIN32
  // Avoid cmd.exe quoting rules for paths and arguments on Windows.
  std::vector<const char*> argv;
  argv.reserve(arguments.size() + 2);
  argv.push_back(executable.c_str());
  for (const auto& argument : arguments) argv.push_back(argument.c_str());
  argv.push_back(nullptr);
  return _spawnv(_P_WAIT, executable.c_str(), argv.data());
#else
  std::string command = quote(executable);
  for (const auto& argument : arguments) command += " " + quote(argument);
  return std::system(command.c_str());
#endif
}
}
int main(int argc, char** argv) {
  try {
    std::map<std::string, std::string> args;
    for (int i = 1; i < argc; ++i) {
      const std::string option = argv[i];
      args[option] = value(i, argc, argv);
    }
    for (const auto* key : {"--output", "--gmsh", "--python", "--extractor"})
      if (!args.count(key)) throw std::invalid_argument(std::string("required option missing: ") + key);
    std::filesystem::path output = std::filesystem::absolute(args["--output"]);
    std::filesystem::path geo = output;
    geo.replace_extension(".geo");
    const std::filesystem::path msh = output.parent_path() / (output.stem().string() + ".surface.msh");
    std::filesystem::create_directories(output.parent_path());
    const double meshSize = args.count("--mesh-size") ? std::stod(args["--mesh-size"]) : 0.03;
    if (!(meshSize > 0.0)) throw std::invalid_argument("mesh size must be positive");
    const double mainRadius = (args.count("--main-diameter") ? std::stod(args["--main-diameter"]) : 0.140) / 2.0;
    const double branchRadius = (args.count("--branch-diameter") ? std::stod(args["--branch-diameter"]) : 0.100) / 2.0;
    const double mainInlet = args.count("--main-inlet-length") ? std::stod(args["--main-inlet-length"]) : 1.220;
    const double mainOutlet = args.count("--main-outlet-length") ? std::stod(args["--main-outlet-length"]) : 1.220;
    const double branchLength = args.count("--branch-inlet-length") ? std::stod(args["--branch-inlet-length"]) : 0.693;
    // Match the Python generator: OCC face bounds need a finite selection slab.
    const double eps = std::max(1.0e-5, meshSize * 1.0e-3);
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
           << "Physical Surface(\"main_inflow\", 2) = {inlet[]};\n"
           << "Physical Surface(\"branch_inflow\", 3) = {branch_inlet[]};\n"
           << "Physical Surface(\"outflow\", 4) = {outlet[]};\n"
           << "Physical Volume(\"fluid_domain\", 1) = {fluid[]};\n";
    script.close();
    if (run_program(args["--gmsh"], {
          geo.string(), "-2", "-format", "msh2", "-o", msh.string(), "-v", "2"
        }) != 0) {
      throw std::runtime_error("Gmsh surface meshing failed");
    }
    if (run_program(args["--python"], {
          args["--extractor"], "--extract-msh", msh.string(), "--output", output.string(),
          "--generator-kind", "cpp", "--generator-source",
          args.count("--generator-source") ? args["--generator-source"] : std::string(__FILE__),
          "--mesh-size", std::to_string(meshSize), "--validate-pipe-tee",
          "--main-diameter", std::to_string(mainRadius * 2.0),
          "--branch-diameter", std::to_string(branchRadius * 2.0),
          "--main-inlet-length", std::to_string(mainInlet),
          "--main-outlet-length", std::to_string(mainOutlet),
          "--branch-inlet-length", std::to_string(branchLength)
        }) != 0) {
      throw std::runtime_error("shared MSH2-to-STL exporter failed");
    }
    std::cout << "PASS: native C++ geometry driver generated " << output << '\n';
  } catch (const std::exception& error) {
    std::cerr << "T-junction C++ geometry failed: " << error.what() << '\n';
    return 2;
  }
  return 0;
}
