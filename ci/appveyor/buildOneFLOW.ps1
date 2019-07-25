function BuildOneFLOW() {
    # Go to build directory and build
    # Build OneFLOW with cmake
    md build
    cd build
    $env:Path += ";C:/Program Files/CMake/bin/"
    cmake --version
    Write-Host "cmake --help..."
    cmake --help
    cmake -G "Visual Studio 14 2015 Win64" -Se:C:/projects/OneFLOW/ -Be:C:/projects/OneFLOW/build/
    cmake --build . --target INSTALL --config release
}

function main() {
    BuildOneFLOW
}

main