function InstallCGNS() {
    Write-Host "Installing CGNS..."
    cd lib
    ls
    Write-Host "ls C:/Program Files/7-zip"
    ls "C:/Program Files/7-zip"
    Start-Process "C:/Program Files/7-zip/7z.exe" -Wait -ArgumentList 'x .\cgns-3.4.0.tar.gz'
    ls
    Start-Process "C:/Program Files/7-zip/7z.exe" -Wait -ArgumentList 'x .\cgns-3.4.0.tar'
    ls
    cd cgns-3.4.0
    md build
    cd build
    $env:Path += ";C:/Program Files/CMake/bin/"
    cmake --version
    $env:HDF5_DIR="C:/Program Files/HDF_Group/HDF5/1.10.5/share/cmake"
    $env:HDF5_DIR
    cmake -G "Visual Studio 14 2015 Win64" -SC:/projects/OneFLOW/lib/cgns-3.4.0 -BC:/projects/OneFLOW/lib/cgns-3.4.0/build/ `
    -DCGNS_ENABLE_64BIT="ON" -DCGNS_ENABLE_HDF5="ON" -DCGNS_BUILD_SHARED="OFF"

    Write-Host "cmake -G finished"
    cmake --build . --target INSTALL --config release
    Write-Host "cmake --build --target INSTALL --config release finished"
    Write-Host "Installing CGNS Finished..."
}


function main() {
    InstallCGNS
}

main