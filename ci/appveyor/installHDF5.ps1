function InstallHDF5() {
    Write-Host "Installing HDF5..."
    cd lib
    ls
    Write-Host "ls C:/Program Files/7-zip"
    ls "C:/Program Files/7-zip"
    Start-Process "C:/Program Files/7-zip/7z.exe" -Wait -ArgumentList 'x .\hdf5-1.10.5.tar.gz'
    ls
    Start-Process "C:/Program Files/7-zip/7z.exe" -Wait -ArgumentList 'x .\hdf5-1.10.5.tar'
    ls
    cd hdf5-1.10.5
    md build
    $env:Path += ";C:/Program Files/CMake/bin/"
    cmake --version
    cmake -G "Visual Studio 14 2015 Win64" -SC:/projects/OneFLOW/lib/hdf5-1.10.5 -BC:/projects/OneFLOW/lib/hdf5-1.10.5/build/ `
    -DBUILD_SHARED_LIBS="OFF"
    Write-Host "cmake -G finished"
    cmake --build . --target INSTALL --config release
    Write-Host "cmake --build --target INSTALL --config release finished"
    Write-Host "Installing HDF5 Finished..."
}


function main() {
    InstallHDF5
}

main