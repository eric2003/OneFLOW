function InstallMETIS() {
    Write-Host "Installing METIS5.1.0..."
    cd lib
    ls
    Write-Host "ls C:/Program Files/7-zip"
    ls "C:/Program Files/7-zip"
    Start-Process "C:/Program Files/7-zip/7z.exe" -Wait -ArgumentList 'x .\metis-5.1.0eric.tar.gz'
    ls
    Start-Process "C:/Program Files/7-zip/7z.exe" -Wait -ArgumentList 'x .\metis-5.1.0eric.tar'
    ls
    cd metis-5.1.0eric
    #md build
    cd build
    $env:Path += ";C:/Program Files/CMake/bin/"
    cmake --version
    cmake -G "Visual Studio 14 2015 Win64" -SC:/projects/OneFLOW/lib/metis-5.1.0eric -BC:/projects/OneFLOW/lib/metis-5.1.0eric/build/ `
    -DMETIS_ENABLE_64BIT="ON" -DCMAKE_INSTALL_PREFIX="C:/Program Files/metis"
    Write-Host "cmake -G finished"
    cmake --build . --target INSTALL --config release
    Write-Host "cmake --build --target INSTALL --config release finished"
    Write-Host "Installing METIS5.1.0 Finished..."
}


function main() {
    InstallMETIS
}

main