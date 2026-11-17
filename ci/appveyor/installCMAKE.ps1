function InstallCMAKE() {
    md CMake_DIR
    cd CMake_DIR
    # install CMake
    Write-Host "Installing CMake v3.17.0-rc2..." -ForegroundColor Cyan
    appveyor DownloadFile https://github.com/Kitware/CMake/releases/download/v3.17.0-rc2/cmake-3.17.0-rc2-win64-x64.msi
    Start-Process msiexec.exe -Wait -ArgumentList '/I cmake-3.17.0-rc2-win64-x64.msi /quiet /qn /l*v log.txt'

    $env:Path += ";C:/Program Files/CMake/bin/"
    cmake --version
    Write-Host "CMake v3.17.0-rc2 installed" -ForegroundColor Green
}


function main() {
    InstallCMAKE
}

main