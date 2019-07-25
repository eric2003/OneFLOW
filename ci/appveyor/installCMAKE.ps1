function InstallCMAKEOld() {
    Write-Host "Installing CMake 3.4.0 ..." -ForegroundColor Cyan
    $exePath = "$($env:USERPROFILE)\cmake-3.4.0-rc2-win32-x86.exe"
    Write-Host "Downloading..."
    (New-Object Net.WebClient).DownloadFile('https://cmake.org/files/v3.4/cmake-3.4.0-rc2-win32-x86.exe', $exePath)
    Write-Host "Installing..."
    cmd /c start /wait $exePath /S
    cmake --version
    Write-Host "CMake 3.4.0 installed" -ForegroundColor Green
}

function InstallCMAKE() {
    md CMake_DIR
    cd CMake_DIR
    # install CMake
    Write-Host "Installing CMake 3.15.0..." -ForegroundColor Cyan
    appveyor DownloadFile https://github.com/Kitware/CMake/releases/download/v3.15.0/cmake-3.15.0-win64-x64.msi
    #Start-Process -FilePath msiexec.exe -ArgumentList "/quiet /qn /i cmake-3.15.0-win64-x64.msi" -Wait
    Start-Process msiexec.exe -Wait -ArgumentList '/I cmake-3.15.0-win64-x64.msi /quiet /qn /l*v log.txt'
    cmake --version
    Write-Host "CMake 3.15.0 installed" -ForegroundColor Green
}


function main() {
    InstallCMAKE
}

main