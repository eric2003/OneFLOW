function InstallCMAKE() {
    Write-Host "Installing CMake 3.4.0 ..." -ForegroundColor Cyan
    $exePath = "$($env:USERPROFILE)\cmake-3.4.0-rc2-win32-x86.exe"
    Write-Host "Downloading..."
    (New-Object Net.WebClient).DownloadFile('https://cmake.org/files/v3.4/cmake-3.4.0-rc2-win32-x86.exe', $exePath)
    Write-Host "Installing..."
    cmd /c start /wait $exePath /S
    cmake --version
    Write-Host "CMake 3.4.0 installed" -ForegroundColor Green
}


function main() {
    InstallCMAKE
}

main