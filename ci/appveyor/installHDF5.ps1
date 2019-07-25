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
    Write-Host "Installing HDF5 Finished..."
}


function main() {
    InstallHDF5
}

main