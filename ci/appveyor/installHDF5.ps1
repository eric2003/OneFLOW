function InstallHDF5() {
    md HDF5_DIR
    cd HDF5_DIR
    Write-Host "Installing HDF5..."
    Write-Host "Installing HDF5 Finished..."
}


function main() {
    InstallHDF5
}

main