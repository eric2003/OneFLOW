function InstallMPI() {
    md mpi
    cd mpi
    # install MPI SDK and Runtime
    Write-Host "Installing Microsoft MPI SDK..."
    appveyor DownloadFile https://download.microsoft.com/download/A/E/0/AE002626-9D9D-448D-8197-1EA510E297CE/msmpisdk.msi
    Start-Process -FilePath msiexec.exe -ArgumentList "/quiet /qn /i msmpisdk.msi" -Wait
    Write-Host "Microsoft MPI SDK installation complete"
    Write-Host "Installing Microsoft MPI Runtime..."
    appveyor DownloadFile https://download.microsoft.com/download/A/E/0/AE002626-9D9D-448D-8197-1EA510E297CE/msmpisetup.exe
    Start-Process -FilePath MSMpiSetup.exe -ArgumentList -unattend -Wait
    Write-Host "Microsoft MPI Runtime installation complete..."
    cd 
    dir
    c:
    cd "C:/Program Files (x86)/Microsoft SDKs/MPI"
    dir
}


function main() {
    InstallMPI
}

main