function AddMachinePath( $varPath ) {
    $Env:path = GetMachineEnvironmentVariable( "path" )
    $Env:path = $Env:Path + ";$varPath" 
    ModifyMachineEnvironmentVariable "path" $Env:path
}

function GetMachineEnvironmentVariable( $varName ) {
    $varValue = [System.Environment]::GetEnvironmentvariable($varName , [System.EnvironmentVariableTarget]::Machine)
    $varValue
}

function ModifyMachineEnvironmentVariable( $varName, $varValue ) {
    $target = "Machine"
    [System.Environment]::SetEnvironmentVariable($varName, $varValue, $target)
}

function MyGetFileName( $filePath ) {
    $file_name_complete = [System.IO.Path]::GetFileName("$filePath")
    $file_name_complete
    #Write-Host "fileNameFull :" $file_name_complete    
}

function MyDownloadFile( $fullFilePath ) {
    $my_filename = MyGetFileName("$fullFilePath")
    $my_location = Get-Location
    $my_local_filename = "$my_location" + "/" + $my_filename
    
    $my_client = new-object System.Net.WebClient
    $my_client.DownloadFile( $fullFilePath, $my_local_filename )    
}

function MyDownloadFile2( $fullFilePath, $my_filename ) {
    Write-Host "MyDownloadFile2 fullFilePath is $fullFilePath"
    Write-Host "MyDownloadFile2 my_filename is $my_filename"
    $my_location = Get-Location
    Write-Host "MyDownloadFile2 my_location is $my_location"
    $my_local_filename = "$my_location" + "/" + $my_filename
    Write-Host "MyDownloadFile2 my_local_filename is $my_local_filename"
    
    $my_client = new-object System.Net.WebClient
    $my_client.DownloadFile( $fullFilePath, $my_local_filename )    
}

function InstallMSMPI() {
    # install MPI SDK and Runtime
    Write-Host "Installing Microsoft MPI SDK..."
    $download_url = "https://download.microsoft.com/download/A/E/0/AE002626-9D9D-448D-8197-1EA510E297CE/"
    $msmpisdk_filename = "msmpisdk.msi"
    $msmpisdk_webfilename = $download_url + $msmpisdk_filename

    MyDownloadFile( $msmpisdk_webfilename )
    
    Start-Process -FilePath msiexec.exe -ArgumentList "/quiet /qn /i msmpisdk.msi" -Wait
    Write-Host "Microsoft MPI SDK installation complete"
    
    Write-Host "Installing Microsoft MPI Runtime..."
    $msmpisetup_filename = "msmpisetup.exe"
    $msmpisetup_webfilename = $download_url + $msmpisetup_filename
    
    MyDownloadFile( $msmpisetup_webfilename )
    
    Start-Process -FilePath MSMpiSetup.exe -ArgumentList -unattend -Wait
    Write-Host "Microsoft MPI Runtime installation complete..."
    
    $msmpi_bin_path = "C:/Program Files/Microsoft MPI/Bin"
    $msmpi_sdk_path = "C:/Program Files (x86)/Microsoft SDKs/MPI"
    
    AddMachinePath( $msmpi_bin_path )
    
    Write-Host "ls $msmpi_sdk_path"
    ls $msmpi_sdk_path
    Write-Host "ls $msmpi_sdk_path"
    ls $msmpi_bin_path
    
    Write-Host "mpiexec"
    mpiexec 
}

function DownloadHDF5() {
    $global:hdf5_major = 1
    $global:hdf5_minor = 14
    $global:hdf5_patch = 3
    $hdf5_main_name = "hdf5"
    $hdf5_dir1 = "$hdf5_main_name-$hdf5_major.$hdf5_minor"
    $hdf5_dir2 = "$hdf5_dir1.$hdf5_patch"
    $global:hdf5_version_name = "$hdf5_dir2"
    Write-Host "Downloading $hdf5_version_name..."
    $hdf5_url = "https://support.hdfgroup.org/ftp/HDF5/releases/$hdf5_dir1/$hdf5_dir2/bin/windows/"
    $hdf5_name         = "$hdf5_version_name-Std-win10_64-vs17"
    $global:hdf5_package_name = $hdf5_name + ".zip"
    $hdf5_web_addr = $hdf5_url + $hdf5_package_name
          
    MyDownloadFile( $hdf5_web_addr )
    ls
    Write-Host "$hdf5_version_name downloading complete"
}

function InstallHDF5() {
    DownloadHDF5
    $hdf5_version_name_upper = $hdf5_version_name.ToUpper()
    $zipexe = "C:/Program Files/7-zip/7z.exe" 
    $arg = "x ./$global:hdf5_package_name"
    Start-Process $zipexe -Wait -ArgumentList $arg
    ls
    cd hdf
    ls
    Write-Host "Installing $hdf5_version_name_upper..."
    $arg1 = "/quiet /qn /i $hdf5_version_name_upper-win64.msi"
    Start-Process -FilePath msiexec.exe -ArgumentList $arg1 -Wait
    Write-Host "$hdf5_version_name_upper installation complete"
    $HDF5_InstallDir = "C:/Program Files/HDF_Group/HDF5/$hdf5_major.$hdf5_minor.$hdf5_patch"
    Write-Host "ls $HDF5_InstallDir"
    ls $HDF5_InstallDir
    
    $hdf5_dir_varName = "HDF5_DIR"
    $hdf5_dir_varValue = "$HDF5_InstallDir/cmake"
    Write-Host "hdf5_dir_varName = $hdf5_dir_varName"
    Write-Host "hdf5_dir_varValue = $hdf5_dir_varValue"
    ModifyMachineEnvironmentVariable $hdf5_dir_varName $hdf5_dir_varValue
    
    Write-Host "Checking Env:HDF5_DIR..."
    Write-Host "Env:HDF5_DIR = $Env:HDF5_DIR"
    $tmp = GetMachineEnvironmentVariable("HDF5_DIR")
    Write-Host "tmp = $tmp"
    $tmp1 = GetMachineEnvironmentVariable("$hdf5_dir_varName")
    Write-Host "tmp1 = $tmp1"


    cd ..
    pwd
    
    Write-Host "$hdf5_version_name_upper installation complete..."
}

function InstallCGNS() {
    $global:cgns_version = "4.4.0"
    DownloadCGNS    
    Write-Host "Installing CGNS..."
    $zipexe = "C:/Program Files/7-zip/7z.exe" 
    $arg = "x ./CGNS-$cgns_version.zip"
    Start-Process $zipexe -Wait -ArgumentList $arg
    ls
    cd CGNS-$cgns_version
    ls
    Write-Host "mkdir build..."
    mkdir build
    Write-Host "ls..."
    ls
    cd build
    $tmp = GetMachineEnvironmentVariable("HDF5_DIR")
    Write-Host "Machine Environment HDF5_DIR = $tmp"
    Write-Host "local Env:HDF5_DIR = $Env:HDF5_DIR"
    $Env:HDF5_DIR = $tmp;
    Write-Host "now Env:HDF5_DIR = $Env:HDF5_DIR"
    $cgns_prefix = "C:/dev/cgns/$cgns_version"
    $cgns_bin = $cgns_prefix + "/bin"
    cmake -DCGNS_ENABLE_64BIT="ON" `
          -DCGNS_ENABLE_HDF5="ON" `
          -DCGNS_BUILD_SHARED="ON" `
          ../
    #      -DCMAKE_INSTALL_PREFIX="C:/dev/cgns" ../
    cmake --build . --parallel 4 --config release
    cmake --install . --prefix $cgns_prefix
    AddMachinePath( $cgns_bin )
    cd ../../
    pwd
    Write-Host "CGNS-$cgns_version installation complete..."
}

function DownloadMETIS() {
    Write-Host "Downloading METIS-5.1.0..."
    git --version
    $metis_project_url      = "https://github.com/eric2003/"
    $metis_project_name     = "METIS-5.1.0-Modified"
    $metis_project_git_name = $metis_project_name + ".git"
    $metis_project_web_addr = $metis_project_url + $metis_project_git_name
    git clone $metis_project_web_addr
    ls
    cd $metis_project_name
    ls
    Write-Host "Downloading METIS-5.1.0 complete..."
}

function InstallMETIS() {
    DownloadMETIS
    Write-Host "Installing METIS-5.1.0..."
    Write-Host "mkdir build..."
    mkdir build
    Write-Host "ls..."
    ls
    cd build
    $metis_prefix = "C:/dev/METIS/METIS-VS2022-STATIC/"
    cmake ../
    cmake --build . --parallel 4 --config release
    cmake --install . --prefix $metis_prefix
    cd ../../
    pwd
    Write-Host "METIS-5.1.0 installation complete..."
}

function DownloadCGNS() {
    Write-Host "Downloading CGNS-$cgns_version..."
    $download_url = "https://github.com/CGNS/CGNS/archive/refs/tags/"
    $cgns_filename = "v$cgns_version.zip"
    $cgns_real_filename = "CGNS-$cgns_version.zip"  
    $cgns_webfilename = $download_url + $cgns_filename
    
    Write-Host "download_url is $download_url"
    Write-Host "cgns_webfilename is $cgns_webfilename"
    Write-Host "cgns_real_filename is $cgns_real_filename"
    Write-Host "calling MyDownloadFile2 cgns_webfilenamec cgns_real_filename"

    #MyDownloadFile( $cgns_webfilename )
    MyDownloadFile2 $cgns_webfilename $cgns_real_filename
    ls
    Write-Host "CGNS-$cgns_version downloading complete"
}

function InitDownload() {
    mkdir download
    cd download
}

function ExitDownload() {
    cd ..
}

function CompileOneFLOW() {
    Write-Host "Compile OneFLOW ..."
    mkdir build
    cd build
    cmake ../
    $oneflow_prefix = "C:/OneFLOW"
    $oneflow_bin = $oneflow_prefix + "/bin"
    cmake --build . --parallel 4 --config release
    cmake --install . --prefix $oneflow_prefix
    AddMachinePath( $oneflow_bin )
    Write-Host "Compile OneFLOW complete..."
}

function main() {
    InitDownload
    InstallMSMPI
    InstallHDF5
    InstallCGNS
    InstallMETIS
    ExitDownload
    CompileOneFLOW
}

main
