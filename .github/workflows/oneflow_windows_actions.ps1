# ============================================================
# OneFLOW Windows dependency versions
# ============================================================

$global:HDF5_VERSION = "1.14.3"
$global:CGNS_VERSION = "4.4.0"
$global:METIS_VERSION = "5.1.0"

# ============================================================
# OneFLOW dependency installation prefixes
# ============================================================

$global:CGNS_PREFIX =
    "$env:GITHUB_WORKSPACE/_deps/cgns/$global:CGNS_VERSION"

$global:METIS_PREFIX =
    "$env:GITHUB_WORKSPACE/_deps/metis/METIS-VS2022-STATIC"

$global:HDF5_PREFIX =
    "$env:GITHUB_WORKSPACE/_deps/hdf5/$global:HDF5_VERSION"
	
# ============================================================
# Environment variable helpers
# ============================================================

function AddMachinePath( $varPath ) {
    $Env:path = GetMachineEnvironmentVariable( "path" )
    $Env:path = $Env:Path + ";$varPath"

    ModifyMachineEnvironmentVariable "path" $Env:path
}


function GetMachineEnvironmentVariable( $varName ) {
    $varValue =
        [System.Environment]::GetEnvironmentvariable(
            $varName,
            [System.EnvironmentVariableTarget]::Machine
        )

    $varValue
}


function ModifyMachineEnvironmentVariable( $varName, $varValue ) {
    $target = "Machine"

    [System.Environment]::SetEnvironmentVariable(
        $varName,
        $varValue,
        $target
    )
}


# ============================================================
# File download helpers
# ============================================================

function MyGetFileName( $filePath ) {
    $file_name_complete =
        [System.IO.Path]::GetFileName("$filePath")

    $file_name_complete

    # Write-Host "fileNameFull :" $file_name_complete
}


function MyDownloadFile( $fullFilePath ) {
    $my_filename = MyGetFileName("$fullFilePath")
    $my_location = Get-Location
    $my_local_filename =
        "$my_location" + "/" + $my_filename

    $my_client = new-object System.Net.WebClient

    $my_client.DownloadFile(
        $fullFilePath,
        $my_local_filename
    )
}


function MyDownloadFile2( $fullFilePath, $my_filename ) {
    Write-Host "MyDownloadFile2 fullFilePath is $fullFilePath"
    Write-Host "MyDownloadFile2 my_filename is $my_filename"

    $my_location = Get-Location

    Write-Host "MyDownloadFile2 my_location is $my_location"

    $my_local_filename =
        "$my_location" + "/" + $my_filename

    Write-Host "MyDownloadFile2 my_local_filename is $my_local_filename"

    $my_client = new-object System.Net.WebClient

    $my_client.DownloadFile(
        $fullFilePath,
        $my_local_filename
    )
}


# ============================================================
# Microsoft MPI
# ============================================================

function InstallMSMPI() {
    # Install MPI SDK and Runtime.
    Write-Host "Installing Microsoft MPI SDK..."

    $download_url =
        "https://download.microsoft.com/download/A/E/0/AE002626-9D9D-448D-8197-1EA510E297CE/"

    $msmpisdk_filename = "msmpisdk.msi"
    $msmpisdk_webfilename =
        $download_url + $msmpisdk_filename

    MyDownloadFile( $msmpisdk_webfilename )

    Start-Process `
        -FilePath msiexec.exe `
        -ArgumentList "/quiet /qn /i msmpisdk.msi" `
        -Wait

    Write-Host "Microsoft MPI SDK installation complete"

    Write-Host "Installing Microsoft MPI Runtime..."

    $msmpisetup_filename = "msmpisetup.exe"
    $msmpisetup_webfilename =
        $download_url + $msmpisetup_filename

    MyDownloadFile( $msmpisetup_webfilename )

    Start-Process `
        -FilePath MSMpiSetup.exe `
        -ArgumentList -unattend `
        -Wait

    Write-Host "Microsoft MPI Runtime installation complete..."

    $msmpi_bin_path =
        "C:/Program Files/Microsoft MPI/Bin"

    $msmpi_sdk_path =
        "C:/Program Files (x86)/Microsoft SDKs/MPI"

    AddMachinePath( $msmpi_bin_path )

    # Update the current PowerShell process environment.
    $Env:Path =
        "$msmpi_bin_path;$Env:Path"

    if ( Test-Path "$msmpi_bin_path/mpiexec.exe" ) {
        Write-Host "MPI executable: FOUND"
    }
    else {
        Write-Error "MPI executable: NOT FOUND"
        exit 1
    }

    # Persist MPI path for subsequent GitHub Actions steps.
    $msmpi_bin_path |
        Out-File `
            -FilePath $env:GITHUB_PATH `
            -Encoding utf8 `
            -Append

    Write-Host "ls $msmpi_sdk_path"
    ls $msmpi_sdk_path

    Write-Host "ls $msmpi_bin_path"
    ls $msmpi_bin_path

    Write-Host "Resolved mpiexec:"
    Get-Command mpiexec

    Write-Host "Testing mpiexec:"
    mpiexec -help
}


# ============================================================
# HDF5
# ============================================================

function DownloadHDF5() {
    $global:hdf5_major = 1
    $global:hdf5_minor = 14
    $global:hdf5_patch = 3

    $hdf5_main_name = "hdf5"

    $hdf5_dir1 =
        "$hdf5_main_name-$hdf5_major.$hdf5_minor"

    $hdf5_dir2 =
        "$hdf5_dir1.$hdf5_patch"

    $global:hdf5_version_name =
        "$hdf5_dir2"

    Write-Host "Downloading $hdf5_version_name..."

    $hdf5_url =
        "https://support.hdfgroup.org/ftp/HDF5/releases/$hdf5_dir1/$hdf5_dir2/bin/windows/"

    $hdf5_name =
        "$hdf5_version_name-Std-win10_64-vs17"

    $global:hdf5_package_name =
        $hdf5_name + ".zip"

    $hdf5_web_addr =
        $hdf5_url + $hdf5_package_name

    MyDownloadFile( $hdf5_web_addr )

    ls

    Write-Host "$hdf5_version_name downloading complete"
}


function InstallHDF5() {
    DownloadHDF5

    $hdf5_version_name_upper =
        $hdf5_version_name.ToUpper()

    $zipexe =
        "C:/Program Files/7-zip/7z.exe"

    $arg =
        "x ./$global:hdf5_package_name"

    Start-Process `
        $zipexe `
        -Wait `
        -ArgumentList $arg

    ls

    cd hdf

    ls

    Write-Host "Installing $hdf5_version_name_upper..."

    $arg1 =
        "/quiet /qn /i $hdf5_version_name_upper-win64.msi"

    Start-Process `
        -FilePath msiexec.exe `
        -ArgumentList $arg1 `
        -Wait

    Write-Host "$hdf5_version_name_upper installation complete"

    $HDF5_InstallDir =
        "C:/Program Files/HDF_Group/HDF5/$hdf5_major.$hdf5_minor.$hdf5_patch"

    Write-Host "ls $HDF5_InstallDir"

    ls $HDF5_InstallDir

    $hdf5_dir_varName =
        "HDF5_DIR"

    $hdf5_dir_varValue =
        "$HDF5_InstallDir/cmake"

    Write-Host "hdf5_dir_varName = $hdf5_dir_varName"
    Write-Host "hdf5_dir_varValue = $hdf5_dir_varValue"

    ModifyMachineEnvironmentVariable `
        $hdf5_dir_varName `
        $hdf5_dir_varValue

    Write-Host "Checking Env:HDF5_DIR..."
    Write-Host "Env:HDF5_DIR = $Env:HDF5_DIR"

    $tmp =
        GetMachineEnvironmentVariable("HDF5_DIR")

    Write-Host "tmp = $tmp"

    $tmp1 =
        GetMachineEnvironmentVariable("$hdf5_dir_varName")

    Write-Host "tmp1 = $tmp1"

    cd ..

    pwd

    Write-Host "$hdf5_version_name_upper installation complete..."
}


# ============================================================
# CGNS
# ============================================================

function InstallCGNS() {
    $global:cgns_version =
        $global:CGNS_VERSION

    $cgns_prefix =
        $global:CGNS_PREFIX

    # Check whether the final CGNS installation already exists.
    if ( (Test-Path "$cgns_prefix/include/cgnslib.h") -and
         (Test-Path "$cgns_prefix/lib/cgnsdll.lib") ) {

        Write-Host "===== CGNS cache hit ====="
        Write-Host "CGNS installation already exists:"
        Write-Host "$cgns_prefix"

        Write-Host "CGNS include directory:"
        ls "$cgns_prefix/include"

        Write-Host "CGNS library directory:"
        ls "$cgns_prefix/lib"

        return
    }

    Write-Host "===== CGNS cache miss ====="
    Write-Host "Installing CGNS-$cgns_version..."

    DownloadCGNS

    Write-Host "Installing CGNS..."

    $zipexe =
        "C:/Program Files/7-zip/7z.exe"

    $arg =
        "x ./CGNS-$cgns_version.zip"

    Start-Process `
        $zipexe `
        -Wait `
        -ArgumentList $arg

    cd CGNS-$cgns_version

    mkdir build

    cd build

    $tmp =
        GetMachineEnvironmentVariable("HDF5_DIR")

    $Env:HDF5_DIR =
        $tmp

    cmake `
        -DCGNS_ENABLE_64BIT="ON" `
        -DCGNS_ENABLE_HDF5="ON" `
        -DCGNS_BUILD_SHARED="ON" `
        ../

    cmake --build . --parallel 4 --config release

    cmake --install . --prefix $cgns_prefix

    Write-Host "CGNS installation directory:"
    Write-Host "$cgns_prefix"

    Write-Host "CGNS include directory:"
    Write-Host "$cgns_prefix/include"

    ls "$cgns_prefix/include"

    cd ../../

    Write-Host "CGNS-$cgns_version installation complete..."

    # Verify the final installation.
    if ( Test-Path "$cgns_prefix/include/cgnslib.h" ) {
        Write-Host "CGNS header: FOUND"
    }
    else {
        Write-Error "CGNS header: NOT FOUND"
        exit 1
    }

    if ( Test-Path "$cgns_prefix/lib/cgnsdll.lib" ) {
        Write-Host "CGNS library: FOUND"
    }
    else {
        Write-Error "CGNS library: NOT FOUND"
        exit 1
    }
}


# ============================================================
# METIS
# ============================================================

function DownloadMETIS() {
    Write-Host "Downloading METIS-$global:METIS_VERSION..."

    git --version

    $metis_project_url =
        "https://github.com/eric2003/"

    $metis_project_name =
        "METIS-$global:METIS_VERSION-Modified"

    $metis_project_git_name =
        $metis_project_name + ".git"

    $metis_project_web_addr =
        $metis_project_url + $metis_project_git_name

    git clone $metis_project_web_addr

    ls

    cd $metis_project_name

    ls

    Write-Host "Downloading METIS-$global:METIS_VERSION complete..."
}


function InstallMETIS() {
    $metis_prefix =
        $global:METIS_PREFIX

    # Check whether the final METIS installation already exists.
    if ( (Test-Path "$metis_prefix/include/metis.h") -and
         (Test-Path "$metis_prefix/lib/metis.lib") ) {

        Write-Host "===== METIS cache hit ====="
        Write-Host "METIS installation already exists:"
        Write-Host "$metis_prefix"

        return
    }

    Write-Host "===== METIS cache miss ====="
    Write-Host "Installing METIS-$global:METIS_VERSION..."

    DownloadMETIS

    Write-Host "mkdir build..."

    mkdir build

    Write-Host "ls..."

    ls

    cd build

    cmake ../

    cmake --build . --parallel 4 --config release

    cmake --install . --prefix $metis_prefix

    cd ../../

    pwd

    Write-Host "METIS-$global:METIS_VERSION installation complete..."

    # Verify the final installation.
    if ( Test-Path "$metis_prefix/include/metis.h" ) {
        Write-Host "METIS header: FOUND"
    }
    else {
        Write-Error "METIS header: NOT FOUND"
        exit 1
    }

    if ( Test-Path "$metis_prefix/lib/metis.lib" ) {
        Write-Host "METIS library: FOUND"
    }
    else {
        Write-Error "METIS library: NOT FOUND"
        exit 1
    }
}


# ============================================================
# CGNS download
# ============================================================

function DownloadCGNS() {
    Write-Host "Downloading CGNS-$global:CGNS_VERSION..."

    $download_url =
        "https://github.com/CGNS/CGNS/archive/refs/tags/"

    $cgns_filename =
        "v$global:CGNS_VERSION.zip"

    $cgns_real_filename =
        "CGNS-$global:CGNS_VERSION.zip"

    $cgns_webfilename =
        $download_url + $cgns_filename

    Write-Host "download_url is $download_url"
    Write-Host "cgns_webfilename is $cgns_webfilename"
    Write-Host "cgns_real_filename is $cgns_real_filename"

    Write-Host "calling MyDownloadFile2..."

    MyDownloadFile2 `
        $cgns_webfilename `
        $cgns_real_filename

    ls

    Write-Host "CGNS-$global:CGNS_VERSION downloading complete"
}


# ============================================================
# Download directory management
# ============================================================

function InitDownload() {
    mkdir download

    cd download
}


function ExitDownload() {
    cd ..
}


# ============================================================
# OneFLOW build
# ============================================================

function CompileOneFLOW() {
    Write-Host "Compile OneFLOW ..."

    mkdir build

    cd build

    $oneflow_prefix =
        "$env:GITHUB_WORKSPACE/install"

    $metis_root =
        $global:METIS_PREFIX

    $cgns_root =
        $global:CGNS_PREFIX

    Write-Host "METIS_ROOT = $metis_root"
    Write-Host "CGNS_ROOT  = $cgns_root"

    Write-Host "Checking CGNS header:"

    ls "$cgns_root/include"
    ls "$cgns_root/include/cgnslib.h"

    Write-Host "===== Dependency Check ====="

    Write-Host "CGNS_ROOT = $cgns_root"
    Write-Host "CGNS include = $cgns_root/include"
    Write-Host "CGNS library = $cgns_root/lib/cgnsdll.lib"

    if ( Test-Path "$cgns_root/include/cgnslib.h" ) {
        Write-Host "CGNS header: FOUND"
    }
    else {
        Write-Error "CGNS header: NOT FOUND"
        exit 1
    }

    if ( Test-Path "$cgns_root/lib/cgnsdll.lib" ) {
        Write-Host "CGNS library: FOUND"
    }
    else {
        Write-Error "CGNS library: NOT FOUND"
        exit 1
    }

    if ( Test-Path "$metis_root/include/metis.h" ) {
        Write-Host "METIS header: FOUND"
    }
    else {
        Write-Error "METIS header: NOT FOUND"
        exit 1
    }

    if ( Test-Path "$metis_root/lib/metis.lib" ) {
        Write-Host "METIS library: FOUND"
    }
    else {
        Write-Error "METIS library: NOT FOUND"
        exit 1
    }

    $start = Get-Date

    # Configure
    cmake `
        -DMETIS_ROOT="$metis_root" `
        -DCGNS_ROOT="$cgns_root" `
        ../

    Write-Host "===== CMake Configure: $((Get-Date) - $start) ====="

    $start = Get-Date

    # Build
    cmake --build . --parallel 8 --config release

    Write-Host "===== CMake Build: $((Get-Date) - $start) ====="

    $start = Get-Date

    # Install
    cmake --install . --prefix $oneflow_prefix

    Write-Host "===== CMake Install: $((Get-Date) - $start) ====="

    Write-Host "Compile OneFLOW complete..."
}


# ============================================================
# Main
# ============================================================

function main() {
    InitDownload

    $start = Get-Date

    InstallMSMPI

    Write-Host "===== InstallMSMPI: $((Get-Date) - $start) ====="

    $start = Get-Date

    InstallHDF5

    Write-Host "===== InstallHDF5: $((Get-Date) - $start) ====="

    $start = Get-Date

    InstallCGNS

    Write-Host "===== InstallCGNS: $((Get-Date) - $start) ====="

    $start = Get-Date

    InstallMETIS

    Write-Host "===== InstallMETIS: $((Get-Date) - $start) ====="

    ExitDownload

    $start = Get-Date

    CompileOneFLOW

    Write-Host "===== CompileOneFLOW: $((Get-Date) - $start) ====="
}


main