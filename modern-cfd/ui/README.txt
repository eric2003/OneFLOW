cmake  ../ -D CMAKE_TOOLCHAIN_FILE="c:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake"

The package python3 is compatible with built-in CMake targets:

    find_package(Python3 COMPONENTS Development REQUIRED)
    target_link_libraries(main PRIVATE Python3::Python)

The package python3 provides a python interpreter that supports virtual environments:

    >tools\python3\python3.10 -m venv c:\path\to\venv
    >set VIRTUAL_ENV=c:\path\to\venv
    >set PATH=c:\path\to\venv\bin;%PATH%
    >set PYTHONHOME=

    See https://docs.python.org/3/library/venv.html for more details.