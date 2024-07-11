cmake  ../ -D CMAKE_TOOLCHAIN_FILE="c:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake"
cmake .. -D CMAKE_PREFIX_PATH:PATH=C:/local/Qt/6.3.2/msvc2019_64/
cmake .. -D CMAKE_PREFIX_PATH:PATH=C:/local/Qt/Qt6.4.0/6.4.0/msvc2019_64/
cmake .. -D CMAKE_PREFIX_PATH:PATH=C:/local/Qt/6.7.2/msvc2019_64/
c:\local\Qt\6.7.2\
c:\local\Qt\Qt6.4.0\
cmake --build .
C:\local\Qt\6.3.2\msvc2019_64\bin\windeployqt.exe .\Debug\testprj.exe
.\Debug\testprj.exe

cmake  ../ -D CMAKE_TOOLCHAIN_FILE="c:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake"

c:\dev\vcpkg\installed\x64-windows\tools\Qt6\bin\windeployqt.debug.bat `
d:\github\OneFLOW\modern-cfd\build\bin\OneFLOW-UI.exe

C:/dev/vcpkg/installed/x64-windows/tools/Qt6/bin/windeployqt.exe `
d:\github\OneFLOW\modern-cfd\build\bin\OneFLOW-UI.exe