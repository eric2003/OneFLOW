0.
cd d:/github/OneFLOW/example/poisson/gauss_seidel/02/build/

1.
cmake -DCMAKE_TOOLCHAIN_FILE="c:/dev/vcpkg/scripts/buildsystems/vcpkg.cmake"  -DCMAKE_BUILD_TYPE=Release ../

2.
cmake --build . --parallel 4 --config Release

3.
 ./Release/testprj.exe

----------------------------------------
Remove-Item * -Recurse -Force
