PS D:\github\OneFLOW\example\eno\enolinear\fortran\01\build> cmd.exe "/K" '"C:\Program Files (x86)\Intel\oneAPI\setvars.bat" && powershell'
ifort ../enolinear.for -o enolinear
cmake ../ -T fortran=ifx