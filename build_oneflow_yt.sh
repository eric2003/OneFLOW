prj_dir=${PWD}
mkdir build
cd build
rm -rf *
module avail
module purge
module load gnu8/8.3.0
module load nvhpc/21.7
module load cmake/gcc-8.3.0/3.21.1
module load metis/nvhpc-21.7/5.1.0
module load hdf5/nvhpc-21.7/1.8.22
module load cgns/nvhpc-21.7/4.2.0
pwd
root_dir=${PWD}
root_dir=${PWD}
echo "root_dir=$root_dir"
install_dir=$root_dir/Run
cmake -DCMAKE_INSTALL_PREFIX=$install_dir ..
#cmake --build . --parallel 4  > log.txt 2>&1
cmake  --build . --parallel 4
make install
export PATH=$install_dir/bin/:$install_dir/lib:$PATH
export LD_LIBRARY_PATH=$install_dir/lib:$LD_LIBRARY_PATH
mpirun -np 1 OneFLOW
cd $install_dir
rm -rf test
cp $prj_dir/test . -rf
cd ./test
python3 --version
python3 test.py "mpirun -np 1" "OneFLOW"
