prj_dir=${PWD}
cd build
rm -rf *
module avail
module purge
module load test/cmake-3.21.1
module load test/metis-5.1.0
module load test/cgns-3.4.0
module load nvhpc/21.7
module load test/nvhpcmpi
pwd
root_dir=${PWD}
echo "root_dir=$root_dir"
install_dir=$root_dir/Run
cmake -DCMAKE_INSTALL_PREFIX=$install_dir ..
#cmake --build . --parallel 4  > log.txt 2>&1
cmake  --build . --parallel 4
make install
export PATH=$install_dir/bin/:$install_dir/lib:$PATH
mpirun -np 1 OneFLOW
cd $install_dir
rm -rf test
cp $prj_dir/test . -rf
cd ./test
python3 --version
python3 test.py "mpirun -np 1" "OneFLOW"