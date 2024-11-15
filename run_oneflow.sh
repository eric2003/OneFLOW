prj_dir=${PWD}
cd build
module avail
module purge
module load test/cmake-3.21.1
module load test/metis-5.1.0
module load test/cgns-3.4.0
module load nvhpc/21.7
module load test/nvhpcmpi

root_dir=${PWD}
echo "root_dir=$root_dir"
install_dir=$root_dir/Run
export PATH=$install_dir/bin/:$install_dir/lib:$PATH
export LD_LIBRARY_PATH=$install_dir/lib:$LD_LIBRARY_PATH
mpirun -np 1 OneFLOW
cd $install_dir
rm -rf test
cp $prj_dir/test . -rf
cd ./test
python3 --version
python3 test.py "mpirun -np 1" "OneFLOW"