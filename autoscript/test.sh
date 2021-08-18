rm test_examples -rf
cp test test_examples -rf
module load hdf5/1.8.22 cgns/4.2.0 metis/5.1.0 mpi/latest
export PATH=${PWD}/oneflow_bin/bin:$PATH
cd test_examples
python3 test.py "mpirun -np 1" "OneFLOW"
