rm -rf build oneflow_bin
mkdir build oneflow_bin
cd build
module load cmake/3.21.1 hdf5/1.8.22 cgns/4.2.0 metis/5.1.0 mpi/latest qt/6.1.2
cmake .. -DCMAKE_C_COMPILER=${MPI_HOME}/bin/mpicc -DCMAKE_CXX_COMPILER=${MPI_HOME}/bin/mpicxx
cmake --build . --parallel 4 --config release
cmake --install . --prefix ../oneflow_bin
cd ..
