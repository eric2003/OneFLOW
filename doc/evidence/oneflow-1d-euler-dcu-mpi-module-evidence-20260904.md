# OneFLOW 本地 MPI module 构建证据

更新时间：2026-09-04

## Module 探测

- 当前 module 系统：Modules 4.8.0
- 当前默认状态：无 module loaded
- 可用候选：`openmpi/4.1.4`、`mpi/latest`、`oneflow/nvhpcmpi`
- 本轮采用：`openmpi/4.1.4`

该 module 提供：

- `mpicxx`、`mpicc`、`mpirun`
- `mpi.h`
- `libmpi.so`
- OpenMPI 4.1.4 wrapper compile/link flags

## 主工程验证

在临时 shell 中执行 `module load openmpi/4.1.4`，并使用 `CXX=mpicxx`；同时显式传入本地 CGNS/METIS 库路径：

```text
MPI_HOME_INC=$MPI_HOME/include
MPI_HOME_LIB=$MPI_HOME/lib/libmpi.so
CGNS_HOME_INC=/usr/include
CGNS_HOME_LIB=/usr/lib/x86_64-linux-gnu/libcgns.so
METIS_HOME_INC=/usr/include
METIS_HOME_LIB=/usr/lib/x86_64-linux-gnu/libmetis.so
```

结果：

- `OneFLOW` 主目标编译：100% 成功
- `hello_test`：构建成功
- `oneflow_euler_backend_cpu_test`：构建成功
- CTest：7/7 通过，0 failed，0 `NOT_BUILT`
- 新增 `EulerDomainStateRegistry.cpp`：已进入并编译于主 `OneFLOW` 目标

## 边界

- module 只解决本地 CPU/MPI 构建环境，不提供 DCU 运行证据。
- 本轮没有改变持久 shell 配置，也没有安装 MPI。
- 当前新增 domain contract/state registry 尚未接入主 solver 时间推进和 CTest 注册。
