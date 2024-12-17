#include <iostream>  
#include <mpi.h>  
#include <vector>  

void thomas_algorithm(double* a, double* b, double* c, double* d, double* x, int n) {  
    // ǰ����Ԫ  
    for (int i = 1; i < n; ++i) {  
        double w = a[i - 1] / b[i - 1];  
        b[i] -= w * c[i - 1];  
        d[i] -= w * d[i - 1];  
    }  

    // �������  
    x[n - 1] = d[n - 1] / b[n - 1];  
    for (int i = n - 2; i >= 0; --i) {  
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];  
    }  
}  

int main(int argc, char** argv) {  
    MPI_Init(&argc, &argv);  

    int rank, size;  
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  
    MPI_Comm_size(MPI_COMM_WORLD, &size);  

    const int n = 6; // ��������СΪ6  
    double a[n - 1] = { -1, -1, -1, -1, -1 }; // �¶Խ���  
    double b[n] = {  2,  2,  2,  2,  2,  2 }; // ���Խ���  
    double c[n - 1] = { -1, -1, -1, -1, -1 }; // �϶Խ���  
    double d[n] = { 1, 0, 0, 0, 0, 1 }; // �ұߵĳ�������  
    double x[n]; // �������  

    // �����⻮�ָ���������  
    int local_n = n / size; // ÿ�����̴������Ŀ����������n���Ա�size����  

    // �������ÿ��������ֻ����̶���һ����  
    // ���в��е�ǰ����Ԫ  
    for (int i = 0; i < size; ++i) {  
        if (rank == i) {  
            thomas_algorithm(a + i * local_n, b + i * local_n, c + i * local_n, d + i * local_n, x + i * local_n, local_n);  
        }  
        MPI_Bcast(x, n, MPI_DOUBLE, i, MPI_COMM_WORLD); // �㲥���  
    }  

    if (rank == 0) {  
        std::cout << "Solution: ";  
        for (int i = 0; i < n; ++i) {  
            std::cout << x[i] << " ";  
        }  
        std::cout << std::endl;  
    }  

    MPI_Finalize();  
    return 0;  
}  