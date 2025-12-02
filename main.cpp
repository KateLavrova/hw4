#include <iostream>
#include <vector>
#include <iomanip>
#include <mpi.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // N=8, must be divisible by the number of processes (so that everyone gets the same number of rows)
    const int N = 8;
    if (N % size != 0) {
        if (rank == 0) {
            std::cerr << "Error: N must be divisible by number of processes" << std::endl;
        }
        MPI_Finalize();
        return 1;
    }
    
    const int block_size = N / size;
    
    // creating a custom data type for a matrix block
    MPI_Datatype block_type;
    MPI_Type_vector(block_size, N, N, MPI_DOUBLE, &block_type);
    MPI_Type_commit(&block_type);
    
    // initialization and data distribution
    std::vector<double> local_matrix(block_size * N);
    std::vector<double> vector(N);
    std::vector<double> local_result(block_size, 0.0);
    
    if (rank == 0) {
        std::cout<<"Номер 31 (84): MPI-паттерн: Пользовательские типы данных, тип чисел double"<<std::endl;

        // create a complete matrix and vector on process 0
        std::vector<double> full_matrix(N * N);
        
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                full_matrix[i * N + j] = i * N + j + 1;
            }
        }

        for (int i = 0; i < N; ++i) {
            vector[i] = i + 1;
        }
        
        std::copy(full_matrix.begin(), full_matrix.begin() + block_size * N, local_matrix.begin());
        
        // sending blocks to other processes using a custom type
        for (int proc = 1; proc < size; ++proc) {
            int start_row = proc * block_size;
            MPI_Send(&full_matrix[start_row * N], 1, block_type, proc, 0, MPI_COMM_WORLD);
        }
        
    } else {
        MPI_Recv(local_matrix.data(), block_size * N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    
    // sending the vector to all processes
    MPI_Bcast(vector.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    for (int i = 0; i < block_size; ++i) {
        for (int j = 0; j < N; ++j) {
            local_result[i] += local_matrix[i * N + j] * vector[j];
        }
    }
    
    // collection of results
    std::vector<double> global_result;
    if (rank == 0) {
        global_result.resize(N);
    }
    MPI_Gather(local_result.data(), block_size, MPI_DOUBLE,global_result.data(), block_size, MPI_DOUBLE,0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        std::cout<<std::endl<<"Результаты умножения матрицы "<<N<<"x"<<N<<" на вектор "<<N<<"x1:"<<std::endl;
        for (int i = 0; i < N; ++i) {
            std::cout<<"y["<<i<<"] = "<<global_result[i]<<std::endl;
        }
    }
    
    MPI_Type_free(&block_type);
    MPI_Finalize();
    
    return 0;
}
