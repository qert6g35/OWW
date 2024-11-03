#include <mpi.h>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <thread>
#include <vector>

static const int NUMBER_OF_THREADS = 2;
static const int MATRIX_SIZE = 4;
using type = int16_t;
//using MPI_type = MPI_INT

// using Matrix = std::vector<std::vector<int>>;

void T(type *m){
    //
    // transponowanie macierzy (no odwracamy się ) otrzymujemy m^T
    //
    static type mT[MATRIX_SIZE * MATRIX_SIZE] = {0};
    for (int i = 0; i < MATRIX_SIZE; i++)
        for (int j = 0; j < MATRIX_SIZE; j++){
            mT[i * MATRIX_SIZE + j] = m[j * MATRIX_SIZE + i];
        }
    for (int i = 0; i < MATRIX_SIZE*MATRIX_SIZE; i++)
        m[i] = mT[i];
}

void generate_data(type *m, const type range_from,const type range_to){
    //
    // generujemy dane i wkładamy je do m
    //
    // dane są z rozkładu jednostajnego z zakresu <range_from,range_to)
    //
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<type>  distr(range_from, range_to);
    //type sum;
    for(int i = 0; i<MATRIX_SIZE*MATRIX_SIZE; i++){
        m[i] = distr(generator);
    }
    //std::cout << "sum:" << sum << "   we sum n values, n:"<< MATRIX_SIZE*MATRIX_SIZE <<std::endl;
}

void matmul_single(type *c, const type *a, const type *b) {
    for (int i = 0; i < MATRIX_SIZE; i++)
        for (int j = 0; j < MATRIX_SIZE; j++)
            for (int k = 0; k < MATRIX_SIZE; k++)
                c[i * MATRIX_SIZE + j] += a[i * MATRIX_SIZE + k] * b[k * MATRIX_SIZE + j];
}

void showMatrix(const type *matrix){
    for (int i = 0; i < MATRIX_SIZE; i++){
        for (int j = 0; j < MATRIX_SIZE; j++){
            std::cout<<matrix[i * MATRIX_SIZE + j]<<" ";
        }
        std::cout<< std::endl;
    }
}


// zakładam że ktos transponował B i teraz indeksy idą tak
// B =  0 3 6
//      1 4 7 
//      2 5 8
// 
// A =  0 1 2 
//      3 4 5
//      6 7 8

void matmul_thread_t(type *c, const type *a, const type *b, int thread_id) {
    int all_operations = (MATRIX_SIZE * MATRIX_SIZE);
    int op_for_one = all_operations/NUMBER_OF_THREADS;
    int rest_op = all_operations%NUMBER_OF_THREADS;

    int start_op = op_for_one * thread_id + rest_op * ( thread_id != 0 );
    int finish_op = op_for_one * (thread_id+1) + rest_op + rest_op;

    type elem_of_a, elem_of_b;
    type elem_of_c;
    int row_a, start_a, finish_a;
    int col_b, start_b, finish_b;

    for(int elem_id = start_op; elem_id < finish_op; elem_id++){
        row_a = (elem_id/MATRIX_SIZE);
        //      0 0 0  id/size dla size = 3 macierzy wyjściowej
        //      1 1 1
        //      2 2 2
        start_a = (row_a*MATRIX_SIZE);
        // finish_a = start_a + MATRIX_SIZE;
        col_b = (elem_id%MATRIX_SIZE); // zakładam że transponowano b i ze indeksowanie teraz jest z gory na dol
        //      0 1 2  id%size  dla size = 3 macierzy wyjściowej
        //      0 1 2
        //      0 1 2
        start_b = (col_b*MATRIX_SIZE); 
        // finish_b = start_b + MATRIX_SIZE;
        elem_of_c = 0;
        for(int i = 0; i < MATRIX_SIZE; i++){
            elem_of_a = a[start_a + i];
            elem_of_b = b[start_b + i];
            elem_of_c += elem_of_a * elem_of_b;
        }
        c[elem_id] = elem_of_c;
        
    }
}

// Funkcja inicjująca wątki
void matmul_threaded(type *c, const type *a, const type *b) {
    std::thread threads[NUMBER_OF_THREADS];
    for (int i = 0; i < NUMBER_OF_THREADS; ++i) {
        threads[i] = std::thread(matmul_thread_t, c, a, b, i);
    }

    for (int i = 0; i < NUMBER_OF_THREADS; ++i) {
        threads[i].join();
    }
    
}

void matmul_MPI(type *c, const type *a, const type *b, int start_op,int finish_op) {
    type elem_of_a, elem_of_b;
    type elem_of_c;
    int row_a, start_a, finish_a;
    int col_b, start_b, finish_b;

    for(int elem_id = start_op; elem_id < finish_op; elem_id++){
        row_a = (elem_id/MATRIX_SIZE);
        //      0 0 0  id/size dla size = 3 macierzy wyjściowej
        //      1 1 1
        //      2 2 2
        start_a = (row_a*MATRIX_SIZE);
        // finish_a = start_a + MATRIX_SIZE;
        col_b = (elem_id%MATRIX_SIZE); // zakładam że transponowano b i ze indeksowanie teraz jest z gory na dol
        //      0 1 2  id%size  dla size = 3 macierzy wyjściowej
        //      0 1 2
        //      0 1 2
        start_b = (col_b*MATRIX_SIZE); 
        // finish_b = start_b + MATRIX_SIZE;
        elem_of_c = 0;
        for(int i = 0; i < MATRIX_SIZE; i++){
            elem_of_a = a[start_a + i];
            elem_of_b = b[start_b + i];
            elem_of_c += elem_of_a * elem_of_b;
        }
        c[elem_id] = elem_of_c;
        
    }
}


int main(int argc,char** argv) {
    MPI_Init(&argc,&argv);// inicjalizacji MPI
    type a[MATRIX_SIZE * MATRIX_SIZE];
    type b[MATRIX_SIZE * MATRIX_SIZE];  //=  { 1,0,0,0,
                                       //     0,1,0,0,
                                       //     0,0,1,0,
                                       //     0,0,0,1};
    type c[MATRIX_SIZE * MATRIX_SIZE] = {0};
    
    const type range_start = 0;// wypadalo by to zmienic na wczytywane z arg
    const type range_end = 10;
    
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    if(world_rank == 0){
      generate_data(a,range_start,range_end);
      generate_data(b,range_start,range_end);
    }
    
    MPI_Bcast(a,MATRIX_SIZE * MATRIX_SIZE,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(b,MATRIX_SIZE * MATRIX_SIZE,MPI_INT,0,MPI_COMM_WORLD);
    
    int start_end_operation[2] = {0,0};
    
    const int all_operations = (MATRIX_SIZE * MATRIX_SIZE);
    const int op_for_one = all_operations/world_size;
    const int rest_op = all_operations%world_size;
    if(world_rank == 0){
        std::cout<<"matrix A"<<std::endl;
        showMatrix(a);
        std::cout<<"matrix B"<<std::endl;
        showMatrix(b);
        for(int thread_id = world_size-1; thread_id >= 0; thread_id--){
            start_end_operation[0] = op_for_one * thread_id;
            start_end_operation[1] = op_for_one * (thread_id+1) + rest_op * ( thread_id == world_size-1 );
            
            if(thread_id != 0){
                MPI_Send(&start_end_operation, 2, MPI_INT,thread_id , 0, MPI_COMM_WORLD);
            }
        }
    }else{
        MPI_Recv(&start_end_operation,2,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    
    //std::cout<<"worker_rank:" << world_rank << ", start_op:"<< start_end_operation[0] << ", end_op:"<< start_end_operation[1] <<std::endl;
    matmul_MPI(c, a, b, start_end_operation[0],start_end_operation[1]);

    if(world_rank == 0){
        //showMatrix(c);
        type c_rcv[op_for_one + rest_op] = {0};
        for(int thread_id = world_size-1; thread_id > 0; thread_id--){
            MPI_Recv(&c_rcv,(op_for_one + rest_op* ( thread_id == world_size-1 )),MPI_INT,thread_id,thread_id,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            //std::cout<<"recived data from worker:" << thread_id << " start_sending"<<std::endl;
            for (int i = 0; i< (op_for_one + rest_op* ( thread_id == world_size-1 )); i++){
              c[thread_id * op_for_one + i] = c_rcv[i];
            }
        }
        std::cout<<"Solution"<<std::endl;
        showMatrix(c);
    }else{
        const int num_of_op = start_end_operation[1] - start_end_operation[0];
        type sender[num_of_op] = {0};
        for(int i = 0;i < num_of_op; i++){
          sender[i] = c[i + start_end_operation[0]];
          //std::cout<<"worker_rank:" << world_rank << " calculated position"<< i + start_end_operation[0] << ", as value equal:"<< sender[i] <<std::endl;
        }
        //std::cout<<"worker_rank:" << world_rank << " start_sending"<<std::endl;
        MPI_Send(&sender,num_of_op,MPI_INT, 0,world_rank, MPI_COMM_WORLD);
    }
    
    MPI_Finalize();

    return 0;
}
