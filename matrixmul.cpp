#include <mpi.h>
#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <thread>
#include <vector>

static const bool MANUAL_TEST = false;
static const int NUMBER_OF_THREADS = 2;
static int MATRIX_SIZE = 5;
static int RANNGE_FROM = -100;
static int RANGE_TO = 100;
static bool collect_communication_time = false;
using type = int;
//using MPI_type = MPI_INT

// using Matrix = std::vector<std::vector<int>>;

void T(type *m){
    //
    // transponowanie macierzy (no odwracamy się ) otrzymujemy m^T
    //
    type mT[MATRIX_SIZE * MATRIX_SIZE] = {0};
    for (int i = 0; i < MATRIX_SIZE; i++)
        for (int j = 0; j < MATRIX_SIZE; j++){
            mT[i * MATRIX_SIZE + j] = m[j * MATRIX_SIZE + i];
        }
    for (int i = 0; i < MATRIX_SIZE*MATRIX_SIZE; i++)
        m[i] = mT[i];
}

void generate_data(type *m,const bool spawn_unitary = false){
    //
    // generujemy dane i wkładamy je do std::cout<<m
    //
    // dane są z rozkładu jednostajnego z zakresu <range_from,range_to)
    //
    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>  distr(RANNGE_FROM, RANGE_TO);
    //type sum;
    if(spawn_unitary == false){
        for(int i = 0; i<MATRIX_SIZE*MATRIX_SIZE; i++){
            m[i] = distr(generator);
        }
    }else{
        for(int i = 0; i<MATRIX_SIZE*MATRIX_SIZE; i++){
            m[i] = 0;
        }
        for(int i = 0; i<MATRIX_SIZE; i++){
            m[i + MATRIX_SIZE * i] = 1;
        }
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

void reed_matsize(int argc,char** argv){
  //std::cout<<"reed_states";
  if(argc > 1){
    MATRIX_SIZE = std::stoi(argv[argc -1]);
  }//std::cout<<"reed_states";
  if(argc > 4 || argc == 3){//To daje możliwość odpalenia Size,comunication/c,rangeFrom,rangeTo lub Size,communication/c
    char flag = argv[argc - 2][0]; // Read the first character of the argument
    collect_communication_time = (flag == 'c'); // Compare it with 'c'
  }
  if(argc > 4){
    RANGE_TO = std::stoi(argv[argc -3]);
  }
  if(argc > 4){ // To daje możliwość odpalenia Size,comunication,rangeFrom,rangeTo
    RANNGE_FROM = std::stoi(argv[argc -4]);
  }
  if(argc == 4){ // To daje możliwość odpalenia Size,rangeFrom,rangeTo
    RANGE_TO = std::stoi(argv[argc -2]);
  }
  if(argc == 4){
    RANNGE_FROM = std::stoi(argv[argc -3]);
  }
  //std::cout<<"reed_states 2";
}


int main(int argc,char** argv) {
    //std::cout<<"reed_states 0";
    MPI_Init(&argc,&argv);// inicjalizacji MPI
    //std::cout<<"reed_states -1";
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    //std::cout<<"reed_states -2";
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //std::cout<<"reed_states -3";
    reed_matsize(argc,argv);

    type *a;
    a = (type*) malloc(sizeof(type)* MATRIX_SIZE * MATRIX_SIZE);
    type *b;
    b = (type*) malloc(sizeof(type)* MATRIX_SIZE * MATRIX_SIZE);
    type *c;
    c = (type*) malloc(sizeof(type)* MATRIX_SIZE * MATRIX_SIZE);
    auto start_time = std::chrono::steady_clock::now();
    auto end_time = std::chrono::steady_clock::now();
    double comm_time = 0.0;
    std::chrono::duration<double,std::milli> exec_time;
    
    const type range_start = 0;
    const type range_end = 10;
    
    if(world_rank == 0){
      generate_data(a);
      generate_data(b,false);
        if(MANUAL_TEST){
            std::cout<<"matrix A"<<std::endl;
            showMatrix(a);
            std::cout<<"matrix B"<<std::endl;
            showMatrix(b);
        }
      start_time = std::chrono::steady_clock::now();
      comm_time -= MPI_Wtime();
    }
    
    MPI_Bcast(a,MATRIX_SIZE * MATRIX_SIZE,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(b,MATRIX_SIZE * MATRIX_SIZE,MPI_INT,0,MPI_COMM_WORLD);
    
    if(world_rank == 0){
      comm_time += MPI_Wtime();
    }

    int start_end_operation[2] = {0,0};
    
    const int all_operations = (MATRIX_SIZE * MATRIX_SIZE);
    const int op_for_one = all_operations/world_size;
    const int rest_op = all_operations%world_size;
    //std::cout<<"all_op:" << all_operations << ", op_for_one:"<<op_for_one <<", rest_op:"<<rest_op <<std::endl;
    if(world_rank == 0){
        for(int thread_id = world_size-1; thread_id >= 0; thread_id--){
            start_end_operation[0] = op_for_one * thread_id;
            start_end_operation[1] = op_for_one * (thread_id+1) + rest_op * ( thread_id == world_size-1 );
            
            if(thread_id != 0){
                comm_time -= MPI_Wtime();
                MPI_Send(&start_end_operation, 2, MPI_INT,thread_id , 0, MPI_COMM_WORLD);
                comm_time += MPI_Wtime();
            }
        }
    }else{
        MPI_Recv(&start_end_operation,2,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
    }
    
    //std::cout<<"worker_rank:" << world_rank << ", start_op:"<< start_end_operation[0] << ", end_op:"<< start_end_operation[1] <<std::endl;
    matmul_MPI(c, a, b, start_end_operation[0],start_end_operation[1]);
    if(world_rank == 0){
        //showMatrix(c);
        type *RecvMatrix;
        RecvMatrix = (type*) malloc(sizeof(type)* (op_for_one + rest_op));
        //std::cout<<"worker 0 allocated memory RecvMatrix = "<<RecvMatrix<<std::endl;
        for(int thread_id = world_size-1; thread_id > 0; thread_id--){
            //std::cout<<"worker 0 strat loop:"<<thread_id<<std::endl;
            comm_time -= MPI_Wtime();
            MPI_Recv(RecvMatrix,(op_for_one + rest_op* ( thread_id == (world_size-1) )),MPI_INT,thread_id,thread_id,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            comm_time += MPI_Wtime();
            //std::cout<<"worker 0 recived data"<<std::endl;
            for (int i = 0; i< (op_for_one + rest_op* ( thread_id == (world_size-1) )); i++){
            //std::cout<<"worker 0|  FROM c["<<thread_id * op_for_one + i<<"] = rcv["<<i<<"]"<<std::endl;
            //std::cout<<"worker 0|  c["<<thread_id * op_for_one + i<<"] = "<<c[thread_id * op_for_one + i]<<std::endl;
            //std::cout<<"worker 0|  rcv["<<i<<"] = "<<RecvMatrix[i]<<std::endl;
              c[thread_id * op_for_one + i] = RecvMatrix[i];
            }
            //std::cout<<"worker 0 copied data"<<std::endl;
            //showMatrix(c);
        }
        free(RecvMatrix);
        end_time = std::chrono::steady_clock::now();
        exec_time = end_time-start_time;
        char input;
        if(MANUAL_TEST){
            //showMatrix(c);
            std::cout<<"Solution achived in: "<< exec_time.count() << "ms" <<std::endl;
            std::cout<<"Communication time: "<< comm_time*1000.0 << "ms"<<std::endl;
        }else{
            if(collect_communication_time){
                std::cout<< comm_time*1000.0;
            }else{
                std::cout<< exec_time.count();
            }
        }
    }else{
        //showMatrix(c);
        const int num_of_op = start_end_operation[1] - start_end_operation[0];
        type *sender;
        sender = (type*) malloc(sizeof(type)* num_of_op);
        for(int i = 0;i < num_of_op; i++){
          sender[i] = c[i + start_end_operation[0]];
          //std::cout<<"worker_rank:" << world_rank <<" calculated position"<< i + start_end_operation[0] << ", as value equal:"<< sender[i] <<std::endl;
        }
        //std::cout<<"worker_rank:" << world_rank << " start_sending"<<std::endl;
        MPI_Send(sender,num_of_op,MPI_INT, 0,world_rank, MPI_COMM_WORLD);
        free(sender);
    }
    free(a);
    free(b);
    free(c);
    MPI_Finalize();

    return 0;
}
