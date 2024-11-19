#include <iostream>
#include <fstream>
#include <random>
#include <chrono>
#include <thread>
#include <vector>
#include <cstdlib>  // For std::atoi

static int MATRIX_SIZE = 5;
static int RANNGE_FROM = -100;
static int RANGE_TO = 100;

using type = int;

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

void read_args(int argc,char** argv){
  if(argc > 1){
    MATRIX_SIZE = std::stoi(argv[argc -1]);
  }
  if(argc > 2){
    RANGE_TO = std::stoi(argv[argc -2]);
  }
  if(argc > 3){
    RANNGE_FROM = std::stoi(argv[argc -3]);
  }
}

void matmul_single(type *c, const type *a, const type *b) {
    for (int i = 0; i < MATRIX_SIZE; i++)
        for (int j = 0; j < MATRIX_SIZE; j++)
            for (int k = 0; k < MATRIX_SIZE; k++)
                c[i * MATRIX_SIZE + j] += a[i * MATRIX_SIZE + k] * b[k * MATRIX_SIZE + j];
}


int main(int argc, char* argv[]) {
    std::chrono::duration<double,std::milli> exec_time;
    read_args(argc,argv);
    type *a;
    a = (type*) malloc(sizeof(type)* MATRIX_SIZE * MATRIX_SIZE);
    type *b;
    b = (type*) malloc(sizeof(type)* MATRIX_SIZE * MATRIX_SIZE);
    type *c;
    c = (type*) malloc(sizeof(type)* MATRIX_SIZE * MATRIX_SIZE);
    // Define dynamic matrices of the specified size
    generate_data(a);
    generate_data(b,false);

    // Measure the time taken for multiplication
    auto start_time = std::chrono::steady_clock::now();

    // Perform matrix multiplication
    matmul_single(c,a,b);
    
    auto end_time = std::chrono::steady_clock::now();
    exec_time = end_time-start_time;

    std::cout<< exec_time.count();

    return 0;
}