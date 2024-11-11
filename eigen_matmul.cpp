#include <Eigen/Dense>
#include <iostream>
#include <cstdlib>  // For std::atoi
#include <chrono>   // For timing
#include <random>

static int MATRIX_SIZE = 5;
static int RANNGE_FROM = -100;
static int RANGE_TO = 100;


Eigen::MatrixXd generate_data(){
    //
    // generujemy dane i wkładamy je do std::cout<<m
    //
    // dane są z rozkładu jednostajnego z zakresu <range_from,range_to)
    //

    Eigen::MatrixXd matrix(MATRIX_SIZE,MATRIX_SIZE);

    std::random_device                  rand_dev;
    std::mt19937                        generator(rand_dev());
    std::uniform_int_distribution<int>  distr(RANNGE_FROM, RANGE_TO);
    //type sum;
    for(int j = 0; j<MATRIX_SIZE;j++){
        for(int i = 0; i<MATRIX_SIZE; i++){
            matrix(i,j) = distr(generator);
        }
    }
    return matrix;
    //std::cout << "sum:" << sum << "   we sum n values, n:"<< MATRIX_SIZE*MATRIX_SIZE <<std::endl;
}

void reed_args(int argc,char** argv){
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


int main(int argc, char* argv[]) {
    std::chrono::duration<double,std::milli> exec_time;
    reed_args(argc,argv);

    // Define dynamic matrices of the specified size
    Eigen::MatrixXd matrixA = generate_data();
    Eigen::MatrixXd matrixB = generate_data();

    // Measure the time taken for multiplication
    auto start_time = std::chrono::steady_clock::now();
    
    // Perform matrix multiplication
    Eigen::MatrixXd result = matrixA * matrixB;
    
    auto end_time = std::chrono::steady_clock::now();
    exec_time = end_time-start_time;

    std::cout<< exec_time.count();

    return 0;
}