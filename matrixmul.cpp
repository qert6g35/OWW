#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

static const int NUMBER_OF_THREADS = 2;
static const int MATRIX_SIZE = 2;
typedef float type;

void matmul(const type *a, const type *b, type *c) {
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

void runThreds(){
    type a[2][2] = {
                    { 1, 2 },
                    { 3, 4 }
                    };
    type b[2][2] = {
                { 1, 0 },
                { 0, 1 }
                };
    type c[2][2];

    std::thread Threads[NUMBER_OF_THREADS];

    for(int i = 0; i<NUMBER_OF_THREADS; i++){
        Threads[i] = std::thread(matmul,*a,*b,*c);
    }

    for(int i = 0; i<NUMBER_OF_THREADS; i++){
        Threads[i].join();
    }

    showMatrix(*c);

}

int main(){
    runThreds();
    return 0;
}

