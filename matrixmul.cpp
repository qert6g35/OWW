#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>

static const int NUMBER_OF_THREADS = 2;
static const int n = 2;
typedef float type;

void matmul(const type *a, const type *b, type *c) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                c[i * n + j] += a[i * n + k] * b[k * n + j];
}

void showMatrix(const type *matrix){
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            std::cout<<matrix[i * n + j]<<" ";
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

