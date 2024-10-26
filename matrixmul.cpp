#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <vector>

static const int NUMBER_OF_THREADS = 4;
static const int MATRIX_SIZE = 2;
using type = int;

// using Matrix = std::vector<std::vector<int>>;

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


// // Funkcja mnożąca fragment macierzy
void matmul_thread(type *c, const type *a, const type *b, int thread_id) {

}

// Funkcja inicjująca wątki
void matmul_threaded(type *c, const type *a, const type *b) {
    std::thread threads[NUMBER_OF_THREADS];
    for (int i = 0; i < NUMBER_OF_THREADS; ++i) {
        threads[i] = std::thread(matmul_thread, c, a, b, i);
    }

    for (int i = 0; i < NUMBER_OF_THREADS; ++i) {
        threads[i].join();
    }
    
}

int main() {
    type a[MATRIX_SIZE * MATRIX_SIZE];
    type b[MATRIX_SIZE * MATRIX_SIZE];
    type c[MATRIX_SIZE * MATRIX_SIZE] = {0};

    for (int i = 0; i < MATRIX_SIZE * MATRIX_SIZE; ++i) {
        a[i] = rand() % 10;
        b[i] = rand() % 10;
    }

    std::cout << "Matrix A:" << std::endl;
    showMatrix(a);
    std::cout << "Matrix B:" << std::endl;
    showMatrix(b);

    matmul_threaded(c, a, b);
    std::cout << "Result Matrix C:" << std::endl;
    showMatrix(c);

    return 0;
}