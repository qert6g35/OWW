#include <iostream>
#include <fstream>
#include <chrono>
#include <thread>
#include <vector>

static const int NUMBER_OF_THREADS = 4;
static const int MATRIX_SIZE = 5;
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
        col_b = (elem_id%MATRIX_SIZE); // zakładam że transponowano b i ze indeksowanie teraz jest z gory na dol

        start_a = (row_a*MATRIX_SIZE);
        // finish_a = start_a + MATRIX_SIZE;

        start_b = (col_b*MATRIX_SIZE); 
        // finish_b = start_b + MATRIX_SIZE;
        elem_of_c = 0;
        for(int i = 0; i < MATRIX_SIZE; i++){
            elem_of_a = start_a + i;
            elem_of_b = start_b + i;
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

int main() {
    type a[MATRIX_SIZE * MATRIX_SIZE];
    type b[MATRIX_SIZE * MATRIX_SIZE];
    type c[MATRIX_SIZE * MATRIX_SIZE] = {0};

    for (int i = 0; i < MATRIX_SIZE * MATRIX_SIZE; ++i) {
        a[i] = i;
        b[i] = i;
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