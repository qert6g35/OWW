#include <iostream>
#include <fstream>

void matmul(const float *a, const float *b, float *c, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                c[i * n + j] += a[i * n + k] * b[k * n + j];
}

int main(){
    float a[2][2] = {
                    { 1, 2 },
                    { 3, 4 }
                    };
    float b[2][2] = {
                { 1, 2 },
                { 3, 4 }
                };
    float c[2][2];
    matmul(*a, *b, *c, 2);
    std::cout << c[1][1] << std::endl;
    return 0;
}