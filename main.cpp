#include <iostream>
#include <sstream>
#include <string>
#include <iostream>
#include "SystemOfEquations.h"

#define JACOBI           1
#define GAUSS_SEIDEL     2
#define LU_FACTORIZATION 3



int main(int argc, char *argv[]) {
    int method;
    printf("This program will solve your system of linear equations (place it in Matrix.txt) for vector b (place it in B.txt).\n");
    printf("Now please choose one of the methods:\n1 - Jacobi\n2 - Gauss-Seidel\n3 - LU factorization\n");

    scanf("%d", &method);

    if (method != 1 && method != 2 && method != 3) {
        printf("Wrong number, try again!\n");
        return 1;
    }

    SystemOfEquations *systemToSolve = new SystemOfEquations();
    printf("Solving...\n");

    switch (method) {
        case JACOBI:
            systemToSolve->jacobi();
            break;
        case GAUSS_SEIDEL:
            systemToSolve->gauss_seidel();
            break;
        case LU_FACTORIZATION:
            systemToSolve->factorizationLU();
            break;
    }

    systemToSolve->print_matrix_A();
    systemToSolve->print_vector(VEC_B);

    systemToSolve->save_x_vector();
    printf("Solved! Check your solution in X.txt\n");
    delete systemToSolve;

    return 0;
}
