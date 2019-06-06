#include "SystemOfEquations.h"
#include <cmath>
#include <cstdio>
#include <chrono>
#include <fstream>



SystemOfEquations::SystemOfEquations() {
    std::fstream matrix_file("Matrix.txt");
    std::fstream b_vector_file("B.txt");

    matrix_file >> this->N;

    //initialization of vector x and residium
    this->vec_x = new double[N];
    this->vec_x1 = new double[N];
    for (int i = 0; i < N; i++) {
        this->vec_x[i] = 0;
        this->vec_x1[i] = 0;
    }

    this->vec_res = new double[N];

    //array initialization
    this->mat_A = new double*[N];
    for (int i = 0; i < N; i++)
        this->mat_A[i] = new double[N];

    //build a custom matrix
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            matrix_file >> this->mat_A[i][j];

    //initialization of vector b
    this->vec_b = new double[N];
    for (int i = 0; i < N; i++) {
        b_vector_file >> this->vec_b[i];
    }

    matrix_file.close();
    b_vector_file.close();
}


void SystemOfEquations::print_matrix_A() {
    printf("\n========= MATRIX A ==========\n");
    for (int i = 0; i < this->N; i++) {
        for (int j = 0; j < this->N; j++)
            printf("%2.0lf ", this->mat_A[i][j]);
        printf("\n");
    }
    printf(  "=============================\n");
}


void SystemOfEquations::print_vector(int what_to_print) {
    double *vector_to_print;

    switch (what_to_print) {
        case VEC_RES:
            printf("\nVector Residium'\n");
            vector_to_print = this->vec_res;
            break;
        case VEC_X:
            printf("\nVector X'\n");
            vector_to_print = this->vec_x;
            break;
        case VEC_B:
            printf("\nVector B'\n");
            vector_to_print = this->vec_b;
            break;
        default:
            printf("\nUndifined vector.\n");
            return;
    }

    printf("< ");
    for (int i = 0; i < this->N; i++)
        printf("%2.2lf ", vector_to_print[i]);
    printf(">\n");
}


void SystemOfEquations::save_x_vector() {
    std::ofstream vector_x_file("X.txt");

    vector_x_file << "Vector X:\n";

    for (int i = 0; i < this->N; i++)
        vector_x_file << this->vec_x[i] << " ";

    vector_x_file.close();
}


double SystemOfEquations::sum_for_iter(int from, int to, int i, double *vector) {
    double result = 0;
    for (; from < to; from++) {
        if (from != i)
            result += this->mat_A[i][from] * vector[from];
    }
    return result;
}


double SystemOfEquations::count_row(int row) {
    double result = 0;
    for (int i = 0; i < this->N; i++)
        result += this->mat_A[row][i] * this->vec_x[i];
    return result;
}


double SystemOfEquations::residium_norm() {
    double result = 0;
    for (int i = 0; i < this->N; i++) {
        this->vec_res[i] = count_row(i) - this->vec_b[i];
    }
    for (int i = 0; i < this->N; i++)
        result += this->vec_res[i] * this->vec_res[i];
    return sqrt(result);
}


void SystemOfEquations::jacobi() {
    int iter;
    std::chrono::high_resolution_clock::time_point start_point = std::chrono::high_resolution_clock::now();

    for (iter = 0; residium_norm() > MIN_RES; iter++) {
        for (int i = 0; i < this->N; i++)
            this->vec_x1[i] = (this->vec_b[i] - sum_for_iter(0, this->N, i, this->vec_x)) / this->mat_A[i][i];

        for (int i = 0; i < this->N; i++)
            this->vec_x[i] = this->vec_x1[i];

        if (iter > 100) {
            printf("Jacobi method: more than 100 iterations. Aborted with norm(residium) = %.2e.\n", residium_norm());
            return;
        }
    }

    std::chrono::high_resolution_clock::time_point end_point = std::chrono::high_resolution_clock::now();
    auto duration = end_point - start_point;

    printf("\n========================\nJacobi iterations = %d\n", iter);
    printf("time = %lf sec\n========================\n", duration / 1000000000.0);
}


void SystemOfEquations::gauss_seidel() {
    int iter;
    std::chrono::high_resolution_clock::time_point start_point = std::chrono::high_resolution_clock::now();

    for (iter = 0; residium_norm() > MIN_RES; iter++) {
        for (int i = 0; i < this->N; i++)
            this->vec_x1[i] = (this->vec_b[i] - sum_for_iter(0, i, i, this->vec_x1) -
                               sum_for_iter(i + 1, this->N, i, this->vec_x)) / this->mat_A[i][i];

        for (int i = 0; i < this->N; i++)
            this->vec_x[i] = this->vec_x1[i];

        if (iter > 100) {
            printf("Gauss-Seidel method: more than 100 iterations. Aborted with norm(residium) = %.2e.\n", residium_norm());
            return;
        }
    }

    std::chrono::high_resolution_clock::time_point end_point = std::chrono::high_resolution_clock::now();
    auto duration = end_point - start_point;

    printf("\n==============================\nGauss-Seidel iterations = %d\n", iter);
    printf("time = %lf sec\n==============================\n", duration / 1000000000.0);
}


void SystemOfEquations::factorizationLU() {
    double tmp;
    double **L;
    double **U;
    double  *y;

    //initialization
    L = new double*[this->N];
    U = new double*[this->N];
    y = new double[this->N];

    for (int i = 0; i < this->N; i++) {
        L[i] = new double[this->N];
        U[i] = new double[this->N];
        y[i] = 0;
    }

    std::chrono::high_resolution_clock::time_point start_point = std::chrono::high_resolution_clock::now();
    //making U and L
    for (int i = 0; i < this->N; i++)
        for (int j = 0; j < this->N; j++) {
            L[i][j] = (i == j) ? 1 : 0;
            U[i][j] = this->mat_A[i][j];
        }

    for (int k = 0; k < this->N - 1; k++)
        for (int j = k + 1; j < this->N; j++) {
            L[j][k] = U[j][k] / U[k][k];
            for (int i = k; i < this->N; i++)
                U[j][i] = U[j][i] - L[j][k] * U[k][i];
        }

    //Ly = b
    for (int i = 0; i < this->N; i++) {
        tmp = this->vec_b[i];
        for (int j = 0; j < i; j++)
            if (j != i)
                tmp -= L[i][j] * y[j];
        y[i] = tmp / L[i][i];
    }

    //Ux = y
    for (int i = this->N - 1; i >= 0; i--) {
        tmp = y[i];
        for (int j = i; j < this->N; j++)
            if (j != i)
                tmp -= U[i][j] * this->vec_x[j];
        this->vec_x[i] = tmp / U[i][i];
    }

    std::chrono::high_resolution_clock::time_point end_point = std::chrono::high_resolution_clock::now();
    auto duration = end_point - start_point;
    printf("\n========================\nLU factorization\n");
    printf("time = %lf sec\n========================\n", duration / 1000000000.0);


    //memory free
    for (int i = 0; i < this->N; i++) {
        delete L[i];
        delete U[i];
    }
    delete L;
    delete U;
    delete y;
}


SystemOfEquations::~SystemOfEquations() {
    delete this->vec_x;
    delete this->vec_x1;
    delete this->vec_res;
    delete this->vec_b;
    for (int i = 0; i < this->N; i++)
        delete this->mat_A[i];
    delete this->mat_A;
}
