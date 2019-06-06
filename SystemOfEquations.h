#ifndef PR2_SYSTEMOFEQUATIONS_H
#define PR2_SYSTEMOFEQUATIONS_H

#define  VEC_RES 0
#define  VEC_X   1
#define  VEC_B   2
#define  MIN_RES 1.0 / 1000000000

#endif //PR2_SYSTEMOFEQUATIONS_H



class SystemOfEquations {
private:
    double *vec_res;
    double **mat_A;
    double *vec_x;
    double *vec_x1;
    double *vec_b;
    int N;

public:
    SystemOfEquations();
    void print_matrix_A();
    void print_vector(int what_to_print);
    void save_x_vector();

    double sum_for_iter(int j, int to, int i, double *vector);
    double count_row(int row);
    double residium_norm();

    //methods
    void jacobi();
    void gauss_seidel();
    void factorizationLU();

    ~SystemOfEquations();
};
