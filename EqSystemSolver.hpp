#ifndef EQSYSTEMSOLVERHEADERDEF
#define EQSYSTEMSOLVERHEADERDEF

#include <iostream>
#include <vector> 

// using matrix = std::vector<std::vector<double>>;


class EqSystemSolver
{
private:
    int size;
    std::vector<std::vector<double>> A;
    // std::vector<double> x;
    std::vector<double> b;

public:
    EqSystemSolver(int size, std::vector<std::vector<double>> A, std::vector<double> b);
    friend double calculateDeterminant(int size, std::vector<std::vector<double>> A);
    double Determinant() const;
    EqSystemSolver Transpose() const;
    friend std::vector<double>  mult(int size, std::vector<std::vector<double>> A, std::vector<double> b);
    friend std::vector<std::vector<double>> reverseMatrix(int size, std::vector<std::vector<double>> A);
    std::vector<std::vector<double>> AugmentedMatrix();
    std::vector<std::vector<double>> Reverse();

    std::vector<double> SolverGauss();
};

#endif