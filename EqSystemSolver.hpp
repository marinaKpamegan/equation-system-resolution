#ifndef EQSYSTEMSOLVERHEADERDEF
#define EQSYSTEMSOLVERHEADERDEF

#include <iostream>
#include <vector> 

// using matrix = std::vector<std::vector<double>>;


class EqSystemSolver
{
private:
    int size;
    int nCols;
    int nRows;
    std::vector<std::vector<double>> A;
    std::vector<double> b;

public:
    EqSystemSolver(int size, std::vector<std::vector<double>> A, std::vector<double> b);
    friend double calculateDeterminant(int size, std::vector<std::vector<double>> A);
    friend double getCofactor(int size, std::vector<std::vector<double>> A);
    double Determinant() const;
    std::vector<std::vector<double>> GetA() const;
    std::vector<std::vector<double>> GetB() const;
    friend std::vector<std::vector<double>> getSubMatrix(std::vector<std::vector<double>> A, int i, int j);
    void SetA(std::vector<std::vector<double>> newA);
    void SetB(std::vector<double> newB);
    std::vector<std::vector<double>> Transpose() const;
    friend std::vector<double> mult(int size, std::vector<std::vector<double>> A, std::vector<double> b);
    friend std::vector<std::vector<double>> inverseMatrix(int size, std::vector<std::vector<double>> A);
    std::vector<std::vector<double>> AugmentedMatrix();
    std::vector<std::vector<double>> Inverse();
    std::vector<double> SolverGauss();
    std::vector<double> SolverInverseMatrix();
};

#endif