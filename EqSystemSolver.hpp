#ifndef EQSYSTEMSOLVERHEADERDEF
#define EQSYSTEMSOLVERHEADERDEF

#include <iostream>

class EqSystemSolver
{
private:
    int size;
    double** A;
    double* x;
    double* b;

public:
    EqSystemSolver(int size, double* A[], double* b);
    double CalculateDeterminant(int size, double** A);
    EqSystemSolver Transpose(EqSystemSolver& A) const;
    EqSystemSolver Triangularize(EqSystemSolver& A, double* b) const;
    EqSystemSolver Mult(EqSystemSolver& A, double b[]) const;
    friend std::ostream& operator<<(std::ostream& output,
    const EqSystemSolver& A);
    friend double** allocMemory();
    friend void freeMemory();
};

#endif