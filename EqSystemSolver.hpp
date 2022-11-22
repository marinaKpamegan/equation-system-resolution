#ifndef EQSYSTEMSOLVERHEADERDEF
#define EQSYSTEMSOLVERHEADERDEF

#include <iostream>
#include <vector> 

using namespace std;
class EqSystemSolver
{
private:
    int size;
    vector<vector<double>> A;
    vector<double> x;
    vector<double> b;

public:
    EqSystemSolver(int size, vector<vector<double>> &A, vector<double> &b);
    friend double calculateDeterminant(int size, vector<vector<double>> A);
    EqSystemSolver Transpose(EqSystemSolver& A) const;
    EqSystemSolver Triangularize(EqSystemSolver& A, double* b) const;
    EqSystemSolver Mult(EqSystemSolver& A, double b[]) const;
    friend std::ostream& operator<<(std::ostream& output,
    const EqSystemSolver& A);
    // friend vector<vector<double>> allocMemory();
    // friend void freeMemory();
};

#endif