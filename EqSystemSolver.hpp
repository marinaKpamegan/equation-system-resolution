#ifndef EQSYSTEMSOLVERHEADERDEF
#define EQSYSTEMSOLVERHEADERDEF

#include <iostream>
#include <vector> 

using matrix = vector<vector<double>>;

using namespace std;
class EqSystemSolver
{
private:
    int size;
    vector<vector<double>> A;
    vector<double> x;
    vector<double> b;

public:
    EqSystemSolver(int size, vector<vector<double>> A, vector<double> b);
    friend double calculateDeterminant(int size, vector<vector<double>> A);
    EqSystemSolver Transpose() const;
    EqSystemSolver Triangularize() const;
    bool isUpperTriangular();
    bool isLowerTriangular();
    friend vector<double>  mult(int size, vector<vector<double>> A, vector<double> b);
    friend std::ostream& operator<<(std::ostream& output,
    const EqSystemSolver& A);
    vector<vector<double>> Concatenate();
    vector<vector<double>> Reverse() const;
    // vector<double> Solver() const;
    // friend vector<vector<double>> allocMemory();
    // friend void freeMemory();
};

#endif