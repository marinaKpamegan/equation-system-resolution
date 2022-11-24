//#include "functions.cpp"
#include "EqSystemSolver.hpp"
#include <iostream>
#include <cmath>
using namespace std;

int main (){
    vector<vector<double>> A = { { 0, 1, 1},
                                {  1, 2, 3 },
                                {  -2, 3, -1 } };
    vector<double> b = {3, 6, 0} ;
    vector<double> x(3);
    double det;
    EqSystemSolver e(3, A, b);
    // det = e.Determinant();
    // det = e.Determinant();
    // std::cout<<det;
    // std::cout<<calculateDeterminant(2, A);
    x = e.SolverGauss();
    std::cout<<x[0]<<" "<<x[1]<<" "<<x[2];
}

