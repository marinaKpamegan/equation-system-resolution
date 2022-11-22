#include "functions.cpp"
#include "EqSystemSolver.hpp"

int main (){
    vector<vector<double>> A = { { 1, 3, 2 },
                       {  3, -2, 3 },
                       { -1, 3, 5 } };
    // double b[3] = {1, 2, 3} ;

    // EqSystemSolver e = EqSystemSolver(3, A, b);
    //printMatrix(3, A);
    std::cout<<calculateDeterminant(3, A);
}


