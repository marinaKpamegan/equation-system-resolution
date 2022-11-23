#include "EqSystemSolver.hpp"
#include <cmath>
#include <cassert>
#include <vector> 
#include "functions.cpp"

using namespace std;
EqSystemSolver::EqSystemSolver(int size, vector<vector<double>> A, vector<double> b){
    
    this->A = A;
    this->b = b;
    this->size = size;
}

vector<vector<double>> EqSystemSolver::Concatenate(){

}

vector<vector<double>> EqSystemSolver::Reverse() const{
    double det;
    vector<vector<double>> matrixReverse(size, vector<double>(size));
    det = calculateDeterminant(size, A);
    assert(det !=0);

    for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
            matrixReverse[i][j] = 1/det * A[i][j];
        }
    }
    return matrixReverse;  // A moins 1 

}




