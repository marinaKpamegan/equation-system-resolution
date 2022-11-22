#include "EqSystemSolver.hpp"
#include <cmath>
#include <cassert>
#include <vector> 

using namespace std;
EqSystemSolver::EqSystemSolver(int size, vector<vector<double>> A, vector<double> b){
    
    this->A = A;
    this->b = b;
    this->size = size;
}

