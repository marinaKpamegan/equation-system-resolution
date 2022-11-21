#include "EqSystemSolver.hpp"
#include <cmath>
#include <cassert>
#include "functions.cpp"


EqSystemSolver::EqSystemSolver(int size, double** A, double* b){
    
    this->A = A;
    this->b = b;
    this->size = size;
}

double EqSystemSolver::CalculateDeterminant(int size, double** A){
    double det = 0.0;
    double** subMatrix = allocMemory(size); // submatrix will contains a part of the matrix
    int sign = 1; // sign of covariants (1 for + and -1 for -)

    assert(subMatrix != NULL); // checking if memory for temp matrix is allocated 
    
    if(size == 2){

        return ((A[0][0] * A[1][1]) - (A[1][0] * A[0][1]));

    }else {

        for(int k=0; k<size; k++){

            int subMatrixI = 1; // submatrix index
            for(int i=1; i<size; i++){

                int subMatrixJ = 0; // submatrix index
                for(int j=0; j<size; j++){
                    if(j==k)
                        continue;
                    subMatrix[subMatrixI][subMatrixJ] = A[i][j]; 
                    subMatrixJ++;                  
                }

                subMatrixI++;
            }
            det = det + (pow(-1, k) * A[0][k] * CalculateDeterminant(size - 1, subMatrix));
        }
    }

    freeMemory(size, subMatrix); // free allocated memory

    return det;
}




