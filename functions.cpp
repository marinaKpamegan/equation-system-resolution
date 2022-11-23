#include <vector> 
#include <cmath>
#include <iostream>

using namespace std;


double** allocMemory(int size){
    double** matrix;

    matrix = new double*[size];
    for (int i=0; i<size; i++){
        matrix[i] = new double[size];
    }
    return matrix;
}

void freeMemory(int size, double** matrix){
    
    for (int i=0; i<size; i++){
        delete[] matrix[i];
    }
    delete[] matrix;
}


/* void printMatrix(int size, vector<vector<double>> A){
    for(int i=0; i<size; i++){
        
        for(int j=0; j<size; j++){
            std::cout<<A[i][j] <<" ";
        }
        std::cout<<"\n";
    }

    std::cout<<"\n";
}*/



double calculateDeterminant(int size, vector<vector<double>> A){
    double det = 0.0;
    std::vector<vector<double>> subMatrix(size, vector<double>(size)); // submatrix will contains a part of the matrix
    int sign = -1; // sign of covariants (1 for + and -1 for -)
    int subMatrix_i , subMatrix_j;
    
    if(size == 2){        
        return ((A[0][0] * A[1][1]) - (A[1][0] * A[0][1]));
    }else {

        for(int k=0; k<size; k++){
            subMatrix_i = 0;
            for(int i=1; i<size; i++){
                subMatrix_j = 0;
                for(int j=0; j<size; j++){

                    // std::cout<< "i="<<i<< " j="<<j <<" k="<<k<<"\n";
                    if(j!=k){
                        subMatrix[subMatrix_i][subMatrix_j] = A[i][j];            
                        subMatrix_j++;                        
                    }else{
                        sign = pow(-1, i+j);                        
                    }
                    
                }
                subMatrix_i++;
            }
            
            det = det + (sign * A[0][k] * calculateDeterminant(size - 1, subMatrix));
        }

    }

    // freeMemory(size, subMatrix); // free allocated memory

    return det;
}

vector<double>  mult(int size, vector<vector<double>> A, vector<double> b){
    vector<double> Ab;

    for(int i=0; i<size; i++){
        int sum = 0;
        for(int j=0; j<size; j++){
            sum += A[i][j] * b[j];
        }
        Ab[i] = sum;
    }

    return Ab;
}



