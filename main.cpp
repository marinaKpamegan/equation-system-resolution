//#include "functions.cpp"
#include "EqSystemSolver.hpp"
#include <iostream>
#include <cmath>
using namespace std;

int main (){
    vector<vector<double>> A = { {  1, 1, 1 },
                                {  1, 2, 3 },
                                {  -2, 3, -1 }};
    vector<double> b = {3, 6, 0} ;
    vector<double> x;
    double det;
    int n;
    
    // det = e.Determinant();
    std::cout<<"Nombres d'inconnues "<<std::endl;
    std::cin>>n;
    // det = e.Determinant();
    // std::cout<<det;
    // std::cout<<calculateDeterminant(2, A);
    std::cout<<"Entrez les valeurs de la matrice A de l'équation:"<<std::endl;
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            std::cout<<"A["<<i<<"]["<<j<<"]=";
            std::cin>> A[i][j];
        }
    }
    
    std::cout<<"Entrez les valeurs du vecteur b de l'équation:"<<std::endl;
    for(int i=0; i<n; i++){
        std::cout<<"b["<<i<<"]=";
        std::cin>> b[i];
    }
    EqSystemSolver e(n, A, b);
    x = e.SolverGauss();
    std::cout<<x[0]<<" "<<x[1]<<" "<<x[2];
}

