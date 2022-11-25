//#include "functions.cpp"
#include "EqSystemSolver.hpp"
#include <iostream>
#include <cassert>
using namespace std;

int main (){
    double det;
    int n;
    int choix = 1;
    
    // det = e.Determinant();
    std::cout<<"Nombres d'inconnues "<<std::endl;
    std::cin>>n;
    // det = e.Determinant();
    // std::cout<<det;
    // std::cout<<calculateDeterminant(2, A);
    assert(n>0);

    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    std::vector<double> b(n);
    vector<double> x;
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

    std::cout<<"Choisir la méthode de résolution."<<std::endl;
    std::cout<<"1. Pivot de Gauss"<<std::endl;
    std::cout<<"2. Matrice inverse"<<std::endl;
    std::cin>> choix;

    if(choix == 1){
        x = e.SolverGauss();
        for (int i = 0; i < n; i++)
        {
            std::cout << x[i] << " ";
        }
        std::cout << "\n";
    }else if(choix== 2){
        x = e.SolverInverseMatrix();
        for (int i = 0; i < n; i++)
        {
            std::cout << x[i] << " ";
        }
        std::cout << "\n";
    }else{
        std::cout<<"Mauvais choix. Réessayez"<<std::endl;
        std::cin>> choix;
    }
   
    
    
}

