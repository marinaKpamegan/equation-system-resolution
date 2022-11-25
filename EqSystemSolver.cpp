

#include "EqSystemSolver.hpp"
#include <cmath>
#include <cassert>
#include <vector>

EqSystemSolver::EqSystemSolver(int n, std::vector<std::vector<double>> matrix, std::vector<double> B)
{
    A = matrix;
    b = B;
    size = n;
}

double calculateDeterminant(int size, std::vector<std::vector<double>> A)
{
    double det = 0.0;
    std::vector<std::vector<double>> subMatrix(size, std::vector<double>(size)); // submatrix will contains a part of the matrix
    int sign = -1;                                                               // sign of covariants (1 for + and -1 for -)
    int subMatrix_i, subMatrix_j;

    if (size == 2)
    {
        return ((A[0][0] * A[1][1]) - (A[1][0] * A[0][1]));
    }
    else
    {

        for (int k = 0; k < size; k++)
        {
            subMatrix_i = 0;
            for (int i = 1; i < size; i++)
            {
                subMatrix_j = 0;
                for (int j = 0; j < size; j++)
                {

                    // std::cout<< "i="<<i<< " j="<<j <<" k="<<k<<"\n";
                    if (j != k)
                    {
                        subMatrix[subMatrix_i][subMatrix_j] = A[i][j];
                        subMatrix_j++;
                    }
                    else
                    {
                        sign = pow(-1, i + j);
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

std::vector<double> mult(int size, std::vector<std::vector<double>> A, std::vector<double> b)
{
    // multiplication entre une matrice et un vecteur
    std::vector<double> Ab(size);

    for (int i = 0; i < size; i++)
    {
        double sum = 0;
        for (int j = 0; j < size; j++)
        {
            sum += A[i][j] * b[j];
        }
        Ab[i] = sum;
    }

    return Ab;
}

std::vector<std::vector<double>> EqSystemSolver::Transpose() const
{
    // les lignes de la matrice A deviennent colonnes et les colonnes deviennent les lignes
    std::vector<std::vector<double>> matrixT(size, std::vector<double>(size));
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrixT[i][j] = A[j][i];
        }
    }

    return matrixT;
}

std::vector<std::vector<double>> getSubMatrix(std::vector<std::vector<double>> A, int row, int col, int size)
{
    int subMatrix_i = 0, subMatrix_j = 0;
    std::vector<std::vector<double>> subMatrix(size - 1, std::vector<double>(size - 1));

    // on retourne les sous matrices de A selon la ligne et la colonne
    for (int i = 0; i < size; i++)
    {
        if (i != row)
        {
            subMatrix_j = 0;
            for (int j = 0; j < size; j++)
            {
                if (j != col)
                {
                    subMatrix[subMatrix_i][subMatrix_j] = A[i][j];
                    subMatrix_j++;
                }
            }
            subMatrix_i++;
        }
    }

    return subMatrix;
}

double EqSystemSolver::Determinant() const
{
    double det = 0.0;
    std::vector<std::vector<double>> subMatrix(size, std::vector<double>(size)); // submatrix will contains a part of the matrix
    int sign = -1;                                                               // sign of covariants (1 for + and -1 for -)
    int subMatrix_i, subMatrix_j;

    if (size == 2)
    {
        return ((A[0][0] * A[1][1]) - (A[1][0] * A[0][1])); //
    }
    else
    {

        for (int k = 0; k < size; k++)
        {
            subMatrix_i = 0;
            for (int i = 1; i < size; i++)
            {
                subMatrix_j = 0;
                for (int j = 0; j < size; j++)
                {
                    if (j != k)
                    {
                        subMatrix[subMatrix_i][subMatrix_j] = A[i][j];
                        subMatrix_j++;
                    }
                    else
                    {
                        sign = pow(-1, i + j);
                    }
                }
                subMatrix_i++;
            }

            det = det + (sign * A[0][k] * calculateDeterminant(size - 1, subMatrix));
        }
    }

    return det;
}

std::vector<double> EqSystemSolver::SolverGauss()
{
    std::vector<double> x(size);
    std::vector<std::vector<double>> matrix(size, std::vector<double>(size + 1));
    double ratio;
    matrix = AugmentedMatrix(); // on détermine la matrix augmenté grace à une concaténation de A et b : on obtient une matrice nxm avec m=n+1
    for (int i = 0; i < size; i++)
    {

        if (matrix[i][i] == 0.0)
        {
            std::cout << "Mathematical Error!"; // On s'assure qu'aucun pivot n'est nul
            exit(0);
        }

        // triangularisation de la matrice par la méthode de pivot
        for (int j = i + 1; j < size; j++)
        {
            ratio = matrix[j][i] / matrix[i][i];

            for (int k = 0; k <= size; k++)
            {
                matrix[j][k] = matrix[j][k] - ratio * matrix[i][k];
            }
        }
    }

    // détermination des solutions par la méthode de substitution
    x[size - 1] = matrix[size - 1][size] / matrix[size - 1][size - 1]; // on commence par la base du triangle et nous obtenons une valeur de la solution
    for (int i = size - 2; i >= 0; i--)
    {
        x[i] = matrix[i][size];
        for (int j = i + 1; j <= size; j++)
        {
            x[i] = x[i] - matrix[i][j] * x[j];
        }
        x[i] = x[i] / matrix[i][i];
    }

    return x;
}

std::vector<std::vector<double>> EqSystemSolver::AugmentedMatrix()
{
    std::vector<std::vector<double>> matrix(size, std::vector<double>(size + 1));
    // concatenation de la matrice A au vecteur b
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size + 1; j++)
        {
            if (j == size)
                matrix[i][j] = b[i];
            else
                matrix[i][j] = A[i][j];
        }
    }

    return matrix;
}

std::vector<double> EqSystemSolver::SolverInverseMatrix()
{
    std::vector<std::vector<double>> invMatrix(size, std::vector<double>(size));
    std::vector<std::vector<double>> matrixT(size, std::vector<double>(size));
    std::vector<double> x(size);

    int det = Determinant(); // etape 1 calculer le déterminant
    int sign = -1;

    if (det != 0)
    {
        matrixT = Transpose(); // etape 2: transposer la matrice
        // etape 3: déterminer la matrice adjointe en calculant les cofacteurs
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                sign = pow(-1, i + j); // le signe du cofacteur
                // std::cout<<"i= "<<i<<" j=" <<j << " sign"<<sign<<std::endl;
                invMatrix[i][j] = (sign * calculateDeterminant(size - 1, getSubMatrix(matrixT, i, j, size))) / det; // les valeurs de la matrice inverse en multipliant les valeurs de la matrice adjointe par l'inverse du determinant;
            }
        }

        x = mult(size, invMatrix, b); // la solution est la multiplication de la matrice inverse par b
        return x;
    }
    else
    {
        std::cout << "Calcul impossible\n";
        exit(0);
    }
}

std::vector<std::vector<double>> EqSystemSolver::GetA() const
{
    return A;
}

std::vector<double> EqSystemSolver::GetB() const
{
    return b;
}