

#include "EqSystemSolver.hpp"
#include <cmath>
#include <cassert>
#include <vector>

EqSystemSolver::EqSystemSolver(int taille, std::vector<std::vector<double>> matrix, std::vector<double> B)
{
    A = matrix;
    b = B;
    size = taille;
}

/*vector<vector<double>> EqSystemSolver::Concatenate(){

}*/
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
    std::vector<double> Ab(size);

    for (int i = 0; i < size; i++)
    {
        int sum = 0;
        for (int j = 0; j < size; j++)
        {
            sum += A[i][j] * b[j];
        }
        Ab[i] = sum;
    }

    return Ab;
}

std::vector<std::vector<double>> reverseMatrix(int size, std::vector<std::vector<double>> A)
{
    double det;
    std::vector<std::vector<double>> matrixReverse(size, std::vector<double>(size));
    det = calculateDeterminant(size, A);

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrixReverse[i][j] = 1 / det * A[i][j];
        }
    }
    return matrixReverse; // A moins 1
}

std::vector<std::vector<double>> EqSystemSolver::Transpose() const
{

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

    // std::cout << "row= " << row << " j= " << col << "\n";
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
                    // std::cout << "value at " << i << " " << j << " --> " << A[i][j] << "\n";
                    subMatrix_j++;
                }
            }
            subMatrix_i++;
        }
    }

    // todo to remove
    for (int i = 0; i < size - 1; i++)
    {
        for (int j = 0; j < size - 1; j++)
        {
            std::cout << subMatrix[i][j] << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
    return subMatrix;
}

/*std::vector<std::vector<double>> EqSystemSolver::Inverse(){
    double det;
    std::vector<std::vector<double>> matrixReverse(size, std::vector<double>(size));
    det = calculateDeterminant(size, A);
    assert(det !=0);

    for (int i=0; i<size; i++){
        for (int j=0; j<size; j++){
            matrixReverse[i][j] = 1/det * A[i][j];
        }
    }
    return matrixReverse;  // A moins 1

}*/

double EqSystemSolver::Determinant() const
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
    matrix = AugmentedMatrix();
    for (int i = 0; i < size; i++)
    {

        if (matrix[i][i] == 0.0)
        {
            std::cout << "Mathematical Error!";
            exit(0);
        }
        for (int j = i + 1; j < size; j++)
        {
            ratio = matrix[j][i] / matrix[i][i];

            for (int k = 0; k <= size; k++)
            {
                matrix[j][k] = matrix[j][k] - ratio * matrix[i][k];
            }
        }
    }

    x[size - 1] = matrix[size - 1][size] / matrix[size - 1][size - 1];
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

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size + 1; j++)
        {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << "\n";
    }
    return matrix;
}

std::vector<double> EqSystemSolver::SolverInverseMatrix()
{
    std::vector<std::vector<double>> invMatrix(size, std::vector<double>(size));
    std::vector<std::vector<double>> matrixT(size, std::vector<double>(size));
    std::vector<double> x(size);

    matrixT = Transpose();
    int det = Determinant();

    if (det != 0)
    {
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                invMatrix[i][j] = (pow(-1, i + j) * calculateDeterminant(size - 1, getSubMatrix(matrixT, i, j, size)))/det;
            }
        }

        x = mult(size, invMatrix, b);

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                std::cout << invMatrix[i][j] << " ";
            }
            std::cout << "\n";
        }
    }
    return x;
}