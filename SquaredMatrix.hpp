#ifndef SQUAREDMATRIXHEADERDEF
#define SQUAREDMATRIXHEADERDEF

#include <iostream>

class SquaredMatrix
{
private:
    int nRows;
    int nCols;
public:
    SquaredMatrix();
    SquaredMatrix(int m, int n);
    int CalculateDeterminant(SquaredMatrix& M) const;
    SquaredMatrix Transpose(SquaredMatrix& M) const;
    SquaredMatrix Triangularize(SquaredMatrix& M) const;
    SquaredMatrix Mult(SquaredMatrix& M, double* b) const;
    friend std::ostream& operator<<(std::ostream& output,
    const SquaredMatrix& M);
}

#endif