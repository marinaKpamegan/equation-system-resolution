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