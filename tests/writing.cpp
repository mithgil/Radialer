#include <iostream>
#include <fstream>
#include <complex>
#include <vector>

void write4DComplexArrayToFile(std::complex<double>**** array, int dim1Size, int dim2Size, int dim3Size, int dim4Size, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return;
    }

    for (int i = 0; i < dim1Size; ++i) {
        for (int j = 0; j < dim2Size; ++j) {
            for (int k = 0; k < dim3Size; ++k) {
                for (int l = 0; l < dim4Size; ++l) {
                    outFile << array[i][j][k][l].real() << "+" << array[i][j][k][l].imag() << "i ";
                }
                outFile << '\n';
            }
            outFile << '\n';
        }
        outFile << '\n';
    }

    outFile.close();
}

int main() {
    int dim1Size = 2, dim2Size = 2, dim3Size = 2, dim4Size = 2;

    // Allocate the 4D array dynamically
    std::complex<double>**** array = new std::complex<double>***[dim1Size];
    for (int i = 0; i < dim1Size; ++i) {
        array[i] = new std::complex<double>**[dim2Size];
        for (int j = 0; j < dim2Size; ++j) {
            array[i][j] = new std::complex<double>*[dim3Size];
            for (int k = 0; k < dim3Size; ++k) {
                array[i][j][k] = new std::complex<double>[dim4Size];
                for (int l = 0; l < dim4Size; ++l) {
                    array[i][j][k][l] = std::complex<double>(i + j + k + l, i - j - k + l); // Example values
                }
            }
        }
    }

    // Write the array to a file
    write4DComplexArrayToFile(array, dim1Size, dim2Size, dim3Size, dim4Size, "4d_complex_array.txt");

    // Free the allocated memory
    for (int i = 0; i < dim1Size; ++i) {
        for (int j = 0; j < dim2Size; ++j) {
            for (int k = 0; k < dim3Size; ++k) {
                delete[] array[i][j][k];
            }
            delete[] array[i][j];
        }
        delete[] array[i];
    }
    delete[] array;

    return 0;
}
