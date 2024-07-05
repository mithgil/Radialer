#include <iostream>
#include <vector>
#include <complex>
#include "cnpy.h"

int main() {
    // Define a 4D array of std::complex<double>
    size_t dim1 = 3, dim2 = 81, dim3 = 81, dim4 = 81;
    std::vector<std::complex<double>> data(dim1 * dim2 * dim3 * dim4);

    // Fill the array with some data
    for (size_t i = 0; i < dim1 * dim2 * dim3 * dim4; ++i) {
        data[i] = std::complex<double>(i, -i);
    }

    // Save to a .npy file
    cnpy::npy_save("data.npy", &data[0], {dim1, dim2, dim3, dim4}, "w");

    // Load the data back from the .npy file
    cnpy::NpyArray arr = cnpy::npy_load("radial_array.npy");
    std::complex<double>* loaded_data = arr.data<std::complex<double>>();

    // Print out some of the loaded data to verify
    for (size_t i = 0; i < 10; ++i) {
        std::cout << loaded_data[i] << std::endl;
    }

    return 0;
}

