#include <iostream>
#include <H5Cpp.h>

int main() {
    // Define the 2D array
    int data[3][3] = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    // Create an HDF5 file for writing
    H5::H5File file("matrix_example.h5", H5F_ACC_TRUNC);

    // Define the dimensions of the dataset
    hsize_t dims[2] = {3, 3};

    // Create a dataspace for the dataset
    H5::DataSpace dataspace(2, dims);

    // Create a dataset in the file
    H5::DataSet dataset = file.createDataSet("matrix_dataset", H5::PredType::STD_I32LE, dataspace);

    // Write the data to the dataset
    dataset.write(data, H5::PredType::STD_I32LE);

    std::cout << "Data written to matrix_dataset in matrix_example.h5" << std::endl;

    return 0;
}
