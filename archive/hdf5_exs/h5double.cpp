#include <iostream>
#include <H5Cpp.h>

int main() {
    // Define the 2D array with double data
    double data[3][3] = {
        {1.1, 2.2, 3.3},
        {4.4, 5.5, 6.6},
        {7.7, 8.8, 9.9}
    };

    // Create an HDF5 file for writing
    H5::H5File file("matrix_double_example.h5", H5F_ACC_TRUNC);

    // Define the dimensions of the dataset
    hsize_t dims[2] = {3, 3};

    // Create a dataspace for the dataset
    H5::DataSpace dataspace(2, dims);

    // Create a dataset in the file with double datatype
    H5::DataSet dataset = file.createDataSet("matrix_dataset", H5::PredType::NATIVE_DOUBLE, dataspace);

    // Write the data to the dataset
    dataset.write(data, H5::PredType::NATIVE_DOUBLE);

    std::cout << "Data written to matrix_dataset in matrix_double_example.h5" << std::endl;

    // Open the HDF5 file for reading
    H5::H5File read_file("matrix_double_example.h5", H5F_ACC_RDONLY);

    // Open the dataset
    H5::DataSet read_dataset = read_file.openDataSet("matrix_dataset");

    // Read the data from the dataset
    double read_data[dims[0]][dims[1]];
    read_dataset.read(read_data, H5::PredType::NATIVE_DOUBLE);

    // Display the read data
    std::cout << "Read data from matrix_dataset:" << std::endl;
    for (hsize_t i = 0; i < dims[0]; ++i) {
        for (hsize_t j = 0; j < dims[1]; ++j) {
            std::cout << read_data[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
