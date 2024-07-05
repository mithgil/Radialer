#include <iostream>
#include <H5Cpp.h>

int main() {
    // Open the HDF5 file for reading
    H5::H5File file("matrix_example.h5", H5F_ACC_RDONLY);

    // Open the dataset
    H5::DataSet dataset = file.openDataSet("matrix_dataset");

    // Get the dataspace of the dataset
    H5::DataSpace dataspace = dataset.getSpace();

    // Get the number of dimensions in the dataspace
    int numDims = dataspace.getSimpleExtentNdims();

    // Get the dimensions of the dataset
    hsize_t dims[numDims];
    dataspace.getSimpleExtentDims(dims, nullptr);

    // Read the data from the dataset
    int data[dims[0]][dims[1]];
    dataset.read(data, H5::PredType::STD_I32LE);

    // Display the read data
    std::cout << "Read data from matrix_dataset:" << std::endl;
    for (hsize_t i = 0; i < dims[0]; ++i) {
        for (hsize_t j = 0; j < dims[1]; ++j) {
            std::cout << data[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}
