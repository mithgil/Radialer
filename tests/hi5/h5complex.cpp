/* An example code to write/read comple of double 3D array into hdf5
2024/6/29
 
 @Tim

*/

#include <iostream>
#include <complex>
#include <H5Cpp.h>

int main() {
    // Define the 3D array with complex<double> data
    std::complex<double> data[2][3][4] = {
        {{{1.1, 2.2}, {3.3, 4.4}, {5.5, 6.6}, {7.7, 8.8}},
        {{9.9, 10.1}, {11.11, 12.12}, {13.13, 14.14}, {15.15, 16.16}},
        {{17.17, 18.18}, {19.19, 20.2}, {21.21, 22.22}, {23.23, 24.24}}},
        {{{1.1, 2.2}, {3.3, 4.4}, {5.5, 6.6}, {7.7, 8.8}},
        {{9.9, 10.1}, {11.11, 12.12}, {13.13, 14.14}, {15.15, 16.16}},
        {{17.17, 18.18}, {19.19, 20.2}, {21.21, 22.22}, {23.23, 24.24}}}

    };

    // Create an HDF5 file for writing
    H5::H5File file("complex_3d_example.h5", H5F_ACC_TRUNC);

    // Define the dimensions of the dataset
    hsize_t dims[3] = {2, 3, 4};

    // Create a dataspace for the dataset
    H5::DataSpace dataspace(3, dims);

    // Create a compound datatype for complex<double>
    H5::CompType complexType(sizeof(std::complex<double>));
    complexType.insertMember("real", 0, H5::PredType::NATIVE_DOUBLE);
    complexType.insertMember("imag", sizeof(double), H5::PredType::NATIVE_DOUBLE);

    // Create a dataset in the file with complex<double> datatype
    H5::DataSet dataset = file.createDataSet("complex_dataset", complexType, dataspace);

    // Write the data to the dataset
    dataset.write(data, complexType);

    std::cout << "Data written to complex_dataset in complex_3d_example.h5" << std::endl;

    // Open the HDF5 file for reading
    H5::H5File read_file("complex_3d_example.h5", H5F_ACC_RDONLY);

    // Open the dataset
    H5::DataSet read_dataset = read_file.openDataSet("complex_dataset");

    // Read the data from the dataset
    std::complex<double> read_data[dims[0]][dims[1]][dims[2]];
    read_dataset.read(read_data, complexType);

    // Display the read data
    std::cout << "Read data from complex_dataset:" << std::endl;
    for (hsize_t i = 0; i < dims[0]; ++i) {
        for (hsize_t j = 0; j < dims[1]; ++j) {
            for (hsize_t k = 0; k < dims[2]; ++k) {
                std::cout << "(" << read_data[i][j][k].real() << ", " << read_data[i][j][k].imag() << ") ";
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    return 0;
}
