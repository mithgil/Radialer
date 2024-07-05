/* An example code to write/read comple of double 3D array into hdf5
2024/7/1
 
 @Tim

*/
#include <iostream>
#include <complex>
#include <H5Cpp.h>

int main() {
    // Open the HDF5 file for reading
    H5::H5File read_file("../build/4d_complex_array.h5", H5F_ACC_RDONLY);

    // Open the dataset
    H5::DataSet read_dataset = read_file.openDataSet("complex_dataset");

    // Define the compound datatype for std::complex<double>
    H5::CompType complexType(sizeof(std::complex<double>));
    complexType.insertMember("real", 0, H5::PredType::NATIVE_DOUBLE);
    complexType.insertMember("imag", sizeof(double), H5::PredType::NATIVE_DOUBLE);

    int dimset = 81;
    hsize_t dims[4] = {3, static_cast<hsize_t>(dimset), static_cast<hsize_t>(dimset), static_cast<hsize_t>(dimset)};
    
    // Read the data from the dataset
    std::complex<double> read_data[3][dimset][dimset][dimset];
    read_dataset.read(read_data, complexType);

    // Display the read data
    std::cout << "Read data from complex_dataset:" << std::endl;
    
    for (hsize_t i = 0; i < 3; ++i) {
        for (hsize_t j = 0; j < 10; ++j) {  // Adjust the range for display
            for (hsize_t k = 0; k < 10; ++k) {
                for (hsize_t l = 0; l < 10; ++l) {
                    std::cout << "(" << read_data[i][j][k][l].real() << ", " << read_data[i][j][k][l].imag() << ") ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    return 0;
}

