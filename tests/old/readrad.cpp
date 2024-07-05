
/*
 To read h5 file in C++ seems buggy

@ Tim
update: 2023/12/5

*/

#include <iostream>
#include <string>
#include <vector>
#include <H5Cpp.h>
#include <complex>

using namespace std;

const std::string FILE_NAME = "complex_rad.h5";
const std::string DATASET_NAME = "complex_dataset";

int main() {
    // Create a compound datatype for complex<double>
    H5::CompType complexType(sizeof(std::complex<double>));
    complexType.insertMember("real", 0, H5::PredType::NATIVE_DOUBLE);
    complexType.insertMember("imag", sizeof(double), H5::PredType::NATIVE_DOUBLE);

   
    hsize_t dims[4] = {3, 81, 81, 81};
    
    // Read the data from the dataset
    try {
        H5::H5File read_file("complex_rad.h5", H5F_ACC_RDONLY);

        // Open the dataset
        H5::DataSet read_dataset = read_file.openDataSet("complex_dataset");
        
        // Read the data from the dataset
        
        std::complex<double> read_data[dims[0]][dims[1]][dims[2]][dims[3]];

        read_dataset.read(read_data, complexType);

        /*
        // Display the read data
        std::cout << "Read data from complex_dataset:" << std::endl;
        
        for (hsize_t i = 0; i < static_cast<hsize_t>(dims[0]); ++i) {
            for (hsize_t j = 0; j < static_cast<hsize_t>(dims[1]); ++j) {
                for (hsize_t k = 0; k < static_cast<hsize_t>(dims[2]); ++k) {
                    std::cout << "(" << read_data[0][i][j][k].real() << ", " << read_data[0][i][j][k].imag() << ") ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
            
        }
        */
        
    } catch (const H5::Exception& e) {
        std::cerr << "HDF5 Exception: " << e.getCDetailMsg() << std::endl;
    }
    

    return 0;
}