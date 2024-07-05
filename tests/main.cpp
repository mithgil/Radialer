/* 
Radialer

This cpp script is for calculating 3D focal field of radial mode

    1. define numerical calculations for complex fields as functions
    2. define 4 dimentisonal array by 4-fold pointers to run calculation according to the 3d volume and other settings
    3. write those data to txt, npy 
    4. delete memoryfor pointers

@ Tim
update: 2024/7/4

*/

#include <iostream> // cin cout-must have!
#include <fstream> // write files
#include <math.h> // math: sqrt() pow()
#include <complex> 
#include <iomanip>      // std::setprecision
#include "radial.h" // header of functions
#include <H5Cpp.h>
#include "cnpy.h" // export npy for numpy

using namespace std;

int main(){

    // Set: int comp, int dim, float wavelengthIn, float amp0, float NAIn, float hwindow
    int dimSet = 81; 
    float wavelengthSet = 0.636; // um
    float fillfactorSet = 1.0;
    float ampSet = 1.0;
    float NASet = 0.998;
    float hwindowSet = 1.0;

    // Radial field of ex, ey, ez components by using a 4D pointer array
    complex<double> ****RadArr = CreateArray(dimSet, wavelengthSet, fillfactorSet,ampSet, NASet, hwindowSet);

    // Flatten the 4D array
    std::vector<std::complex<double>> flattenedArray;
    flatten4DArray(RadArr, flattenedArray, 3, dimSet, dimSet, dimSet);

    // write to npy
    cnpy::npy_save("radial_array.npy", &flattenedArray[0], {3, dimSet, dimSet, dimSet}, "w");

    //write to txt file
    writeTotxt(RadArr,3,dimSet,dimSet,dimSet,"radial_array.txt");    
  	
    // write to hdf5 file: buggy
    //writeToHDF5(RadArr,3,dimSet,dimSet,dimSet);    
    
    //free the allocated memory
    DeleteArray(RadArr,3,dimSet,dimSet);
    
    return 0;
}
