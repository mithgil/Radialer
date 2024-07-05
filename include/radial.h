/*
Names of collected functions in the calculations in a header

@ Tim

1. factorial
2. Bessel of the first kind
3. Numerical form of the Radial Mode
4. Create an array of Radial Mode

update: 2024/7/1

*/
#ifndef RADIAL_H
#define RADIAL_H

#include <cmath> // math: sqrt() pow()
#include <complex> 
#include <vector>
#include <string>
#include <iostream> // std
#include <iomanip>      // std::setprecision
using namespace std;
//factorial
long fact(int n);

// nth order besselJ at x position with nonzero integer n & n>=0
double besselj(int n, double x);

// total number of steps for the following numerical integration 
extern const int NumStep; // will be determined in the radial.cpp

typedef std::complex<double> Complex;

struct Params {
    float wavelength;
    float fillfactor;
    float NA;
    float e0;
    float x;
    float y;
    float z;
};

// Radial mode complex form of x,y,z components: (ex,ey,ez)
Complex midpointRad(char comp, const Params& params);

double xyrad_i(double alpha, const Params& params);

double xyrad_i(double alpha, const Params& params);

double zrad_r(double alpha, const Params& params);

double zrad_i(double alpha, const Params& params);

Complex gaussRad(char comp, const Params& params);

// Creating Radial mode by pointers
Complex ****CreateArray(int dim, float wavelengthIn, float fillfactorIn, float amp0, float NAIn, float hwindow);
// Delete pointer array
void DeleteArray(Complex ****array, int dim1Size, int dim2Size, int dim3Size);
// flatten array for npy
void flatten4DArray(std::complex<double>**** source, std::vector<std::complex<double>>& dest, size_t dim1, size_t dim2, size_t dim3, size_t dim4);
// writing array to hdf5
void writeToHDF5(Complex **** array, int dim1Size, int dim2Size, int dim3Size, int dim4Size);
// writing array to txt
void writeTotxt(Complex **** array, int dim1Size, int dim2Size, int dim3Size, int dim4Size, const std::string& filename);


#endif
