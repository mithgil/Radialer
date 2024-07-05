/*
Collected functions in the calculations in a header
@ Tim
update: 2024/7/3

1. factorial
2. Bessel of the first kind
3. Numerical form of the Radial Mode
    - Mid-point rule
    - Gauss quadrature
4. create array for 
5. exporting array 
6. delete the array

*/

#include <vector>
#include "radial.h"
#include <cmath> // math: sqrt() pow()
#include <complex> 
#include <iostream> // std
#include <iomanip>  // std::setprecision
#include <fstream> //output

#include "H5Cpp.h"

#include <boost/math/quadrature/gauss.hpp> // gauss quadrature

using namespace std;

//factorial
long fact(int n){

	long factorial = 1.0;

	for(int i=1; i<=n ; ++i){
		factorial*=i;
	}

	return factorial;
}

// bessel order of the first kind
// nth order besselJ at x position with nonzero integer n & n>=0
double besselj(int n, double x){
	double sum = 0;
	int M=20; // note this will make factorial overflow

	for (int k=0; k<=M; ++k){
		sum += pow(-1,k)*pow(0.5*x,2*k+n)/fact(k)/fact(k+n);
	}

	return sum;
}

// Radial mode complex form of x,y,z components: (ex,ey,ez) using midpoint method Riemann integral
// the error estimation will be (b-a)/24N^2*f''(eta) where eta in [a,b]

const int NumStep = 60; 

Complex midpointRad(char comp, const Params& params){
    
    float k = 2*M_PI/params.wavelength;
    // filling factor
    double f0 = params.fillfactor;

    float thetaM = asin(params.NA);

    // initialize two integrals
    double integRe = 0.0, integIm = 0.0;
    // For calculating numberical integration 
    double step = thetaM/ NumStep;  

    // `comp` dnotes x, y, or z
    switch(comp) {

        case 'x':
            // using Midpoint rule Riemann Sum
            
            for (int i = 0; i <= NumStep; i ++) {

                // Apodization function: (pupil filter) exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))

                integRe += step* params.e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*sin(2*(i + 0.5) * step)*besselj(1, k*sqrt(params.x*params.x+params.y*params.y)*sin((i + 0.5) * step))*cos(k*params.z*cos((i + 0.5) * step));
                integIm += step* params.e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*sin(2*(i + 0.5) * step)*besselj(1, k*sqrt(params.x*params.x+params.y*params.y)*sin((i + 0.5) * step))*sin(k*params.z*cos((i + 0.5) * step));
            }

            integRe *= cos(atan2(params.y,params.x));
            integIm *= cos(atan2(params.y,params.x));
            
            return {integRe, integIm}; // To create an object you must use curly-braces {}

            break;
        case 'y':
            for (int i = 0; i <= NumStep; i ++) {

                integRe += step* params.e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*sin(2*(i + 0.5) * step)*besselj(1, k*sqrt(params.x*params.x+params.y*params.y)*sin((i + 0.5) * step))*cos(k*params.z*cos((i + 0.5) * step));
                integIm += step* params.e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*sin(2*(i + 0.5) * step)*besselj(1, k*sqrt(params.x*params.x+params.y*params.y)*sin((i + 0.5) * step))*sin(k*params.z*cos((i + 0.5) * step));
            }

            integRe *= sin(atan2(params.y,params.x));
            integIm *= sin(atan2(params.y,params.x));
            
            return {integRe, integIm}; 

            break;
        case 'z':
            
            for (int i = 0; i <= NumStep; i ++) {

                integRe += -2*step*params.e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*besselj(0, k*sqrt(params.x*params.x+params.y*params.y)*sin((i + 0.5) * step))*sin(k*params.z*cos((i + 0.5) * step));
                integIm += 2*step*params.e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*besselj(0, k*sqrt(params.x*params.x+params.y*params.y)*sin((i + 0.5) * step))*cos(k*params.z*cos((i + 0.5) * step));
            }
            
            return {integRe, integIm};
            break;

        default:
            cout<<"Incorrect component! (x,y,z).";
            return {0,0};
    }
}

// Radial mode complex form of x,y,z components: (ex,ey,ez) using Gauss quadrature
// form for exy

double xyrad_r(double alpha, const Params& params) {

    float k = 2*M_PI/params.wavelength;
    // filling factor
    float f0 = params.fillfactor;

    float thetaM = asin(params.NA);

    //use Bessel-Gauss beam 
    return params.e0*exp(-1*pow(f0,-2)*pow(sin(alpha)/sin(thetaM),2))*besselj(1,2*f0*sin(alpha)/sin(thetaM))*
            sqrt(cos(alpha))*sin(2*alpha)*besselj(1, k*sqrt(params.x*params.x+params.y*params.y)*sin(alpha))*cos(k*params.z*cos(alpha));
}


double xyrad_i(double alpha, const Params& params) {
    float k = 2*M_PI/params.wavelength;
    // filling factor
    float f0 = params.fillfactor;

    float thetaM = asin(params.NA);

    return params.e0*exp(-1*pow(f0,-2)*pow(sin(alpha)/sin(thetaM),2))*besselj(1,2*f0*sin(alpha)/sin(thetaM))*
            sqrt(cos(alpha))*sin(2*alpha)*besselj(1, k*sqrt(params.x*params.x+params.y*params.y)*sin(alpha))*sin(k*params.z*cos(alpha));
}

// form fo ez
double zrad_r(double alpha, const Params& params) {
    float k = 2*M_PI/params.wavelength;
    // filling factor
    float f0 = params.fillfactor;

    float thetaM = asin(params.NA);

    return -2*params.e0*exp(-1*pow(f0,-2)*pow(sin(alpha)/sin(thetaM),2))*besselj(1,2*f0*sin(alpha)/sin(thetaM))*
            sqrt(cos(alpha))*pow(sin(alpha),2)*besselj(0, k*sqrt(params.x*params.x+params.y*params.y)*sin(alpha))*sin(k*params.z*cos(alpha));
}

double zrad_i(double alpha, const Params& params) {
    float k = 2*M_PI/params.wavelength;
    // filling factor
    float f0 = params.fillfactor;

    float thetaM = asin(params.NA);

    return 2*params.e0*exp(-1*pow(f0,-2)*pow(sin(alpha)/sin(thetaM),2))*besselj(1,2*f0*sin(alpha)/sin(thetaM))*
            sqrt(cos(alpha))*pow(sin(alpha),2)*besselj(0, k*sqrt(params.x*params.x+params.y*params.y)*sin(alpha))*cos(k*params.z*cos(alpha));
}


using namespace boost::math::quadrature;

Complex gaussRad(char comp, const Params& params){
        
    const std::size_t n = 10;

    gauss<double, n> integrator_r;
    gauss<double, n> integrator_i;

    double integRe = 0.0, integIm = 0.0;
    float thetaM = asin(params.NA);

    // Create a lambda function that captures the parameters a and b
    auto txyrad_i = [&params](double x) {
        return xyrad_i(x, params);
    };
    
    auto txyrad_r = [&params](double x) {
        return xyrad_r(x, params);
    };

    auto tzrad_r = [&params](double x) {
        return zrad_r(x, params);
    };

    auto tzrad_i = [&params](double x) {
        return zrad_i(x, params);
    };

    // `comp` dnotes x, y, or z
    switch(comp) {

        case 'x':
                                           
            integRe = integrator_r.integrate(txyrad_r, 0, thetaM)*cos(atan2(params.y,params.x));
            integIm  = integrator_i.integrate(txyrad_i, 0, thetaM)*cos(atan2(params.y,params.x));
            
            return {integRe, integIm}; // To create an object you must use curly-braces {}

            break;
        case 'y':

                                
            integRe = integrator_r.integrate(txyrad_r, 0, thetaM)*sin(atan2(params.y,params.x));
            integIm  = integrator_i.integrate(txyrad_i, 0, thetaM)*sin(atan2(params.y,params.x));
            
            return {integRe, integIm}; 

            break;
        case 'z':

            integRe = integrator_r.integrate(tzrad_r, 0, thetaM);
            integIm = integrator_i.integrate(tzrad_i, 0, thetaM);

            return {integRe, integIm};
            break;

        default:
            cout<<"Incorrect component! (x,y,z).";
            return {0,0};
    }
}

//  to create a 4d array (3 components)
Complex ****CreateArray(int dim, float wavelengthIn, float fillfactorIn, float amp0, float NAIn, float hwindow){

    int cube = 3;
    int width = dim;
    int length = dim;
    int height = dim;
    // create a 4d array
    complex<double> ****arr = new complex<double>***[cube];

    for(int l = 0; l < cube; ++l){
        arr[l] = new complex<double>**[width];
        for(int i = 0; i < width; ++i){
            arr[l][i] = new complex<double>*[length];
            for(int j = 0; j < length; ++j){
                for(int j = 0; j < length; ++j){
                    arr[l][i][j] = new complex<double>[height];
                    
                }
            }
        }
    }

    cout << "Radial Mode calculations: " ;
    
    clock_t start, end;
    start = clock();

    cout << std::setprecision(10) << fixed;
    
            cout << "ex ey & ez ..."<<'\n';
            for(int i=0; i < dim; i++){
                for(int j=0; j < dim; j++){
                    for(int k=0; k < dim; k++){
                        
                        Params params = {wavelengthIn, fillfactorIn, NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim};
                        //  return ex, ey, ez as a 4d array
                        arr[0][i][j][k] = gaussRad('x',params); // gaussain-quadrature
                        arr[1][i][j][k] = gaussRad('y',params); 
                        arr[2][i][j][k] = gaussRad('z',params);

                    }
                }
            }
            
    end = clock();      
    
    cout << " is finished!" <<endl;
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << " Time for the calculation is : " << fixed 
            << time_taken << std::setprecision(6);
    cout << " sec " << endl;

    return arr;

}

void DeleteArray(Complex ****array, int dim1Size,int dim2Size,int dim3Size) {
    
    // Clean up allocated memory
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

    cout << "Array by pointers is cleaned!" <<endl;
}

void flatten4DArray(std::complex<double>**** source, std::vector<std::complex<double>>& dest, \
                    size_t dim1, size_t dim2, size_t dim3, size_t dim4) {
    dest.resize(dim1 * dim2 * dim3 * dim4);
    size_t index = 0;
    for (size_t i = 0; i < dim1; ++i) {
        for (size_t j = 0; j < dim2; ++j) {
            for (size_t k = 0; k < dim3; ++k) {
                for (size_t l = 0; l < dim4; ++l) {
                    dest[index++] = source[i][j][k][l];
                }
            }
        }
    }
}

const std::string FILE_NAME = "4d_complex_array.h5";

void writeToHDF5(Complex**** array, int dim1Size, int dim2Size, int dim3Size, int dim4Size) {

   H5::H5File file(FILE_NAME, H5F_ACC_TRUNC);

    // Define the dimensions of the dataset
    hsize_t dims[4] = {static_cast<hsize_t>(dim1Size), static_cast<hsize_t>(dim2Size),
                       static_cast<hsize_t>(dim3Size), static_cast<hsize_t>(dim4Size)};

    // Create a dataspace for the dataset
    H5::DataSpace dataspace(4, dims);

    // Create a compound datatype for complex<double>
    H5::CompType complexType(sizeof(std::complex<double>));
    complexType.insertMember("real", 0, H5::PredType::NATIVE_DOUBLE);
    complexType.insertMember("imag", sizeof(double), H5::PredType::NATIVE_DOUBLE);

    // Create a dataset in the file with complex<double> datatype
    H5::DataSet dataset = file.createDataSet("complex_dataset", complexType, dataspace);

    // Write the data to the dataset
    dataset.write(array, complexType);
    
    std::cout << "Data written into h5." << std::endl;

}

void writeTotxt(Complex**** array, int dim1Size, int dim2Size, int dim3Size, int dim4Size, const std::string& filename) {
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

    std::cout << "Data written into txt." << std::endl;
}
