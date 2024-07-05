/* 

This cpp script is for calculating 3D focal field of radial mode

    1. define numerical calculations for complex fields
    2. Run calculation according to the 3d volume
    3. write those data into files

@ Tim
update: 2022/5/22

*/

#include <iostream> // cin cout
#include <fstream> // write files
#include <bits/stdc++.h> //
#include <math.h> // math: sqrt() pow()
#include <complex> 
// #include <mpi.h>

using namespace std;

// Initialize arrays for restoring ex, ey, ez

/*

#define dim1 81
#define dim2 81
#define dim3 81

complex<double> Xarr [dim1][dim2][dim3];
complex<double> Yarr [dim1][dim2][dim3];
complex<double> Zarr [dim1][dim2][dim3];

float wvl = 0.636; // um

float n1 = 1.0;
float n2 = 1.0;

float k = 2*M_PI/wvl; // 1/um

float NA = 0.998; 

float thetaM = asin(NA)*(M_PI /180); // maximal angle according to the NA

// double hwindow = 1.0; // in um for half window 

// Apodization function: pupil filter

double fw(float thetaM, double x){
    double f0 = 1.0;
    
    double apod = exp(-1*pow(f0,-2)*pow(sin(x)/sin(thetaM),2));
    return apod;
}
*/

// total number of steps forthe following numerical integration 
int NumStep = 250;

// Radial mode complex form of x,y,z components: (ex,ey,ez)

complex<double> ex(float wavelength, float NA, float e0, float x, float y, float z){

    float k = 2*M_PI/wavelength;

    double f0 = 1.0;

    float thetaM = asin(NA);

    // initialize two integrals
    double integRe = 0.0, integIm = 0.0;
    
    // For calculating numberical integration (mid-point method)
    double step = thetaM/ NumStep;  // width of each small rectangle

    for (int i = 0; i < NumStep; i ++) {

        // area += f((i + 0.5) * step) * step; // sum up each small rectangle
        // Apodization function: (pupil filter) exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))
        integRe += e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
        integIm += e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
    }

    integRe = integRe*cos(atan2(y,x));
    integIm = integIm*cos(atan2(y,x));
    

    return {integRe, integIm}; // To create an object you must use curly-braces {}
    
}

complex<double> ey(float wavelength, float NA, float e0, float x, float y, float z){

    float k = 2*M_PI/wavelength;

    double f0 = 1.0;

    float thetaM = asin(NA);
    
    double integRe = 0.0, integIm = 0.0;
    
    double step = thetaM/ NumStep;  // width of each small rectangle

    for (int i = 0; i < NumStep; i ++) {

        integRe += e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
        integIm += e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
    }
    
    integRe = integRe*sin(atan2(y,x));
    integIm = integIm*sin(atan2(y,x));

    return {integRe, integIm};
    
}   

complex<double> ez(float wavelength, float NA, float e0, float x, float y, float z){
    
    float k = 2*M_PI/wavelength;

    double f0 = 1.0;

    float thetaM = asin(NA);

    double integRe = 0.0, integIm = 0.0;
    
    double step = thetaM/ NumStep;  // width of each small rectangle

    for (int i = 0; i < NumStep; i ++) {

        integRe += -2*e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*cyl_bessel_j(0, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
        integIm += 2*e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*cyl_bessel_j(0, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
    }
    
    return {integRe, integIm};

}


complex<double>*** CreateArray(int comp, int dim, float wavelengthIn, float amp0, float NAIn, float hwindow){


    int width = dim;
    int length = dim;
    int height = dim;
    // create a 3d array
    complex<double>*** arr = new complex<double>**[width];

    for(int i = 0; i < width; ++i){
        arr[i] = new complex<double>*[length];
        for(int j = 0; j < length; ++j){
            arr[i][j] = new complex<double>[height];
            
        }
    }

    cout << "Calculation for " ;
    
    clock_t start, end;
    start = clock();

    cout << setprecision(10) << fixed;
    switch (comp) {
        case 1:
            cout << "ex ..."<<'\n';
            for(int i=0; i < dim; i++){
                for(int j=0; j < dim; j++){
                    for(int k=0; k < dim; k++){
                        //  return ex
                        arr[i][j][k] = ex(wavelengthIn,NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim);
                        
                    }
                }
            }
            break;

        case 2:
            cout << "ey ..."<<'\n';     
            for(int i=0; i < dim; i++){
                for(int j=0; j < dim; j++){
                    for(int k=0; k < dim; k++){
                        //  return ey
                        arr[i][j][k] = ey(wavelengthIn,NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim);
                        
                    }
                }
            }
            break;
        case 3:
            cout << "ez ..."<<'\n';
            for(int i=0; i < dim; i++){
                for(int j=0; j < dim; j++){
                    for(int k=0; k < dim; k++){
                        //  return ez
                        arr[i][j][k] = ez(wavelengthIn,NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim);
                        
                    }
                }
            }
            break;
               

    }
    
    end = clock();

    cout << " is finished!" <<endl;
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << " Time taken by the calculation is : " << fixed 
            << time_taken << setprecision(6);
    cout << " sec " << endl;

    return arr;

}



int main(){
    /* clock_t clock(void) returns the number of clock ticks
       elapsed since the program was launched.To get the number 
       of seconds used by the CPU, you will need to divide by 
       CLOCKS_PER_SEC.where CLOCKS_PER_SEC is 1000000 on typical
       32 bit system.  */
    
    cout << setprecision(15) << fixed;

    // Set: int comp, int dim, float wavelengthIn, float amp0, float NAIn, float hwindow
    int dimSet = 81;
    float wavelengthSet = 0.636; // um
    float ampSet = 1.0;
    float NASet = 0.998;
    float hwindowSet = 1.0;

    complex<double>*** exArr = CreateArray(1, dimSet, wavelengthSet,ampSet, NASet, hwindowSet);
    complex<double>*** eyArr = CreateArray(2, dimSet, wavelengthSet,ampSet, NASet, hwindowSet);
    complex<double>*** ezArr = CreateArray(3, dimSet, wavelengthSet,ampSet, NASet, hwindowSet);

    cout<< "ex="<<exArr[20][5][5] <<endl;
    cout<< "ey="<<eyArr[0][5][1] <<endl;
    cout<< "ez="<<ezArr[1][5][2] <<endl;


    // Write those arrays into a txt file:
    
    ofstream myEx;
    ofstream myEy;
    ofstream myEz;

    myEx.open ("FocalField Radial Ex.txt");
    myEy.open ("FocalField Radial Ey.txt");
    myEz.open ("FocalField Radial Ez.txt");
    
    // myEx << "#Focal field of the radial polarization.\n";
    cout << "Writing into txt files ..." <<endl;
    for(int i=0; i < dimSet; i++){
        for(int j=0; j < dimSet; j++){
            for(int k=0; k < dimSet; k++){
                
                // writing Xarr Yarr Zarr
                
                myEx << exArr[i][j][k];
                myEx << "  "; // spacing between different complex numbers
                myEy << eyArr[i][j][k]; 
                myEy << "  ";
                myEz << ezArr[i][j][k]; 
                myEz << "  ";

            }
            
        }
        
    }

    myEx.close();
    myEy.close();
    myEz.close();
    cout << "done!." <<endl;

    //allocated memory must be deleted:
    for(int i = 0; i < dimSet; ++i){
        
        for( int j = 0; j<dimSet;j++){
            delete[] exArr[i][j];
            delete[] eyArr[i][j];
            delete[] ezArr[i][j];
        }
        delete[] exArr[i];
        delete[] eyArr[i];
        delete[] ezArr[i];
    }

    delete[] exArr;
    delete[] eyArr;
    delete[] ezArr;

    exArr = NULL;
    eyArr = NULL;
    ezArr = NULL;

    cout << "Memory for storing arrays is deleted." <<endl;

    return 0;
}

