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
*/

int NumStep = 250;

// double hwindow = 1.0; // in um for half window 

float e0 = 1.0;


// Apodization function: pupil filter

double fw(float thetaM, double x){
    double f0 = 1.0;
    
    double apod = exp(-1*pow(f0,-2)*pow(sin(x)/sin(thetaM),2));
    return apod;
}


// Radial mode complex form of x,y,z components: (ex,ey,ez)

complex<double> ex(float wavelength, float NA, float x, float y, float z){

    float k = 2*M_PI/wavelength;

    float thetaM = asin(NA)*(M_PI /180);

    // initialize two integrals
    double integRe = 0.0, integIm = 0.0;
    
    // For calculating numberical integration (mid-point method)
    double step = thetaM/ NumStep;  // width of each small rectangle

    for (int i = 0; i < NumStep; i ++) {

        // area += f((i + 0.5) * step) * step; // sum up each small rectangle
        integRe += e0*fw(thetaM,(i + 0.5) * step)*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
        integIm += e0*fw(thetaM,(i + 0.5) * step)*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
    }

    integRe = integRe*cos(atan2(y,x));
    integIm = integIm*cos(atan2(y,x));
    

    return {integRe, integIm}; // To create an object you must use curly-braces {}
    
}

complex<double> ey(float wavelength, float NA, float x, float y, float z){

    float k = 2*M_PI/wavelength;

    float thetaM = asin(NA)*(M_PI /180);
    
    double integRe = 0.0, integIm = 0.0;
    
    double step = thetaM/ NumStep;  // width of each small rectangle

    for (int i = 0; i < NumStep; i ++) {

        integRe += e0*fw(thetaM, (i + 0.5) * step)*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
        integIm += e0*fw(thetaM, (i + 0.5) * step)*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
    }
    
    integRe = integRe*sin(atan2(y,x));
    integIm = integIm*sin(atan2(y,x));

    return {integRe, integIm};
    
}   

complex<double> ez(float wavelength, float NA, float x, float y, float z){
    
    float k = 2*M_PI/wavelength;

    float thetaM = asin(NA)*(M_PI /180);

    double integRe = 0.0, integIm = 0.0;
    
    double step = thetaM/ NumStep;  // width of each small rectangle

    for (int i = 0; i < NumStep; i ++) {

        integRe += -2*e0*fw(thetaM, (i + 0.5) * step)*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*cyl_bessel_j(0, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
        integIm += 2*e0*fw(thetaM, (i + 0.5) * step)*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*cyl_bessel_j(0, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
    }
    
    return {integRe, integIm};

}

class RadField
{
    public:
    
    float wavelengthIn, NAIn, hwindow;
    int dim1, dim2, dim3;
    complex<double> Xarr [dim1][dim2][dim3];
    complex<double> Yarr [dim1][dim2][dim3];
    complex<double> Zarr [dim1][dim2][dim3];
    
    void run(){

        for(int i=0; i < dim1; i++){
            for(int j=0; j < dim2; j++){
                for(int k=0; k < dim3; k++){
                    // calculate from -hwindow to +hwindow of dimX sampling points
                    Xarr[i][j][k] = ex(wavelengthIn,NAIn,-hwindow+i*2*hwindow/dim1,-hwindow+j*2*hwindow/dim2,-hwindow+k*2*hwindow/dim3);
                    Yarr[i][j][k] = ey(wavelengthIn,NAIn,-hwindow+i*2*hwindow/dim1,-hwindow+j*2*hwindow/dim2,-hwindow+k*2*hwindow/dim3);
                    Zarr[i][j][k] = ez(wavelengthIn,NAIn,-hwindow+i*2*hwindow/dim1,-hwindow+j*2*hwindow/dim2,-hwindow+k*2*hwindow/dim3);

                }
            }
        }   

    }

};


int main(){
    /* clock_t clock(void) returns the number of clock ticks
       elapsed since the program was launched.To get the number 
       of seconds used by the CPU, you will need to divide by 
       CLOCKS_PER_SEC.where CLOCKS_PER_SEC is 1000000 on typical
       32 bit system.  */
    clock_t start, end;

    cout << setprecision(15) << fixed;
    
    cout << "Calculation starts..." <<'\n';

    /* Recording the starting clock tick.*/
    start = clock();

    // MPI 
    // MPI_Init(NULL, NULL);

    // computation for Xarr, Yarr,Zarr
    

    // MPI_Finalize();
    // Recording the end clock tick.
    end = clock();
  
    // Calculating total time taken by the program.
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);

    cout << "Congratulations! You've done it." <<endl;

    cout << "Time taken by the calculation is : " << fixed 
         << time_taken << setprecision(6);
    cout << " sec " << endl;

    // Write those arrays into a txt file:
    
    ofstream myEx;
    ofstream myEy;
    ofstream myEz;

    myEx.open ("FocalField Radial Ex.txt");
    myEy.open ("FocalField Radial Ey.txt");
    myEz.open ("FocalField Radial Ez.txt");
    
    // myEx << "#Focal field of the radial polarization.\n";
    cout << "Writing into txt files." <<endl;
    for(int i=0; i < dim1; i++){
        for(int j=0; j < dim2; j++){
            for(int k=0; k < dim3; k++){
                
                // writing Xarr Yarr Zarr
                
                myEx << Xarr[i][j][k];
                myEx << "  "; // spacing between different complex numbers
                myEy << Yarr[i][j][k]; 
                myEy << "  ";
                myEz << Zarr[i][j][k]; 
                myEz << "  ";

            }
            
        }
        
    }

    myEx.close();
    myEy.close();
    myEz.close();
    cout << "Done!." <<endl;

    return 0;
}

