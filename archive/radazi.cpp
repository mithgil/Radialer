/* 

This cpp script is for calculating 3D focal field of radial mode

    1. define numerical calculations for complex fields as functions
    2. define 4 dimentisonal array by 4-fold pointers to run calculation according to the 3d volume and other settings
    3. write those data into files
    4. delete memoryfor pointers

@ Tim

update: 2022/6/17

*/

#include <iostream> // cin cout
#include <fstream> // write files
#include <bits/stdc++.h> //
#include <math.h> // math: sqrt() pow()
#include <complex> 

using namespace std;

// total number of steps forthe following numerical integration 
const int NumStep = 250;

// Functions of Radial mode complex form of x,y,z components: (ex,ey,ez)

complex<double> ex_rad(float wavelength, float NA, float e0, float x, float y, float z){

    float k = 2*M_PI/wavelength;

    double f0 = 1.0; // ratio of pupil radius and beam waist, 1 for overfilling parabolic mirror/objectice

    float thetaM = asin(NA);

    // initialize two integrals
    double integRe = 0.0, integIm = 0.0;
    
    // For calculating numberical integration (mid-point method)
    double step = thetaM/ NumStep;  // width of each small rectangle

    for (int i = 0; i < NumStep; i ++) {

        // area += f((i + 0.5) * step) * step; // sum up each small rectangle
        // Apodization function: (pupil filter) exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))
        integRe += exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
        integIm += exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
    }
    // e0: amplitude
    integRe = e0*integRe*cos(atan2(y,x));
    integIm = e0*integIm*cos(atan2(y,x));
    

    return {integRe, integIm}; // To create an object you must use curly-braces {}
    
}

complex<double> ey_rad(float wavelength, float NA, float e0, float x, float y, float z){

    float k = 2*M_PI/wavelength;

    double f0 = 1.0;

    float thetaM = asin(NA);
    
    double integRe = 0.0, integIm = 0.0;
    
    double step = thetaM/ NumStep;  // width of each small rectangle

    for (int i = 0; i < NumStep; i ++) {

        integRe += exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
        integIm += exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
    }
    
    integRe = e0*integRe*sin(atan2(y,x));
    integIm = e0*integIm*sin(atan2(y,x));

    return {integRe, integIm};
    
}   

complex<double> ez_rad(float wavelength, float NA, float e0, float x, float y, float z){
    
    float k = 2*M_PI/wavelength;

    double f0 = 1.0;

    float thetaM = asin(NA);

    double integRe = 0.0, integIm = 0.0;
    
    double step = thetaM/ NumStep;  // width of each small rectangle

    for (int i = 0; i < NumStep; i ++) {

        integRe += exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*cyl_bessel_j(0, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
        integIm += exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*cyl_bessel_j(0, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
    }
    integRe = -2*e0*integRe;
    integIm = 2*e0*integIm;

    return {integRe, integIm};

}

// Functions of  azimuthal mode

/*
def Ex_a(x, y, z):
    intg_r = lambda s: -2*E0*fw(s)*(cos(s))**(0.5)*sin(s)*sp.jv(1,k*sqrt(x**2+y**2)*sin(s))*cos(k*z*cos(s))
    intg_i = lambda s: -2*E0*fw(s)*(cos(s))**(0.5)*sin(s)*sp.jv(1,k*sqrt(x**2+y**2)*sin(s))*sin(k*z*cos(s))
    return quad(intg_r, 0, thetaM)[0]*sin(arctan2(y,x))+1j*quad(intg_i, 0, thetaM)[0]*sin(arctan2(y,x))
    
def Ey_a(x, y, z):
    intg_r = lambda s: 2*E0*fw(s)*(cos(s))**(0.5)*sin(s)*sp.jv(1,k*sqrt(x**2+y**2)*sin(s))*cos(k*z*cos(s))
    intg_i = lambda s: 2*E0*fw(s)*(cos(s))**(0.5)*sin(s)*sp.jv(1,k*sqrt(x**2+y**2)*sin(s))*sin(k*z*cos(s))
    return quad(intg_r, 0, thetaM)[0]*cos(arctan2(y,x))+1j*quad(intg_i, 0, thetaM)[0]*cos(arctan2(y,x))

*/

complex<double> ex_azi(float wavelength, float NA, float e0, float x, float y, float z){
    
    float k = 2*M_PI/wavelength;

    double f0 = 1.0;

    float thetaM = asin(NA);

    double integRe = 0.0, integIm = 0.0;
    
    double step = thetaM/ NumStep;  // width of each small rectangle

    for (int i = 0; i < NumStep; i ++) {

        integRe += exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*sin((i + 0.5) * step)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
        integIm += exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*sin((i + 0.5) * step)*cyl_bessel_j(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
    }
    integRe = -2*e0*integRe*sin(atan2(y,x));
    integIm = 2*e0*integIm*sin(atan2(y,x));

    return {integRe, integIm};

}


complex<double> ey_azi(float wavelength, float NA, float e0, float x, float y, float z){
    
    float k = 2*M_PI/wavelength;

    double f0 = 1.0;

    float thetaM = asin(NA);

    double integRe = 0.0, integIm = 0.0;
    
    double step = thetaM/ NumStep;  // width of each small rectangle

    for (int i = 0; i < NumStep; i ++) {

        integRe += exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*cyl_bessel_j(0, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
        integIm += exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*cyl_bessel_j(0, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
    }
    integRe = -2*e0*integRe*cos(atan2(y,x));
    integIm = 2*e0*integIm*cos(atan2(y,x));

    return {integRe, integIm};

}

//  to create a 4d array (3 components)
complex<double> ****CreateArray(string ModeCall, int dim, float wavelengthIn, float amp0, float NAIn, float hwindow){
    
    // define size of array
    int width = dim;
    int length = dim;
    int height = dim;
    // create a 4d array
    complex<double> ****arr = new complex<double> ***[3];// 3 components ex,ey,ez

    for(int l = 0; l < 3; ++l){
        arr[l] = new complex<double> **[width];
        for(int i = 0; i < width; ++i){
            arr[l][i] = new complex<double> *[length];
            for(int j = 0; j < length; ++j){
                for(int j = 0; j < length; ++j){
                    arr[l][i][j] = new complex<double>[height];
                    
                }
            }
        }
    }

    cout << "Calculation for " ;
    
    clock_t start, end;
    start = clock();

    cout << setprecision(10) << fixed;
    
            cout << "ex ey & ez ..."<<'\n';
            for(int i=0; i < dim; i++){
                for(int j=0; j < dim; j++){
                    for(int k=0; k < dim; k++){
                        //  return ex, ey, ez as a 4d array
                        if (ModeCall == "Radial"){
                            
                            arr[0][i][j][k] = ex_rad(wavelengthIn,NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim);
                            arr[1][i][j][k] = ey_rad(wavelengthIn,NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim);
                            arr[2][i][j][k] = ez_rad(wavelengthIn,NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim);
                        } else if (ModeCall == "Azimuthal") {

                            arr[0][i][j][k] = ex_azi(wavelengthIn,NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim);
                            arr[1][i][j][k] = ey_azi(wavelengthIn,NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim);
                            arr[2][i][j][k] = 0;

                        } else {
                            cout <<"Mode name is wrong, please select correct one!"<<endl;
                        }
                        
                    }
                }
            }
 
    
    end = clock();

    cout << " is finished!" <<endl;
    /* clock_t clock(void) returns the number of clock ticks
       elapsed since the program was launched.To get the number 
       of seconds used by the CPU, you will need to divide by 
       CLOCKS_PER_SEC.where CLOCKS_PER_SEC is 1000000 on typical
       32 bit system.  */

    // calculate the time for the calculation
    double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
    cout << " Time taken by the calculation is : " << fixed 
            << time_taken << setprecision(6);
    cout << " sec " << endl;

    return arr;

}

int main(){
    

    /* Optical Setup:
     int comp, int dim, float wavelengthIn, float amp0, float NAIn, float hwindow */

    int dimSet = 81;
    float wavelengthSet = 0.636; // um
    float ampSet = 1.0;
    float NASet = 0.998;
    float hwindowSet = 1.0;
    
    // choosing a mode with a user input
    string ModeName = "Radial";

    
    cout<<"3 components (ex,ey,ez) will be calculated"<<endl;
    // Radial field of 3 components by using a 4-fold pointer
    complex<double>**** CalArr = CreateArray(ModeName, dimSet, wavelengthSet,ampSet, NASet, hwindowSet);
    /*
    } else {
        // azimuthal mode 
        cout<<"2 components (ex,ey,with ez=0) will be calculated"<<endl;
        complex<double>**** CalArr = CreateArray(ModeName, dimSet, wavelengthSet,ampSet, NASet, hwindowSet);

    }
    */
    
    // Write those arrays into a txt file:
       

    if (ModeName == "Radial"){
        ofstream myEz;
        myEz.open ("FocalField Radial Ez.txt");
    }

    ofstream myEx;
    ofstream myEy;

    myEx.open ("FocalField Radial Ex.txt");
    myEy.open ("FocalField Radial Ey.txt");
    
  
    
    cout << "Writing into txt files ... ";
    // int dimIn = ;??

    for(int i=0; i < dimSet; i++){
        for(int j=0; j < dimSet; j++){
            for(int k=0; k < dimSet; k++){
                
                // writing Xarr Yarr Zarr
                
                myEx << CalArr[0][i][j][k];
                myEx << "  "; // spacing between different complex numbers
                myEy << CalArr[1][i][j][k]; 
                myEy << "  ";
                if (ModeName == "Radial"){
                    myEz << CalArr[2][i][j][k]; 
                    myEz << "  ";
                }
            }
            
        }
        
    }

    myEx.close();
    myEy.close();
    if (ModeName == "Radial"){
        myEz.close();
    }

    cout << "..done!." <<endl;
    

    for (int l = 0; l < 3; ++l){

        for(int i = 0; i < dimSet; ++i){
        
            for( int j = 0; j<dimSet;j++){
                delete[] CalArr[l][i][j];
                
            }
            delete[] CalArr[l][i];
        
        }
        delete[] CalArr[l];
    }
    
    delete[] CalArr;
    

    CalArr = NULL;
    
    cout << "Memory for storing arrays is deleted." <<endl;

    return 0;
}
