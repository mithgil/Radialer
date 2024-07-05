/* 
This cpp script is for calculating 3D focal field of radial mode

    1. define numerical calculations for complex fields as functions
    2. define 4 dimentisonal array by 4-fold pointers to run calculation according to the 3d volume and other settings
    3. write those data into files
    4. delete memoryfor pointers

@ Tim
update: 2023/2/10

*/


#include <iostream> // cin cout
#include <fstream> // write files
#include <bits/stdc++.h> //
#include <math.h> // math: sqrt() pow()
#include <complex> 

using namespace std;

//factorial
long fact(int n){

	long factorial = 1.0;
	for(int i=1; i<=n ; ++i){
		factorial*=i;
	}

	return factorial;
}

// nth order besselJ at x position with nonzero integer n & n>=0
double besselj(int n, double x){
	double sum = 0;
	int M=12; // note this will make fact overflow

	for (int k=0; k<=M; ++k){
		//cout<<fact(k)<<endl;
		sum += pow(-1,k)*pow(0.5*x,2*k+n)/fact(k)/fact(k+n);
	}

	return sum;
}

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
        integRe += e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*besselj(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
        integIm += e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*besselj(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
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

        integRe += e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*besselj(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
        integIm += e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*pow(cos((i + 0.5) * step),1.5)*pow(sin(2*(i + 0.5) * step),2)*besselj(1, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
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

        integRe += -2*e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*besselj(0, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*sin(k*z*cos((i + 0.5) * step))*step;
        integIm += 2*e0*exp(-1*pow(f0,-2)*pow(sin((i + 0.5) * step)/sin(thetaM),2))*sqrt(cos((i + 0.5) * step))*pow(sin((i + 0.5) * step),2)*besselj(0, k*sqrt(x*x+y*y)*sin((i + 0.5) * step))*cos(k*z*cos((i + 0.5) * step))*step;
    }
    
    return {integRe, integIm};

}

//  to create a 4d array (3 components)
complex<double>**** CreateArray(int dim, float wavelengthIn, float amp0, float NAIn, float hwindow){

    int cube = 3;
    int width = dim;
    int length = dim;
    int height = dim;
    // create a 4d array
    complex<double>**** arr = new complex<double>***[cube];

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

    cout << "Calculation for " ;
    
    clock_t start, end;
    start = clock();

    cout << setprecision(10) << fixed;
    
            cout << "ex ey & ez ..."<<'\n';
            for(int i=0; i < dim; i++){
                for(int j=0; j < dim; j++){
                    for(int k=0; k < dim; k++){
                        //  return ex, ey, ez as a 4d array
                        arr[0][i][j][k] = ex(wavelengthIn,NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim);
                        arr[1][i][j][k] = ey(wavelengthIn,NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim);
                        arr[2][i][j][k] = ez(wavelengthIn,NAIn, amp0, -hwindow+i*2*hwindow/dim,-hwindow+j*2*hwindow/dim,-hwindow+k*2*hwindow/dim);
                    }
                }
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

    // Radial field of 3 components by using a 4-fold pointer
    complex<double>**** RadArr = CreateArray(dimSet, wavelengthSet,ampSet, NASet, hwindowSet);
    

    // Write those arrays into a txt file:
    
    ofstream myEx;
    ofstream myEy;
    ofstream myEz;

    myEx.open ("FocalField Radial Ex-1.txt");
    myEy.open ("FocalField Radial Ey-1.txt");
    myEz.open ("FocalField Radial Ez-1.txt");
    
    // myEx << "#Focal field of the radial polarization.\n";
    cout << "Writing into txt files ... ";
    for(int i=0; i < dimSet; i++){
        for(int j=0; j < dimSet; j++){
            for(int k=0; k < dimSet; k++){
                
                // writing Xarr Yarr Zarr
                
                myEx << RadArr[0][i][j][k];
                myEx << "  "; // spacing between different complex numbers
                myEy << RadArr[1][i][j][k]; 
                myEy << "  ";
                myEz << RadArr[2][i][j][k]; 
                myEz << "  ";

            }
            
        }
        
    }

    myEx.close();
    myEy.close();
    myEz.close();
    cout << "done!." <<endl;

    //allocated memory must be deleted:

    for (int l = 0; l < 3; ++l){

        for(int i = 0; i < dimSet; ++i){
        
            for( int j = 0; j<dimSet;j++){
                delete[] RadArr[l][i][j];
                
            }
            delete[] RadArr[l][i];
        
        }
        delete[] RadArr[l];
    }
    
    delete[] RadArr;
    

    RadArr = NULL;
    
    cout << "Memory for storing arrays is deleted." <<endl;

    return 0;
}

