/* 
This cpp script is for calculating 3D focal field of radial mode

    1. define numerical calculations for complex fields as functions
    2. define 4 dimentisonal array by 4-fold pointers to run calculation according to the 3d volume and other settings
    3. write those data into files
    4. delete memoryfor pointers

Note: single quote 'x' for char and double quote "x" for string (a char and \n)

Export data into hdf5 format


@ Tim

update: 2023/10/2
*/


#include <iostream> // cin cout-must have!
#include <fstream> // write files
#include <math.h> // math: sqrt() pow()
#include <complex> 
#include <iomanip>      // std::setprecision
#include "radial.h" // header of functions

// #include <bits/stdc++.h> // non typical lib - all above except for iostream

using namespace std; // requires iostream



int main(){
    /* clock_t clock(void) returns the number of clock ticks
       elapsed since the program was launched.To get the number 
       of seconds used by the CPU, you will need to divide by 
       CLOCKS_PER_SEC.where CLOCKS_PER_SEC is 1000000 on typical
       32 bit system.  */
    
    // cout << setprecision(15) << fixed;

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

    //the allocated memory must be deleted:

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