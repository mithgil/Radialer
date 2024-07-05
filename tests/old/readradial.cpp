#include <fstream>
#include <iostream>
using namespace std;
 
int main(){

    char data[1000];

    // open a file in read mode.
    ifstream infile; 
    infile.open("FocalField Radial Ex.txt"); 

    cout << "Reading from the file" << endl; 
    infile >> data; 

    // write the data at the screen.
    cout << data << endl;

    // again read the data from the file and display it.
    infile >> data; 
    cout << data << endl; 

    // close the opened file.
    infile.close();

    return 0;
}