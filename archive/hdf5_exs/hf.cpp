/* An example code for creating hdf5
 https://officeguide.cc/c-language-hdf5-program-compilation-tutorial/
 2023/11/28

 compile with:
 $ h5c++ hf.cpp -o hf
*/


#include <iostream>
#include <string>
#include <H5Cpp.h>
using namespace H5;

int main (void) {
  // 以預設屬性建立新的 HDF5 檔案
  H5File file("example.h5", H5F_ACC_TRUNC);

  // 建立資料空間
  hsize_t dims[2];
  dims[0] = 4;
  dims[1] = 6;
  DataSpace dataspace(2, dims);

  // 建立 HDF5 資料集
  DataSet dataset = file.createDataSet("my_dataset", PredType::STD_I32BE, dataspace);

  return 0;
}

