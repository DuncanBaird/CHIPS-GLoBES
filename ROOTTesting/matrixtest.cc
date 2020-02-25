#include <fstream>
#include <iostream>
#include <cmath>
#include <stdlib.h>

using namespace std;

#include <TMatrixD.h>

int setIdentitiy(TMatrixD &dummy1,int len){
  for(int i = 0;i<len;i++){
    for(int j = 0;j<len;j++){
        dummy1[i][j] = 0;
      }
    }
  for(int i = 0;i<len;i++){

    dummy1[i][i]=1;

  }
  return 0;
}

int setValue(TMatrixD &dummy1,int len_r,int len_c, int value){
  for(int i = 0;i<len_r;i++){
    for(int j = 0;j<len_c;j++){
        dummy1[i][j] = value;
      }
    }

  return 0;
}

void matrixtest(){

  const int len = 2;

  TMatrixD test1(len,len);
  TMatrixD test2(len,len);

  TMatrixD vertical(len,1);
  TMatrixD horizontal(1,len);

  // for(int i = 0;i<len;i++){
  //   for(int j = 0;j<len;j++){
  //       test1[i][j] = 0;
  //     }
  //   }
  
  // for(int i = 0;i<len;i++){

  // test1[i][i]=1;

  // }

  setIdentitiy(test1,len);

  setValue(test2,len,len,4);

  setValue(horizontal,1,len,4);
  setValue(vertical,len,1,4);
  

  // test1.Print();
  test2.Print();

  // TMatrixD test3 = test2 * test2;

  // test3.Print();

  vertical.Print();
  horizontal.Print();

  TMatrixD decomp = horizontal * test2 * vertical;

  decomp.Print();
}