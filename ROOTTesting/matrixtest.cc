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

  const int len = 4;

  TMatrixD test1(len,len);
  TMatrixD test2(len,len);

  TMatrixD vertical(len,1);
  TMatrixD vertical2(len,1);
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
  // test2.Print();

  // TMatrixD test3 = test2 * test2;

  // test3.Print();

  // vertical.Print();
  // horizontal.Print();

  vertical2.T().Print();

  test2[0][0] = 12.0;
  test2.Print();


  double det1;

  test2.Invert(&det1);
  test2.Print();

  
  TMatrixD decomp = horizontal * test2 * vertical;

  decomp.Print();

  double x = decomp[0][0];

  cout << "chi is " << x << "\n";

  TMatrixD delta1(200,1);
  TMatrixD delta2(200,1);

  TMatrixD covariance(200,200);
  setValue(covariance,200,200,6);
  covariance[0][0] = 100.0;

  for(int i=0;i<200;i++){
    delta1[i][0] = 1;
    delta2[i][0] = 1;
  }
  delta1.T();
  (delta1*covariance.Invert(&det1)*delta2).Print();
}