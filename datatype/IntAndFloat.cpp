
#include <iostream>
#include <cmath>

int main(){

  /*
  // Let's find how many bits are in an unsigned integer, using "underflow":
  unsigned int A = 0;
  A = A - 1;
  std::cout << "unsigned 0 - 1: " << A << std::endl;
  std::cout << "Number of bits: " << log(float(A))/log(2.) << std::endl;
  */

  /*
  // Let's figure out how many bits are in a regular integer, using "overflow":
  int B = 1;
  for(int i=0; i<35; i++){
    std::cout << "i: " << i << ", 2^i: " << B << std::endl;
    B *= 2;
  }
  */
  //std::cout << sqrt(-1.) << std::endl;

  /*
  // Let's figure out how many bits are given to the exponent of a float, using overflow:
  float C = 1.;
  for(int i=0; i<130; i++){
    std::cout << "i: " << i << ", 2^i: " << C << std::endl;
    C = C * 2.;
  }
  */

  /*
  // Let's figure out how many bits are given to the exponent of a double, using overflow:
  double D = 1.;
  for(int i=0; i<1030; i++){
    std::cout << "i: " << i << ", 2^i: " << D << std::endl;
    D *= 2.;
  }
  */

  
  // Let's figure out how many bits are given to the mantissa of a double, using roundoff error:
  double E = 1.;
  for(int i=0; i<60; i++){
    std::cout << "i = " << i << ", Does 2^i + 1 = 2^i?    " << (E + 1. == E) << std::endl;
    E *= 2.;
  }
  

  
  return 0;
}
