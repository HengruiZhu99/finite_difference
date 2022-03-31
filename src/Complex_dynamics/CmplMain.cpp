
#include <iostream>
#include "MyComplex.hpp"


int main(){

  complex z(1.,2.);
  complex a(3.);

  std::cout << "z = " << z.real() << " + i*" << z.imag() << std::endl;
  std::cout << "a = " << a.real() << " + i*" << a.imag() << std::endl;

  a.set_imag(4.);
  std::cout << "a = " << a.real() << " + i*" << a.imag() << std::endl;

  std::cout << "a* = " << a.conj().real() << " + i*(" << a.conj().imag() << ")" << std::endl;

  std::cout << "|a| = " << a.abs() << std::endl;

  complex b = a*2.;
  std::cout << "b = " << b.real() << " + i*" << b.imag() << std::endl;

  b = 2.*a;
  std::cout << "b = " << b.real() << " + i*" << b.imag() << std::endl;

  complex d = 2.*a - z;
  std::cout << "d = " << d.real() << " + i*(" << d.imag() << ")" << std::endl;
  
  return 0;
}
