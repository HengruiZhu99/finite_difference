
#pragma once

#include <cmath>

class complex{
  // Define the "member variables" (in this case, real and imaginary parts), which
  // by tradition are kept "private" meaning that they can't be accessed except
  // via member functions.
private:
  double re, im;

  // Now, let's start declaring member functions. Member functions are typically
  // defined to be public.
public:
  // A "constructor", essentially instructions to create a complex:
  complex(double r=0., double i=0.){re = r; im = i;}

  // Some "getter" functions, simply returning the values of the member
  // variables:
  double real(){ return re; }
  double imag(){ return im; }

  // Some "setter" functions, that let you modify the member variables:
  void set_real(double r){ re = r; }
  void set_imag(double i){ im = i; }
  // (You might note that if you have "getter" and "setter" functions for
  //  each variable, then you might as well have just made the member variables
  //  public to start with. This is true. The culture of programming paradigms
  //  is silly sometimes.)

  // Now, how about some functions for complex algebra stuff:
  complex conj(){ return complex(re, -im);}

  double abs(){ return sqrt(re*re + im*im); }

  // How about defining some arithmetic operations... Let's say we might want
  // to multiply a complex by a double:
  complex operator*(double s){ return complex(s*re, s*im); }
  
};


// Let's define arithmetic operations outside the class declaration:
complex operator*(double s, complex z) { return complex(s*z.real(), s*z.imag());}

complex operator*(complex u, complex v){
  double realpart = u.real()*v.real() - u.imag()*v.imag();
  double imagpart = u.real()*v.imag() + u.imag()*v.real();
  return complex(realpart, imagpart);
}

complex operator+(complex u, complex v){ return complex(u.real()+v.real(), u.imag()+v.imag()); }

complex operator+(double r, complex z){ return complex(z.real()+r, z.imag()); }
complex operator+(complex z, double r){ return complex(z.real()+r, z.imag()); }

//The minus operator to "negate" a single complex number:
complex operator-(complex z){ return complex(-z.real(), -z.imag()); }
//The minux operator to subtract a second complex number from a first one:
complex operator-(complex u, complex v){return complex(u.real()-v.real(), u.imag()-v.imag());}
