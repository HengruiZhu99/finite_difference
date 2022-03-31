
#include <iostream>
#include <vector>
#include <cmath>

double sqr_by_value(double);
double sqr_by_reference(double&);

std::vector<double> cos_by_ref(std::vector<double>&);
std::vector<double> cos_by_val(std::vector<double>);

int main(){

  // A demo of references and dereferences with basic types (here, double):
  std::cout << "Floating-point number:" << std::endl;
  double x = 3.14159;
  std::cout << "x = " << x << "   (value of x)" << std::endl;
  std::cout << "&x = " << &x << "    (address of [aka, \'reference to\'] x)" << std::endl;
  std::cout << "*(&x) = " << *(&x) << "    (value [aka \'de-reference\'] of address [reference] of x)" << std::endl;

  /*
  // A demo of references and dereferences with std::vector:
  std::cout << std::endl;
  std::cout << "vector<int>:" << std::endl;
  std::vector<int> V{3,1,4,1,5,9};
  std::cout << "Address of vector V: " << &V << std::endl;
  std::cout << "Address of element 0 of V: " << &V[0] << std::endl;
  std::cout << "Address of element 1 of V: " << &V[1] << std::endl;
  std::cout << "Address of element 2 of V: " << &V[2] << std::endl;
  std::cout << "Value of element 0 of V: " << *&V[0] << std::endl;
  std::cout << "Value of element 1 of V: " << *&V[1] << std::endl;
  std::cout << "Value of element 2 of V: " << *&V[2] << std::endl;
  */

  /*
  // A demo of the kind of way we dealt with arrays before std::vector existed:
  std::cout << std::endl;
  std::cout << "raw memory allocated with \'new\': " << std::endl;
  double *W = new double [6];
  std::cout << "W = " << W << std::endl;
  std::cout << "&W[0] = " << &(W[0]) << std::endl;
  std::cout << "W[0] = " << W[0] << std::endl;
  std::cout << "*W = " << *W << std::endl;
  std::cout << "&W[1] = " << &W[1] << " (address of W[1])" << std::endl;
  delete[] W;
  */

  /*
  // A demo of how you can use references in functions, and the main risk of the process:
  std::cout << std::endl;
  std::cout << "Passing to a function \'by value\' versus \'by reference\'" << std::endl;
  std::cout << "Declaring a double: y=5" << std::endl;
  double y = 5;
  std::cout << "  Computing the square using a standard function (passing by value): sqr_by_value(y) = " << sqr_by_value(y) << std::endl;
  std::cout << "     And just to double-check, the value of y is: y = " << y << std::endl;
  std::cout << std::endl;
  std::cout << "  Computing the square using a function that passes by reference: sqr_by_reference(y) = " << sqr_by_reference(y) << std::endl;
  std::cout << "     But now let\'s check what y is: y = " << y << "   (uh oh!)" << std::endl;
  */  

  /*
  // A demo of one reason why you might want to pass references to functions anyway:
  std::cout << std::endl;
  std::cout << "Testing a function on a large vector:" << std::endl;
  int N = int(1.e8);
  std::vector<double> W(N, 2.);
  std::cout << "Testing pass-by-reference or pass-by-value function..." << std::endl;
  std::cout << "   Check memory usage (for instance, using 'top' in linux/unix) to see effect." << std::endl;
  for(int i=0; i<10; i++){
    std::cout << "Iteration " << i << std::endl;
    W = cos_by_ref(W);
    //W = cos_by_val(W);
  }
  std::cout << "W[0] = " << W[0] << std::endl;
  */
  
  return 0;
}

double sqr_by_value(double x){
  x *= x;
  return x;
}

double sqr_by_reference(double& x){
  x *= x;
  return x;
}


std::vector<double> cos_by_ref(std::vector<double>& V){
  for(int j=0; j<V.size(); j++){
    V[j] = cos(V[j]);
  }
  return V;
}

std::vector<double> cos_by_val(std::vector<double> W){
  for(int j=0; j<W.size(); j++){
    W[j] = cos(W[j]);
  }
  return W;
}

