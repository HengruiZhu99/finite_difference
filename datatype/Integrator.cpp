
#include <iostream>
#include <cmath>


double Square(double x){
  return x*x;
}


double RectangleRule(double xmin, double xmax, int N){
  double dx = (xmax-xmin)/N;
  double Sum = 0.;
  for(int i=0; i<N; i++){
    double x = xmin + i*dx;
    Sum += dx*Square(x);
  }
  return Sum;
}

double RectangleRule(double xmin, double xmax, int N, double (*f)(double) ){
  double dx = (xmax-xmin)/N;
  double Sum = 0.;
  for(int i=0; i<N; i++){
    double x = xmin + i*dx;
    Sum += dx*f(x);
  }
  return Sum;
}

double TrapezoidalRule(double xmin, double xmax, int N){
  double dx = (xmax-xmin)/N;
  double Sum = 0.;
  for(int i=0; i<N; i++){
    double x = xmin + i*dx;
    Sum += dx*.5*(Square(x)+Square(x+dx));
  }
  return Sum;
}

double TrapezoidalRule(double xmin, double xmax, int N, double (*f)(double) ){
  double dx = (xmax-xmin)/N;
  double Sum = 0.;
  for(int i=0; i<N; i++){
    double x = xmin + i*dx;
    Sum += dx*.5*(f(x)+f(x+dx));
  }
  return Sum;
}



int main(){

  int N = int(1.e9);
  double IRect = RectangleRule(0., 1., N);
  double ITrap = TrapezoidalRule(0., 1., N);
  std::cout << "Rectangle integral: " << IRect << std::endl;
  std::cout << "             error: " << IRect - 1./3. << std::endl;
  std::cout << "Trapezoidal integral: " << ITrap << std::endl;
  std::cout << "               error: " << ITrap - 1./3. << std::endl;
  
  return 0;
}