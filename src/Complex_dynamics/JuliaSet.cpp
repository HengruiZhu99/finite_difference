
#include <iostream>
#include <fstream>
#include "MyComplex.hpp"

complex QuadIterate(complex z){ return(z*z + complex(-.8,.156)); }

int main(){

  int Nx = 3200, Ny = 3200;
  double xmin = -2., xmax = 2.;
  double ymin = -2., ymax = 2.;
  //double xmin = -.001, xmax = .001;
  //double ymin = -.001, ymax = .001;
  double dx = (xmax-xmin)/Nx;
  double dy = (ymax-ymin)/Ny;
  double x, y;
  complex z;
  
  std::ofstream outfile;
  outfile.open("JuliaSet.dat");

  for(int j=0; j<Ny; j++){
    for(int i=0; i<Nx; i++){
      x = xmin+i*dx;
      y = ymax-j*dy;
      z = complex(x, y);
      int k=0;
      for(k=0; k<300; k++) {
        z = QuadIterate(z);
        if(z.abs()>1000.) {break;}
      }
      outfile << log(1.+k) << " ";
    }
    outfile << std::endl;
  }
  

  return 0.;
}
