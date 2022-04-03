#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "DataGrid.hpp"


double AdvectionRightFluxFunc(double x){
  return x;
}

inline double sqr(double x) { return x*x;}

int main(){

  int N = 1000;
  double x1 = 0.;
  double x2 = 10.;
  double xc = 5.;
  double w = .7;
  std::vector<double> values(N, 0.);
  for(int i=0; i<N; i++){
    double x = x1 + (x2-x1)*i/(N-1);
    values[i] = exp(-sqr(x-xc)/sqr(w));
  }

  DataGrid u(x1, x2, N, 1, values, 1);

  double dx = (x2-x1)/(N-1);
  double dt = .2*dx;

  std::ofstream outfile;
  outfile.open("advection.dat");

  double tfinal = 10.;
  int tindfinal = int(tfinal/dt);
    
  for(int tindex = 0; tindex<tindfinal; tindex++){
    DataGrid k = u + .5*dt*u.Dot(dx, AdvectionRightFluxFunc);
    u = u + dt*k.Dot(dx, AdvectionRightFluxFunc);
    if(tindex%50==0){
      std::cout << "t = " << tindex*dt << " of " << tfinal << std::endl;
      for(int i=0; i<u.Ncells; i+=1){
	outfile << u.Cell[i].value << " ";
      }
      outfile << std::endl;
    }
  }
  

  return 0;
}
