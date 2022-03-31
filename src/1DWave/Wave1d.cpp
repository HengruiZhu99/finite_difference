
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "DataGridMultiField.hpp"


std::vector<double> WaveFluxFunction(std::vector<double> u){
  return std::vector<double>{0.,u[1],-u[2]};
}
std::vector<double> WaveSourceFunction(std::vector<double> u){
  return std::vector<double>{.5*(u[1]+u[2]), 0., 0.};
}

int main(){

  int N = 1000;
  double x1 = 0.;
  double x2 = 10.;
  double xc = 5.;
  double w = .7;
  std::vector<double> values(N, 0.);
  for(int i=0; i<N; i++){
    double x = x1 + (x2-x1)*i/(N-1);
    values[i] = exp(-(x-xc)*(x-xc)/(w*w));
    //values[i] = 5.*(x-xc)*exp(-(x-xc)*(x-xc)/(w*w));
  }
  std::vector<std::vector<double> > uinit;
  for(int i=0; i<N; i++){
    uinit.push_back(std::vector<double>{0.,values[i],values[i]});
  }
  DataGrid u(x1, x2, N, 1, 3, uinit, std::vector<double>{0.,1.,-1.});

  std::cout << "u.Nfields: " << u.Nfields << std::endl;
  std::cout << "u.Cell[0].Nfields: " << u.Cell[0].Nfields << std::endl;
  
  // for(int alpha=0; alpha<3; alpha++){
  //   std::cout << u.Cell[500].values[alpha] << std::endl;
  // }
  
  // Let's initialize the fields such that the basic wave field Psi
  // has initial value zero at all grid points:
  // DataGrid Psi(x1, x2, N, 1, std::vector<double>(N), 0); 
  // The above means that the "phi" field (dPsi/dx) has initial value
  // of zero as well, so u+ = u-, and I may as well set them both to
  // the gaussian just defined above.
  // DataGrid uplus(x1, x2, N, 1, values, 1);
  // DataGrid uminus(x1, x2, N, 1, values, -1);
  // Note that in the above, the last argument of these constructors
  // gives the field's propagation direction, + for uplus, - for uminus,
  // and zero for Psi itself, as it doesn't even have an advection term.
  
  double dx = (x2-x1)/(N-1);
  double dt = .5*dx;

  std::ofstream outfile;
  outfile.open("wave.dat");

  double tfinal = 10.;
  int tindfinal = int(tfinal/dt);
    
  for(int tindex = 0; tindex<tindfinal; tindex++){
    //DataGrid kplus = uplus + .5*dt*uplus.Dot(dx, AdvectionRightFluxFunc);
    //DataGrid kminus = uminus + .5*dt*uminus.Dot(dx, AdvectionLeftFluxFunc);
    DataGrid k = u + .5*dt*u.Dot(dx, WaveFluxFunction, WaveSourceFunction);
    u  = u + dt*k.Dot(dx, WaveFluxFunction, WaveSourceFunction);
    
    if(tindex%20==0){
      std::cout << "t = " << tindex*dt << " of " << tfinal << std::endl;
      for(int i=0; i<u.Ncells; i+=1){
    	outfile << u.Cell[i].values[0] << " ";
      }
      outfile << std::endl;
    }
  }
  

  return 0;
}
