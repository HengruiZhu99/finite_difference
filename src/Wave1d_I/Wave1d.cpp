
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "DataGridMultiField.hpp"


std::vector<double> WaveFluxFunction(std::vector<double> u){
  return std::vector<double>{0.,-u[2],-u[1]};
}
std::vector<double> WaveSourceFunction(std::vector<double> u){
  // Here are source functions for the standard wave equation:
  return std::vector<double>{u[1], 0., 0.};

  // Alternatively, we could comment out the above to trade it
  // for this source function, which adds the simplest kind of
  // nonlinear "self-interaction" for which the PDE still has 
  // a solution over long spans of time.
  // It's been a while since I've looked at this, but
  // I believe the classical limit of the Higgs field has
  // this kind of self-interaction.
  //return std::vector<double>{u[1], -u[0]*u[0]*u[0], 0.};

  // Or, we could include nonlinear source terms analogous
  // to those in GR, which are quadratic in the first spatial
  // derivative of the field (the so-called "Christoffel
  // symbols."
  //return std::vector<double>{u[1], -u[0]*u[2]*u[2], 0.};

}

int main(){

  int N = 600;
  double x1 = 0.;
  double x2 = 10.;
  double xc = 5.;

  // Set initial data:
  // We'll let the wave field "psi" be a gaussian.
  // Its initial time derivative "pi" will be zero.
  // But if psi is nonzero and nonconstant in space,
  // then its derivative "phi" must be set to the
  // derivative of psi's gaussian.
  double w = .7;
  std::vector<double> psivalues(N, 0.);
  std::vector<double> phivalues(N, 0.);
  for(int i=0; i<N; i++){
    double x = x1 + (x2-x1)*i/(N-1);
    psivalues[i] = exp(-(x-xc)*(x-xc)/(w*w));
    phivalues[i] = -2.*((x-xc)/(w*w))*exp(-(x-xc)*(x-xc)/(w*w));
  }
  std::vector<std::vector<double> > uinit;
  for(int i=0; i<N; i++){
    uinit.push_back(std::vector<double>{psivalues[i],0.,phivalues[i]});
  }
  DataGrid u(x1, x2, N, 1, 3, uinit);

  double dx = (x2-x1)/(N-1);
  double dt = 1.*dx;

  std::ofstream outfile;
  outfile.open("wave.dat");

  double tfinal = 20.;
  int tindfinal = int(tfinal/dt);
    
  for(int tindex = 0; tindex<tindfinal; tindex++){
    // Time-stepping via the midpoint method:
    // DataGrid k = u + .5*dt*u.TimeDeriv(dx, WaveFluxFunction, WaveSourceFunction);
    // u  = u + dt*k.TimeDeriv(dx, WaveFluxFunction, WaveSourceFunction);

    // Time-stepping via "Runge-Kutta 4", a 4th-order accurate (in time) method that will allow
    // longer timesteps with less timestepper error.
    std::vector<std::vector<double> > k1 = u.TimeDeriv(dx, WaveFluxFunction, WaveSourceFunction);
    std::vector<std::vector<double> > k2 = (u+.5*dt*k1).TimeDeriv(dx, WaveFluxFunction, WaveSourceFunction);
    std::vector<std::vector<double> > k3 = (u+.5*dt*k2).TimeDeriv(dx, WaveFluxFunction, WaveSourceFunction);
    std::vector<std::vector<double> > k4 = (u+dt*k3).TimeDeriv(dx, WaveFluxFunction, WaveSourceFunction);
    u = u + (1./6.)*dt*(k1+2.*k2+2.*k3+k4);
    
    if(tindex%10==0){
      std::cout << "t = " << tindex*dt << " of " << tfinal << std::endl;
      for(int i=0; i<u.Ncells; i+=1){
    	outfile << u.Cell[i].values[0] << " ";
      }
      outfile << std::endl;
    }
  }
  

  return 0;
}
