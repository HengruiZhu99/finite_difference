
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include "DataGrid2d.hpp"
#include <omp.h>


std::vector<std::vector<double> > WaveFluxFunction(std::vector<double> Fields){
  std::vector<std::vector<double> > result;
  for(int i=0; i<2; i++){
    result.push_back(std::vector<double>(4));
  }
  result[0][0] = 0.;
  result[1][0] = 0.;
  result[0][1] = -Fields[2];
  result[1][1] = -Fields[3];
  result[0][2] = -Fields[1];
  result[1][2] = 0.;
  result[0][3] = 0.;
  result[1][3] = -Fields[1];
    
  return result;
}

std::vector<double> WaveSourceFunction(std::vector<double> Fields){
  std::vector<double> result(4);
  result[0] = Fields[1];
  result[1] = 0.-0.*Fields[0]*Fields[0]*Fields[0];
  result[2] = 0.;
  result[3] = 0.;
  
  return result;
}



int main(){

  
  int Nx = 100;
  int Ny = 100;
  int N = Nx*Ny;
  double x1 = 0.;
  double x2 = 14.;
  double xc = 2.5;
  double y1 = 0.;
  double y2 = 10.;
  double yc = 3.5;

  // Set initial data:
  // We'll let the wave field "psi" be a gaussian.
  // Its initial time derivative "pi" will be zero.
  // But if psi is nonzero and nonconstant in space,
  // then its gradient "phi" must be set to the
  // derivatives of psi's gaussian.
  double wx = .7;
  double wy = .7;
  std::vector<double> psivalues(N, 0.);
  std::vector<double> phixvalues(N, 0.);
  std::vector<double> phiyvalues(N, 0.);
  int I=0;
  for(int i=0; i<Nx; i++){
    for(int j=0; j<Ny; j++){
      double x = x1 + (x2-x1)*i/(Nx-1);
      double y = y1 + (y2-y1)*j/(Ny-1);
      double amp = 1.;
      I = i + Nx*j;
      psivalues[I] = amp*exp(-(x-xc)*(x-xc)/(wx*wx) - (y-yc)*(y-yc)/(wy*wy));
      phixvalues[I] = -amp*2.*((x-xc)/(wx*wx))*exp(-(x-xc)*(x-xc)/(wx*wx) - (y-yc)*(y-yc)/(wy*wy));
      phiyvalues[I] = -amp*2.*((y-yc)/(wy*wy))*exp(-(x-xc)*(x-xc)/(wx*wx) - (y-yc)*(y-yc)/(wy*wy));
    }
  }
  std::vector<std::vector<double> > uinit;
  for(int i=0; i<N; i++){
    uinit.push_back(std::vector<double>{psivalues[i],0.,phixvalues[i],phiyvalues[i]});
  }
  DataGrid u(x1, x2, y1, y2, Nx, Ny, 1, 4, uinit);


  // The following commented-out sections are not meant to be part of the code. Ordinarily I'd delete them.
  // But these code snippets were key to my figuring out what was going wrong with my indexing, so I may uncomment
  // some of it in today's discussion to show how it worked. 
  
  //double dx = (x2-x1)/(Nx-1);
  //double dy = (y2-y1)/(Ny-1);
  // int ind = 3;
  // for(int i=0; i<Nx; i++){
  //   for(int j=0; j<Ny; j++){
  //     I = i+Nx*j;
  //     std::cout << "(" << u.Cell[I].outDirection[ind][0] << ", " << u.Cell[I].outDirection[ind][1] << ") ";
  //   }
  //   std::cout << std::endl;
  // }

  // std::cout << std::endl;
  
  // for(int i=0; i<Nx; i++){
  //   for(int j=0; j<Ny; j++){
  //     I = i+Nx*j;
  //     std::cout << "(" << (u.Cell[u.Cell[I].neighbors[ind]].centerx - u.Cell[I].centerx)/dx << ", " << (u.Cell[u.Cell[I].neighbors[ind]].centery - u.Cell[I].centery)/dy << ") ";
  //   }
  //   std::cout << std::endl;
  // }

  // std::cout << std::endl;
  
  // for(int i=0; i<Nx; i++){
  //   for(int j=0; j<Ny; j++){
  //     I = i+Nx*j;
  //     std::cout << "(" << u.Cell[I].centerx << ", " << u.Cell[I].centery << ") ";
  //   }
  //   std::cout << std::endl;
  // }
  
  // std::cout << std::endl;
  
  // for(int i=0; i<Nx; i++){
  //   for(int j=0; j<Ny; j++){
  //     I = i+Nx*j;
  //     std::cout << "< " << I << ", " << u.Cell[I].neighbors[0] << " >  ";
  //   }
  //   std::cout << std::endl;
  // }

  //std::cout << std::endl;

  // int II = 20;
  // for(int nei=0; nei<4; nei++){
  //   std::cout << "I = " << II << ", index = " << u.Cell[II].index << ", neighbors[0] = " << u.Cell[II].neighbors[nei] << std::endl;
  // }

  
  double dx = (x2-x1)/(Nx-1);
  double dy = (y2-y1)/(Ny-1);
  double Courant = .9;
  double dt = Courant*dx;
  if(dy<dx){dt = Courant*dy;}
  
  std::ofstream outfile;
  outfile.open("wave.dat");

  double tfinal = 60.;
  int tindfinal = int(tfinal/dt);
  
  for(int tindex = 0; tindex<tindfinal; tindex++){
    // Time-stepping via the midpoint method:
    DataGrid k = u + .5*dt*u.TimeDeriv(dx, dy, WaveFluxFunction, WaveSourceFunction);
    u  = u + dt*k.TimeDeriv(dx, dy, WaveFluxFunction, WaveSourceFunction);

    // Time-stepping via "Runge-Kutta 4", a 4th-order accurate (in time) method that will allow
    // longer timesteps with less timestepper error. For today I'm commenting it out in favor of the
    // midpoint method above, as that method (as implemented here) is a little faster.
    //std::vector<std::vector<double> > k1 = u.TimeDeriv(dx, dy, WaveFluxFunction, WaveSourceFunction);
    //std::vector<std::vector<double> > k2 = (u+.5*dt*k1).TimeDeriv(dx, dy, WaveFluxFunction, WaveSourceFunction);
    //std::vector<std::vector<double> > k3 = (u+.5*dt*k2).TimeDeriv(dx, dy, WaveFluxFunction, WaveSourceFunction);
    //std::vector<std::vector<double> > k4 = (u+dt*k3).TimeDeriv(dx, dy, WaveFluxFunction, WaveSourceFunction);
    //u = u + (1./6.)*dt*(k1+2.*k2+2.*k3+k4);
    
    if(tindex%2==0){
      std::cout << "t = " << tindex*dt << " of " << tfinal << std::endl;
      for(int i=0; i<u.Ncells; i+=1){
    	outfile << u.Cell[i].values[0] << " ";
      }
      outfile << std::endl;
    }
  }


  return 0;
}
