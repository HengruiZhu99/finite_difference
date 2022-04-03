
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cassert>

std::vector<double> operator*(double factor, std::vector<double> vect){
  for(int i=0; i<vect.size(); i++){
    vect[i] *= factor;
  }
  return vect;
}

std::vector<double> operator*(std::vector<double> vect, double factor){
  for(int i=0; i<vect.size(); i++){
    vect[i] *= factor;
  }
  return vect;
}


std::vector<double> operator*(int factor, std::vector<double> vect){
  for(int i=0; i<vect.size(); i++){
    vect[i] *= factor;
  }
  return vect;
}

std::vector<double> operator*(std::vector<double> vect1, std::vector<double> vect2){
  assert(vect1.size() == vect2.size());
  std::vector<double> vR(vect1.size());
  for(int i=0; i<vect1.size(); i++){
    vR[i] = vect1[i]*vect2[i];
  }
  return vR;
}

std::vector<double> operator+(std::vector<double> vect1, std::vector<double> vect2){
  assert(vect1.size() == vect2.size());
  std::vector<double> vR(vect1.size());
  for(int i=0; i<vR.size(); i++){
    vR[i] = vect1[i] + vect2[i];
  }
  return vR;
}

std::vector<double> operator/(std::vector<double> vect, double factor){
  for(int i=0; i<vect.size(); i++){
    vect[i] /= factor;
  }
  return vect;
}

std::vector<double> operator/(std::vector<double> vect, int factor){
  for(int i=0; i<vect.size(); i++){
    vect[i] /= factor;
  }
  return vect;
}

std::string datastring(double u){
  std::stringstream s;
  s << u;
  return s.str();
}

std::string datastring(std::vector<double> u){
  std::stringstream s;
  for(int i=0; i<u.size(); i++){
    s << u[i] << " ";
  }
  return s.str();
}

template <typename T>
void EulerMethod(double tmin, double tmax, double dt, T u0, T (*f)(T)){
  int Nsteps = int((tmax-tmin)/dt);
  T u = u0;

  std::ofstream outfile;
  outfile.open("EulerStepping.dat");

  for(int i=0; i<Nsteps; i++){
    double t = tmin + i*dt;
    u = u + dt*f(u);
    if(i%100==0){
      outfile << t << " " << datastring(u) << std::endl;
    }
  }

  outfile.close();  
}

template <typename T>
void MidpointMethod(double tmin, double tmax, double dt, T u0, T (*f)(T)){
  int Nsteps = int((tmax-tmin)/dt);
  T u = u0;

  std::ofstream outfile;
  outfile.open("MidpointStepping.dat");

  for(int i=0; i<Nsteps; i++){
    double t = tmin + i*dt;
    T k = u + .5*dt*f(u);
    u = u + dt*f(k);
    if(i%100==0){
      outfile << t << " " << datastring(u) << std::endl;
    }
  }

  outfile.close();
}


template <typename T>
void RK4(double tmin, double tmax, double dt, T u0, T (*f)(T)){
  int Nsteps = int((tmax-tmin)/dt);
  T u = u0;

  std::ofstream outfile;
  outfile.open("RK4.dat");

  outfile << std::setprecision(15);
  
  for(int i=0; i<Nsteps; i++){
    double t = tmin + i*dt;
    T k1 = f(u);
    T k2 = f(u+dt*k1*0.5);
    T k3 = f(u+dt*k2*0.5);
    T k4 = f(u+dt*k3);
    u = u + dt*(k1+(2*k2)*(2*k3)+k4)/6;
    if(i%100==0){
    outfile << t << " " << datastring(u) << std::endl;
    }
  }
}


double logistic(double u){
  return u*(1.-u);
}

std::vector<double> Hooke(std::vector<double> u){
  double x = u[0];
  double v = u[1];

  double xdot = v;
  double vdot = -x;

  return std::vector<double>{xdot, vdot};
}


int main(){

  double u0 = 0.1;
  //EulerMethod(0.,10.,1.e-3, u0, logistic);
  MidpointMethod(0.,10.,1.e-3, u0, logistic);
  RK4(0.,10.,1.e-3, u0, logistic);

  //std::vector<double> u0{0.,1.};
  //EulerMethod(0.,50.,1.e-3, u0, Hooke);
  //MidpointMethod(0.,50.,1.e-3, u0, Hooke);
  //RK4(0.,50.,1.e-3, u0, Hooke);

  return 0;
}
