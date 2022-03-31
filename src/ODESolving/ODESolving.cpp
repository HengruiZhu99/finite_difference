
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cassert>

void EulerMethod(double tmin, double tmax, double dt, double u0, double (*f)(double)){
  int Nsteps = int((tmax-tmin)/dt);
  double u = u0;

  std::ofstream outfile;//ofstream allows us to output variable
  outfile.open("EulerStepping.dat");

  outfile << std::setprecision(15);
  
  for(int i=0; i<Nsteps; i++){
    double t = tmin + i*dt;
    u += dt*f(u);
    //if(i%100==0){
    outfile << t << " " << u << std::endl;
    //}
  }

  outfile.close();
  
}


std::vector<double> operator*(double factor, std::vector<double> vect){
  for(int i=0; i<vect.size(); i++){
    vect[i] *= factor;
  }
  return vect;
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





void EulerMethod(double tmin, double tmax, double dt, std::vector<double> u0, std::vector<double> (*f)(std::vector<double>)){
  int Nsteps = int((tmax-tmin)/dt);
  std::vector<double> u = u0;

  std::ofstream outfile;
  outfile.open("EulerStepping.dat");

  outfile << std::setprecision(15);
  
  for(int i=0; i<Nsteps; i++){
    double t = tmin + i*dt;
    u = u + dt*f(u);
    if(i%100==0){
      outfile << t << " " << u[0] << std::endl;
    }
  }

  outfile.close();  
}


void MidpointMethod(double tmin, double tmax, double dt, double u0, double (*f)(double)){
  int Nsteps = int((tmax-tmin)/dt);
  double u = u0;

  std::ofstream outfile;
  outfile.open("MidpointStepping.dat");

  outfile << std::setprecision(15);
  
  for(int i=0; i<Nsteps; i++){
    double t = tmin + i*dt;
    double k = u + .5*dt*f(u);
    u += dt*f(k);
    //if(i%100==0){
    outfile << t << " " << u << std::endl;
    //}
  }

  outfile.close();
}

void RK4(double tmin, double tmax, double dt, double u0, double (*f)(double)){
  int Nsteps = int((tmax-tmin)/dt);
  double u = u0;

  std::ofstream outfile;
  outfile.open("RK4.dat");

  outfile << std::setprecision(15);
  
  for(int i=0; i<Nsteps; i++){
    double t = tmin + i*dt;
    double k1 = f(u);
    double k2 = f(u+dt*k1*0.5);
    double k3 = f(u+dt*k2*0.5);
    double k4 = f(u+dt*k3);
    u += dt*(k1+2*k2*2*k3+k4)/6;
    //if(i%100==0){
    outfile << t << " " << u << std::endl;
    //}
  }

  outfile.close();
}

void MidpointMethod(double tmin, double tmax, double dt, std::vector<double> u0, std::vector<double> (*f)(std::vector<double>)){
  int Nsteps = int((tmax-tmin)/dt);
  std::vector<double> u = u0;

  std::ofstream outfile;
  outfile.open("MidpointStepping.dat");

  outfile << std::setprecision(15);
  
  for(int i=0; i<Nsteps; i++){
    double t = tmin + i*dt;
    std::vector<double> k = u + .5*dt*f(u);
    u = u + dt*f(k);
    if(i%100==0){
      outfile << t << " " << u[0] << std::endl;
    }
  }

  outfile.close();
}

void RK4(double tmin, double tmax, double dt, std::vector<double> u0, std::vector<double> (*f)(std::vector<double>)){
  int Nsteps = int((tmax-tmin)/dt);
  std::vector<double> u = u0;

  std::ofstream outfile;
  outfile.open("RK4.dat");

  outfile << std::setprecision(15);
  
  for(int i=0; i<Nsteps; i++){
    double t = tmin + i*dt;
    std::vector<double> k1 = f(u);
    std::vector<double> k2 = f(u+dt*0.5*k1);
    std::vector<double> k3 = f(u+dt*0.5*k2);
    std::vector<double> k4 = f(u+dt*k3);
    u = u + dt*(k1+(2*k2)*(2*k3)+k4)/6;
    if(i%100==0){
      outfile << t << " " << u[0] << std::endl;
    }
  }

  outfile.close();
}




double logistic(double u){
  return u*(1.-u);//exponential growth and levels off
}

std::vector<double> Hooke(std::vector<double> u){
  double x = u[0];
  double v = u[1];

  double xdot = v;
  double vdot = -x;

  return std::vector<double>{xdot, vdot};
}


int main(){

  //double u0 = 0.1;

  //EulerMethod(0.,10.,1.e-1, u0, logistic);
  //MidpointMethod(0.,10.,1.e-3, u0, logistic);
  //RK4(0.,10.,1.e-3, u0, logistic);

  std::vector<double> u0{0.,1.};
  EulerMethod(0.,50.,1.e-6, u0, Hooke);
  MidpointMethod(0.,50.,1.e-3, u0, Hooke);
  RK4(0.,50.,1.e-3, u0, Hooke);

  return 0;
}
