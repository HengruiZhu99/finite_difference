
#pragma once

#include <vector>


// Class for data in a single finite-volume cell:
class cell{
 public:
  // Constructor prototype (I'll actually define the constructor in DataGrid.cpp):
  cell(int Ind, double Cent_x, double Cent_y, int Nfields, std::vector<double> Vals, std::vector<int> Neighbors, std::vector<std::vector<double> > OutDir, std::vector<bool> IsInterior, std::vector<double> BdryLocationx, std::vector<double> BdryLocationy);

  // Member variables:
  //private: // (for simplicity, I'll break protocol and let all member variables be public)
  int index; // A label for the cell in a grid.
  int Nfields; // How many fields are being evolved in this (or any) cell?
  double centerx; // Coordinate of the cell center.
  double centery; // Coordinate of the cell center.
  std::vector<double> values;  // Cell's average field value.
  std::vector<int> neighbors; // Cells that neighbor this one in a grid (2 sides).
  std::vector<std::vector<double> > outDirection; // Outward-pointing normal vector (2 sides).
  std::vector<bool> isInterior; // Are these cell boundaries interior to the grid, or overall boundaries?
  std::vector<double> bdryLocationx; // Coordinates of the 4 cell boundaries.
  std::vector<double> bdryLocationy; // Coordinates of the 4 cell boundaries.

};


// Class for a grid made up of many cells:
class DataGrid{
 public:
  // Constructor (again, just prototyped here) that takes overall left and right boundaries, number of cells, a flag
  // to specify whether the domain is periodic, and (optionally) data values in each cell:
  DataGrid(double x1, double x2, double y1, double y2, int Nx, int Ny, bool periodic, int Nfields, std::vector<std::vector<double> > values);

  // Member function to compute the net flux out of each cell:
  std::vector<std::vector<double> > NetFlux(double, double, std::vector<std::vector<double> > (*f)(std::vector<double>));

  // Member function to use that flux to find the rate of change of the cell averages:
  std::vector<std::vector<double> > TimeDeriv(double dx, double dy, std::vector<std::vector<double> > (*f)(std::vector<double>), std::vector<double> (*s)(std::vector<double>));
  
  // Define the member variables:
  std::vector<cell> Cell;
  double Bdry1x;
  double Bdry2x;
  double Bdry1y;
  double Bdry2y;
  int Ncellsx;
  int Ncellsy;
  int Ncells;
  int Nfields;
  bool IsPeriodic;
  
};

// Simple class for packaging wave fields:
class WaveFields{
 public:
  WaveFields(){
    psi = 0.;
    pi = 0.;
    phi = std::vector<double>{0.,0.};
  }
  
  double psi;
  double pi;
  std::vector<double> phi;
  
};

// Simple class for a set of characteristic fields and associated speeds:
class CharFields{
 public:
  CharFields(int); // Simple constructor; the int is the number of fields.

  // Member variables:
  int numfields;
  std::vector<double> fields;
  std::vector<double> speeds;
};
  
// Prototypes of functions introduced for characteristic decomposition and upwinding:
std::vector<double> CharUpwindRecon(const std::vector<double>& vals, const std::vector<double>& neighborvals, std::vector<double> nhat, double upwindparam);
CharFields BasicToChar(const std::vector<double>& fields, std::vector<double> nhat);
std::vector<double> CharToBasic(const std::vector<double>& u, std::vector<double> nhat);

  

// Finally, let's prototype some binary operators and related convenience functions:

// Scalar product of two std::vector<double>'s:
double DotProd(std::vector<double>, std::vector<double>);

// Addition operation to add a vector to datagrid values:
DataGrid operator+(const DataGrid& Grid, std::vector<double> Addend);
DataGrid operator+(const DataGrid& Grid, std::vector<std::vector<double> > Addend);

// Operation to multiply an array by a double:
std::vector<double> operator*(double Multiplier, const std::vector<double>& Array);
std::vector<std::vector<double> > operator*(double Multiplier, const std::vector<std::vector<double> >& Array);


// Operation to add two std::vector<double>'s:
std::vector<double> operator+(const std::vector<double>&, const std::vector<double>&);
std::vector<double> operator-(const std::vector<double>&, const std::vector<double>&);
std::vector<double> operator*(const std::vector<double>&, const std::vector<double>&);

// Operation to add two std::vector<std::vector<double> >'s:
std::vector<std::vector<double> > operator+(const std::vector<std::vector<double> >&, const std::vector<std::vector<double> >&);

