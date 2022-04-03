
#pragma once

#include <vector>

// Class for data in a single finite-volume cell:
class cell{
 public:
  // Constructor prototype (I'll actually define the constructor in DataGrid.cpp):
  cell(int Ind, double Cent, double Val, std::vector<int> Neighbors, 
  std::vector<double> OutDir, std::vector<bool> IsInterior, std::vector<double> BdryLocation,
  double flowDirection=1.);

  // Member variables:
  //private: // (for simplicity, I'll break protocol and let all member variables be public)
  int index; // A label for the cell in a grid.
  double center; // Coordinate of the cell center.
  double value;  // Cell's average field value.
  double flowDirection; // Whether the field flows in the positive or negative direction (we'll improve this implementation next week)
  std::vector<int> neighbors; // Cells that neighbor this one in a grid (2 sides).
  std::vector<double> outDirection; // Outward-pointing normal vector (2 sides).
  std::vector<bool> isInterior; // Are these cell boundaries interior to the grid, or overall boundaries?
  std::vector<double> bdryLocation; // Coordinates of the 2 cell boundaries.

};


// Class for a grid made up of many cells:
class DataGrid{
 public:
  // Constructor (again, just prototyped here) that takes overall left and right boundaries, number of cells, a flag
  // to specify whether the domain is periodic, and (optionally) data values in each cell:
  DataGrid(double x1, double x2, int N, bool periodic, std::vector<double> values, 
  double flowDirection=1.);

  std::vector<double> Data();

  // Member function to compute the net flux out of each cell:
  std::vector<double> NetFlux(double (*f)(double));

  // Member function to use that flux to find the rate of change of the cell averages:
  std::vector<double> Dot(double dx, double (*f)(double));
  
  // Define the member variables:
  std::vector<cell> Cell;
  double Bdry1;
  double Bdry2;
  int Ncells;
  bool IsPeriodic;
  
};

// Finally, let's prototype a couple of binary operators that'll be useful:

// Addition operation to add a vector to datagrid values:
DataGrid operator+(DataGrid& Grid, std::vector<double> Addend);
// Note: the ampersand here simply tells the function to "pass by reference", which can 
//be slightly faster (or much faster, in some cases).

// Operation to multiply an array by a double:
std::vector<double> operator*(double Multiplier, const std::vector<double>& Array);


// Operation to add two std::vector<double>'s:
std::vector<double> operator+(std::vector<double>, std::vector<double>);
