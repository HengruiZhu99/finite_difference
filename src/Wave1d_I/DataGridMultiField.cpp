
#include <vector>
#include <iostream>
#include <cassert>
#include "DataGridMultiField.hpp"


// Constructor for "Characteristic fields" class:
CharFields::CharFields(int N){
  // Define set the member variable for the number of fields:
  numfields = N;
  // Initialize the member variable vectors for field values and associated speeds:
  fields = std::vector<double>(N);
  speeds = std::vector<double>(N);
}


// The constructor for the "cell" class:
cell::cell(int Ind, double Cent, int Nfieldsgiven, std::vector<double> Vals, std::vector<int> Neighbors, std::vector<double> OutDir, std::vector<bool> IsInterior, std::vector<double> BdryLocation){
  index = Ind;
  center = Cent;
  Nfields = Nfieldsgiven;
  for(int alpha=0; alpha<Nfields; alpha++){
    values.push_back(Vals[alpha]);
  }
  
  for(int i=0; i<2; i++){
    neighbors.push_back(Neighbors[i]);
    outDirection.push_back(OutDir[i]);
    isInterior.push_back(IsInterior[i]);
    bdryLocation.push_back(BdryLocation[i]);
  }
  // I'm hardcoding the assumption that each cell has 2 boundaries.
  // This will have to change when upgrading to higher dimensions.
}



// Constructor for the "DataGrid" class that takes overall left and right boundaries, number of cells, a flag
// to specify whether the domain is periodic, and (optionally) data values in each cell:
DataGrid::DataGrid(double x1, double x2, int N, bool periodic=0, int Nfieldsgiven=3, std::vector<std::vector<double> > values = std::vector<std::vector<double> >{}){//, std::vector<double> flowdirection = std::vector<double>{0.,1.,-1.}){
  // Set the simplest member variables immediately:
  Bdry1 = x1;
  Bdry2 = x2;
  Ncells = N;
  Nfields = Nfieldsgiven;
  IsPeriodic = periodic;
  
  // The boundary locations for each cell:
  std::vector<double> BdryPoints(N+1, 0.); // (this initializes to an array with N+1 elements (bdry's of N neighboring cells, all initially set to zero)
  // But let's set the actual cell boundary locations:
  for(int i=0; i<N+1; i++){
    BdryPoints[i] = x1 + (x2-x1)*double(i)/double(N);
  }
  // Cycle over cells:
  for(int i=0; i<N; i++){
    // Set cell value either to given value, or zero if nothing was given:

    std::vector<double> value{}; // Start with a vector with no elements.
    for(int alpha = 0; alpha<Nfields; alpha++){
      value.push_back((values.size()==N)?values[i][alpha]:0.);
      // Append either the given values, or zero by default.
    }
    
    double centerpoint = .5*(BdryPoints[i]+BdryPoints[i+1]);

    // Define indices for neighboring cells:
    std::vector<int> Neighbors{i-1,i+1};

    // If it's a periodic domain, the left neighbor of the leftmost cell is the rightmost cell:
    if(periodic==1 && i==0) Neighbors[0] = N-1;
    // And vice versa:
    if(periodic==1 && i==N-1) Neighbors[1] = 0;

    // The outgoing direction on both boundaries of the given cell:
    std::vector<double> OutDir{-1., 1.};

    // Most boundaries of most cells will be interior:
    std::vector<bool> IsInterior{1,1};

    // But on the first and last cell, if non-periodic, the boundaries will be exterior:
    if(periodic==0 && i==0) IsInterior[0] = 0;
    if(periodic==0 && i==N-1) IsInterior[1] = 0;

    // Locations of this cell's two boundaries:
    std::vector<double> BdryLocation{BdryPoints[i], BdryPoints[i+1]};

    // Finally, define a cell with these attributes and append it to the array of Cell values:
    cell TheCell(i, centerpoint, Nfields, value, Neighbors, OutDir, IsInterior, BdryLocation);
    Cell.push_back(TheCell); // The vector member function "push_back" should be called "append", in my opinion.
  }
}


// Datagrid member function to compute the net flux out of each cell:
std::vector<std::vector<double> > DataGrid::NetFlux(std::vector<double> (*fluxfunc)(std::vector<double>)){
  std::vector<std::vector<double> > flux{};
  // Initialize and set size:
  for(int i=0; i<Ncells; i++){
    flux.push_back(std::vector<double>(Nfields));
  }

  // Loop over every cell:
  for(int i=0; i<Ncells; i++){
    for(int bdry = 0; bdry<2; bdry++){
      if(Cell[i].isInterior[bdry]==1){

	// Upwind field reconstruction using characteristic decomposition:
	double upwindparam = 0.1;
       	std::vector<double> ReconField = CharUpwindRecon(Cell[i].values, Cell[Cell[i].neighbors[bdry]].values, Cell[i].outDirection[bdry], upwindparam);

	// Add the contribution to the flux from this boundary:
	flux[i] = flux[i] + Cell[i].outDirection[bdry]*fluxfunc(ReconField);
	

	  
      }// End conditional for internal bdry.
      if(Cell[i].isInterior[bdry]==0){
	std::cout << "Currently unable to handle exterior boundaries!" << std::endl;
	assert(Cell[i].isInterior[bdry]==1);
      }
    }
  }
  return flux;
}

// The function that uses characteristic decomposition to do boundary reconstruction for all fields:
std::vector<double> CharUpwindRecon(const std::vector<double>& vals, const std::vector<double>& neighborvals, double nhat, double upwindparam){
  // Convert from basic fields to characteristic fields:
  CharFields u = BasicToChar(vals, nhat);
  CharFields u_neighbor = BasicToChar(neighborvals, nhat);

  // First reconstruct with a simple average of neighboring cell values:
  std::vector<double> u_recon = .5*(u.fields + u_neighbor.fields);

  // Now, use characteristic field speeds to bias the average toward the upwind side:
  u_recon = u_recon + upwindparam*0.5*(u.speeds*(u.fields - u_neighbor.fields));
    
  // Finally, convert back to basic fields:
  return CharToBasic(u_recon, nhat);
}

// The function that converts the basic fields (psi, pi, phi) into the propagating fields (psi, u+, u-) across some surface:
CharFields BasicToChar(const std::vector<double>& fields, double nhat){
  // Reminder: CharFields is a class that holds a vector for the
  // field values and another vector for the associated speeds.
  CharFields u(fields.size()); // (Initialize with zero-vectors.)
  // Hardcoded char decomposition for wave equation:
  u.fields[0] = fields[0];  // u[0] = psi
  u.speeds[0] = 0.;  // zero-speed field
  u.fields[1] = fields[1] - nhat*fields[2]; // u[+] = pi - nhat . phi
  u.speeds[1] = 1.;  // outgoing field
  u.fields[2] = fields[1] + nhat*fields[2]; // u[-] = pi + nhat . phi
  u.speeds[2] = -1.; // ingoing field
  return u;
}

std::vector<double> CharToBasic(const std::vector<double>& ufields, double nhat){
  // Note that, for convenience, I'm passing a std::vector<double>, just the characteristic fields.
  // There's a cleaner way to do this, I'm sure.
  std::vector<double> fields(ufields.size());
  // psi:
  fields[0] = ufields[0];
  // pi = .5*(u[+] + u[-]):
  fields[1] = .5*(ufields[1] + ufields[2]);
  // phi = .5*nhat*(u[-] - u[+]):
  fields[2] = .5*nhat*(ufields[2] - ufields[1]);
  return fields;
}



// Member function to calculate the rate of change of some DataGrid:
std::vector<std::vector<double> > DataGrid::TimeDeriv(double dx, std::vector<double> (*fluxfunc)(std::vector<double>), std::vector<double> (*sourcefunc)(std::vector<double>)){
  std::vector<std::vector<double> > Flux = NetFlux(fluxfunc);
  std::vector<std::vector<double> > Result;
  
  for(int i=0; i<Ncells; i++){
    Result.push_back(std::vector<double>(Nfields));
    std::vector<double> Source = sourcefunc(Cell[i].values);
    for(int alpha=0; alpha<Nfields; alpha++){
      Result[i][alpha] = -(1./dx)*Flux[i][alpha] + Source[alpha];
    }
  }
  return Result;
}




  

// Addition operation to add a vector to datagrid values:
DataGrid operator+(const DataGrid& Grid, std::vector<std::vector<double> > Addend){
  DataGrid result(Grid); // Construct "result" as a copy of "grid"
  for(int i=0; i<Grid.Ncells; i++){
    for(int alpha = 0; alpha<Grid.Nfields; alpha++){
      result.Cell[i].values[alpha] = Grid.Cell[i].values[alpha] + Addend[i][alpha];
    }
  }
  return result;
}

// Operation to multiply an array by a double:
std::vector<double> operator*(double Multiplier, const std::vector<double>& Array){
  std::vector<double> result;
  for(int i=0; i<Array.size(); i++){
    result.push_back(Multiplier*Array[i]);
  }
  return result;
}

// Operation to multiply a multi-field array by a double:
std::vector<std::vector<double> > operator*(double Multiplier, const std::vector<std::vector<double> >& Array){
  std::vector<std::vector<double> > result;
  for(int i=0; i<Array.size(); i++){
    std::vector<double> tmp{};
    for(int alpha = 0; alpha<Array[0].size(); alpha++){
      tmp.push_back(Multiplier*Array[i][alpha]);
    }
    result.push_back(tmp);
  }
  return result;
}


// Operation to add two std::vector<double>'s:
std::vector<double> operator+(const std::vector<double>& V1, const std::vector<double>& V2){
  std::vector<double> Sum(V1.size());
  for(int i=0; i<V1.size(); i++){
    Sum[i] = V1[i]+V2[i];
  }
  return Sum;
}

// Operation to subtract two std::vector<double>'s:
std::vector<double> operator-(const std::vector<double>& V1, const std::vector<double>& V2){
  std::vector<double> Diff(V1.size());
  for(int i=0; i<V1.size(); i++){
    Diff[i] = V1[i]-V2[i];
  }
  return Diff;
}

// Operation to elementwise multiply two std::vector<double>'s:
std::vector<double> operator*(const std::vector<double>& V1, const std::vector<double>& V2){
  std::vector<double> Prod(V1.size());
  for(int i=0; i<V1.size(); i++){
    Prod[i] = V1[i]*V2[i];
  }
  return Prod;
}

// Operation to add two std::vector<std::vector<double> >'s:
std::vector<std::vector<double> > operator+(const std::vector<std::vector<double> >& V1, const std::vector<std::vector<double> >& V2){
  std::vector<std::vector<double> > Sum{};
  for(int i=0; i<V1.size(); i++){
    std::vector<double> tmp{};
    for(int alpha=0; alpha<V1[0].size(); alpha++){
      tmp.push_back(V1[i][alpha] + V2[i][alpha]);
    }
    Sum.push_back(tmp);
  }
  return Sum;
}
