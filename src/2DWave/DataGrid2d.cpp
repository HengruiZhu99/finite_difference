
#include <vector>
#include <iostream>
#include <cassert>
#include "DataGrid2d.hpp"
#include <omp.h>


// Constructor for "Characteristic fields" class:
CharFields::CharFields(int N){
  // Define set the member variable for the number of fields:
  numfields = N;
  // Initialize the member variable vectors for field values and associated speeds:
  fields = std::vector<double>(N);
  speeds = std::vector<double>(N);
}


// The constructor for the "cell" class:
cell::cell(int Ind, double Cent_x, double Cent_y, int Nfieldsgiven, std::vector<double> Vals, std::vector<int> Neighbors, std::vector<std::vector<double> > OutDir, std::vector<bool> IsInterior, std::vector<double> BdryLocationx, std::vector<double> BdryLocationy){
  index = Ind;
  //std::cout << "    The index we got was " << index << std::endl;
  centerx = Cent_x;
  centery = Cent_y;
  Nfields = Nfieldsgiven;
  for(int alpha=0; alpha<Nfields; alpha++){
    values.push_back(Vals[alpha]);
  }

  //std::cout << "Ind = " << Ind << ", Neighbors[0] = " << Neighbors[0] << std::endl;
  
  for(int i=0; i<4; i++){
    // Hardcoding 4 neighbors and 4 boundaries for each cell, assuming 2d.
    neighbors.push_back(Neighbors[i]);
    outDirection.push_back(OutDir[i]);
    isInterior.push_back(IsInterior[i]);
    bdryLocationx.push_back(BdryLocationx[i]);
    bdryLocationy.push_back(BdryLocationy[i]);
  }

  //std::cout << "Ind = " << Ind << ", neighbors[0] = " << neighbors[0] << std::endl;

}



// Constructor for the "DataGrid" class that takes overall left and right boundaries, number of cells, a flag
// to specify whether the domain is periodic, and (optionally) data values in each cell:
DataGrid::DataGrid(double x1, double x2, double y1, double y2, int Nx, int Ny, bool periodic=0, int Nfieldsgiven=4, std::vector<std::vector<double> > values = std::vector<std::vector<double> >{}){//, std::vector<double> flowdirection = std::vector<double>{0.,1.,-1.}){
  // Set the simplest member variables immediately:
  Bdry1x = x1;
  Bdry2x = x2;
  Bdry1y = y1;
  Bdry2y = y2;
  Ncellsx = Nx;
  Ncellsy = Ny;
  Ncells = Nx*Ny;
  Nfields = Nfieldsgiven;
  IsPeriodic = periodic;
  
  // The boundary locations for each cell:
  //std::vector<double> BdryPointsx((Nx+1)*(Ny+1), 0.);
  //std::vector<double> BdryPointsy((Nx+1)*(Ny+1), 0.);
  std::vector<double> BdryPointsx(Nx+1, 0.);
  std::vector<double> BdryPointsy(Ny+1, 0.);
  // for(int i=0; i<Nx+1; i++){
  //   for(int j=0; j<Ny+1; j++){
  //     int I = i+Nx*j;
  //     BdryPointsx[I] = x1 + (x2-x1)*double(i)/double(Nx);
  //     BdryPointsy[I] = y1 + (y2-y1)*double(j)/double(Ny);
  //   }
  // }
  for(int i=0; i<Nx+1; i++){
    BdryPointsx[i] = x1 + (x2-x1)*double(i)/double(Nx);
  }
  for(int j=0; j<Ny+1; j++){
    BdryPointsy[j] = y1 + (y2-y1)*double(j)/double(Ny);
  }

  
  // Cycle over cells:
  for(int j=0; j<Ny; j++){
    for(int i=0; i<Nx; i++){
      //for(int i=0; i<Nx; i++){
      //for(int j=0; j<Ny; j++){

      // A master-index that I'll repeatedly use to label all the grid points.
      int I = i+Nx*j;

      // Set cell value either to given value, or zero if nothing was given:
            std::vector<double> value{}; // Start with a vector with no elements.
      for(int alpha = 0; alpha<Nfields; alpha++){
	value.push_back((values.size()==Ncells)?values[I][alpha]:0.);
	// Append either the given values, or zero by default.
      }
      
      double centerpoint_x = .5*(BdryPointsx[i]+BdryPointsx[i+1]);
      double centerpoint_y = .5*(BdryPointsy[j]+BdryPointsy[j+1]);
     
      // Define indices for neighboring cells:
      std::vector<int> Neighbors(4);
      Neighbors[0] = (i+1) + Nx*j;
      Neighbors[1] = i + Nx*(j+1);
      Neighbors[2] = (i-1) + Nx*j;
      Neighbors[3] = i + Nx*(j-1);
      
      
      // If it's a periodic domain, the left neighbor of the leftmost cell is the rightmost cell:
      if(periodic==1 && i==0) Neighbors[2] = Nx-1 + Nx*j;
      // And vice versa:
      if(periodic==1 && i==Nx-1) Neighbors[0] = 0 + Nx*j;
      if(periodic==1 && j==0) Neighbors[3] = i + Nx*(Ny-1);
      if(periodic==1 && j==Ny-1) Neighbors[1] = i + Nx*0;
      
      // The outgoing direction on both boundaries of the given cell:
      std::vector<std::vector<double> > OutDir{std::vector<double>{1.,0.},
	  std::vector<double>{0.,1.},
	    std::vector<double>{-1.,0.},
	      std::vector<double>{0.,-1.}};


      
      // Most boundaries of most cells will be interior:
      std::vector<bool> IsInterior{1,1,1,1};
      
      // But on the first and last cell, if non-periodic, the boundaries will be exterior:
      if(periodic==0 && i==0) IsInterior[2] = 0;
      if(periodic==0 && i==Nx-1) IsInterior[0] = 0;
      if(periodic==0 && j==0) IsInterior[3] = 0;
      if(periodic==0 && j==Ny-1) IsInterior[1] = 0;
     
      // Locations of this cell's four boundaries:
      std::vector<double> BdryLocationx{BdryPointsx[I], BdryPointsx[I+1]};
      std::vector<double> BdryLocationy{BdryPointsy[I], BdryPointsy[I+Nx]};
      // There could be errors in this BdryLocation infrastructure, but for now
      // we wouldn't see it, because I don't think it's used for anything.
      //
      // Note: this is bad coding practice. I should either test it to make sure
      // it works, or I should take it out. 

      // Finally, define a cell with these attributes and append it to the array of Cell values:
      cell TheCell(I, centerpoint_x, centerpoint_y, Nfields, value, Neighbors, OutDir, IsInterior, BdryLocationx, BdryLocationy);
      Cell.push_back(TheCell); // The vector member function "push_back" should be called "append", in my opinion.
    }
  }
}



// Datagrid member function to compute the net flux out of each cell:
std::vector<std::vector<double> > DataGrid::NetFlux(double dx, double dy, std::vector<std::vector<double> > (*fluxfunc)(std::vector<double>)){
  std::vector<std::vector<double> > flux{};
  // Initialize and set size:
  for(int i=0; i<Ncells; i++){
    flux.push_back(std::vector<double>(Nfields));
  }

  // Loop over every cell:
  #pragma omp parallel for
  for(int i=0; i<Ncells; i++){
    for(int bdry = 0; bdry<4; bdry++){
      if(Cell[i].isInterior[bdry]==1){

	// Upwind field reconstruction using characteristic decomposition:
	double upwindparam = 0.1;
       	std::vector<double> ReconField = CharUpwindRecon(Cell[i].values, Cell[Cell[i].neighbors[bdry]].values, Cell[i].outDirection[bdry], upwindparam);

	// Add the contribution to the flux from this boundary:
	std::vector<double> n = Cell[i].outDirection[bdry];
	std::vector<std::vector<double> > FF = fluxfunc(ReconField);
	//flux[i] = flux[i] + Cell[i].outDirection[bdry]*fluxfunc(ReconField);
	double facesize = dy;
	if(bdry==1 || bdry == 3){facesize = dx;}
	flux[i] = flux[i] + facesize*(n[0]*FF[0] + n[1]*FF[1]);

	  
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
std::vector<double> CharUpwindRecon(const std::vector<double>& vals, const std::vector<double>& neighborvals, std::vector<double> nhat, double upwindparam){
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
CharFields BasicToChar(const std::vector<double>& fields, std::vector<double> nhat){
  // Reminder: CharFields is a class that holds a vector for the
  // field values and another vector for the associated speeds.
  CharFields u(fields.size()+1); // (Initialize with zero-vectors.)
  std::vector<double> phi{fields[2],fields[3]};
  double phidotn = DotProd(nhat, phi);
  std::vector<double> phitrans = phi - phidotn*nhat;
  // Hardcoded char decomposition for wave equation:
  u.fields[0] = fields[0];  // u[0] = psi
  u.speeds[0] = 0.;  // zero-speed field
  u.fields[1] = fields[1] - DotProd(nhat,phi); // u[+] = pi - nhat . phi
  u.speeds[1] = 1.;  // outgoing field
  u.fields[2] = fields[1] + DotProd(nhat,phi); // u[-] = pi + nhat . phi
  u.speeds[2] = -1.; // ingoing field
  u.fields[3] = phitrans[0];
  u.speeds[3] = 0.;
  u.fields[4] = phitrans[1];
  u.speeds[4] = 0.;
  return u;
}

std::vector<double> CharToBasic(const std::vector<double>& ufields, std::vector<double> nhat){
  // Note that, for convenience, I'm passing a std::vector<double>, just the characteristic fields.
  // There's a cleaner way to do this, I'm sure.
  std::vector<double> fields(ufields.size()-1);
  // psi:
  fields[0] = ufields[0];
  // pi = .5*(u[+] + u[-]):
  fields[1] = .5*(ufields[1] + ufields[2]);
  // phi = .5*nhat*(u[-] - u[+]) + phitrans:
  fields[2] = .5*nhat[0]*(ufields[2] - ufields[1])+ufields[3];
  fields[3] = .5*nhat[1]*(ufields[2] - ufields[1])+ufields[4];

  return fields;
}



// Member function to calculate the rate of change of some DataGrid:
std::vector<std::vector<double> > DataGrid::TimeDeriv(double dx, double dy, std::vector<std::vector<double> > (*fluxfunc)(std::vector<double>), std::vector<double> (*sourcefunc)(std::vector<double>)){
//std::vector<std::vector<double> > DataGrid::TimeDeriv(double dx, std::vector<WaveFields> (*fluxfunc)(WaveFields), WaveFields (*sourcefunc)(WaveFields)){
  std::vector<std::vector<double> > Flux = NetFlux(dx, dy, fluxfunc);
  std::vector<std::vector<double> > Result;
  
  for(int i=0; i<Ncells; i++){
    Result.push_back(std::vector<double>(Nfields));
    std::vector<double> Source = sourcefunc(Cell[i].values);
    for(int alpha=0; alpha<Nfields; alpha++){
      Result[i][alpha] = -(1./(dx*dy))*Flux[i][alpha] + Source[alpha];
    }
  }
  return Result;
}





// Scalar-product function:
double DotProd(std::vector<double> V1, std::vector<double> V2){
  double result = 0.;
  for(int i=0; i<V1.size(); i++){
    result += V1[i]*V2[i];
  }
  return result;
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
