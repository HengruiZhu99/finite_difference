
#include <vector>
#include <iostream>
#include "DataGridMultiField.hpp"


// The constructor for the "cell" class:
cell::cell(int Ind, double Cent, int Nfieldsgiven, std::vector<double> Vals, 
std::vector<int> Neighbors, std::vector<double> OutDir, std::vector<bool> IsInterior, 
std::vector<double> BdryLocation, std::vector<double> flowDirection){
  
  index = Ind;
  center = Cent;
  Nfields = Nfieldsgiven;
  for(int alpha=0; alpha<Nfields; alpha++){
    values.push_back(Vals[alpha]);
    flowdirection.push_back(flowDirection[alpha]);
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
DataGrid::DataGrid(double x1, double x2, int N, bool periodic=0, int Nfieldsgiven=3, 
std::vector<std::vector<double> > values = std::vector<std::vector<double> >{}, 
std::vector<double> flowdirection = std::vector<double>{0.,1.,-1.}){
  // Set the simplest member variables immediately:
  Bdry1 = x1;
  Bdry2 = x2;
  Ncells = N;
  Nfields = Nfieldsgiven;
  IsPeriodic = periodic;
  for(int alpha = 0; alpha<Nfields; alpha++){
    flowdirection.push_back(flowdirection[alpha]); 
  }
  
  // The boundary locations for each cell:
  std::vector<double> BdryPoints(N+1, 0.); // (this initializes to an array with N+1 elements (bdry's of N neighboring cells, all initially set to zero)
  // But let's set the actual cell boundary locations:
  for(int i=0; i<N+1; i++){
    BdryPoints[i] = x1 + (x2-x1)*double(i)/double(N);
  }
  // Cycle over cells:
  for(int i=0; i<N; i++){
    // Set cell value either to given value, or zero if nothing was given:

    //int Nfields = flowdirection.size();
    std::vector<double> value{};
    for(int alpha = 0; alpha<Nfields; alpha++){
      value.push_back((values.size()==N)?values[i][alpha]:0.);
    }
    //double value = (values.size()==N)?values[i]:0.; 
    
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
    cell TheCell(i, centerpoint, Nfields, value, Neighbors, OutDir, IsInterior, BdryLocation, flowdirection);
    Cell.push_back(TheCell); // The vector member function "push_back" should be called "append", as far as I'm concerned.
  }
}


// Datagrid member function to compute the net flux out of each cell:
std::vector<std::vector<double> > DataGrid::NetFlux(std::vector<double> (*fluxfunc)(std::vector<double>)){
  std::vector<std::vector<double> > flux{};
  // Initialize and set size:
  //std::cout << "Nfields: " << Nfields << std::endl;
  for(int i=0; i<Ncells; i++){
    flux.push_back(std::vector<double>(Nfields));
  }

  // Loop over every cell:
  for(int i=0; i<Ncells; i++){
    for(int alpha = 0; alpha<Nfields; alpha++){
      for(int bdry = 0; bdry<2; bdry++){
	if(Cell[i].isInterior[bdry]==1||1){
	  // For the advection equation that we're solving for now, the flux is simply the field value,
	  // evaluated on the boundary between any two cells. The simplest approximation of this is
	  // simply the average of the value in "this" cell and the neighboring cell. This is called
	  // the "centered flux," and, as we'll see, it's a little problematic.


	  if(Cell[i].isInterior[bdry]==1) flux[i][alpha] += .5*(fluxfunc(Cell[i].values)[alpha] + fluxfunc(Cell[Cell[i].neighbors[bdry]].values)[alpha])*Cell[i].outDirection[bdry];
	  // Note that we also multiply by the outgoing unit vector on both boundaries. On the left
	  // boundary of the cell, the flux is inward. On the right boundary, the flux is outward.
	  
	  // For stability reasons, it turns out to be better to weight the above average more
	  // towards the side where data is coming in. This is called the "upwind flux."
	  double upwinding = 1.;// 1.0 for true upwinding. 0 for FTCS. -1 for downwinding.
	  if(Cell[i].isInterior[bdry]==1) flux[i][alpha] += upwinding*.5*(fluxfunc(Cell[i].values)[alpha] - fluxfunc(Cell[Cell[i].neighbors[bdry]].values)[alpha])*Cell[i].flowdirection[alpha];
	  if(Cell[i].isInterior[bdry]==0 && Cell[i].outDirection[bdry]*Cell[i].flowdirection[alpha]<0.){
	    flux[i][alpha] = 0.;

	    //The way I'm dealing with nonperiodic boundaries here is hacky and still broken!
	    
	  }
	}// End conditional for internal bdry.
      }
    }
  }
  return flux;
}

// Member function to calculate the rate of change of some DataGrid:
std::vector<std::vector<double> > DataGrid::Dot(double dx, std::vector<double> (*fluxfunc)(std::vector<double>), std::vector<double> (*sourcefunc)(std::vector<double>)){
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

// // Operation to multiply an array by a double:
// std::vector<double> operator*(double Multiplier, const std::vector<double>& Array){
//   std::vector<double> result;
//   for(int i=0; i<Array.size(); i++){
//     result.push_back(Multiplier*Array[i]);
//   }
//   return result;
// }

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
std::vector<double> operator+(std::vector<double> V1, std::vector<double> V2){
  std::vector<double> Sum(V1.size());
  for(int i=0; i<V1.size(); i++){
    Sum[i] = V1[i]+V2[i];
  }
  return Sum;
}
// Operation to add two std::vector<std::vector<double> >'s:
std::vector<std::vector<double> > operator+(std::vector<std::vector<double> > V1, std::vector<std::vector<double> > V2){
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
