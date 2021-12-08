/***********************************************************************
/
/  GRID CLASS (RETURN THE NUMBER OF MONTE CARLO TRACER PARTICLES IN THE
/              FIRST GHOST CELLS)
/
/  written by: Corey Brummel-Smith
/  date:       Noveber, 2021
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"

int grid::CountMonteCarloTracerParticlesInFirstGhostCells()
{

  int icell, nmc = 0;
  MonteCarloTracerParticle *mctp;
  
  const int NUMBER_OF_FACE_CELLS = 2 * (GridDimension[0] * GridDimension[1] +
                                        GridDimension[0] * GridDimension[2] +
                                        GridDimension[1] * GridDimension[2]);
  
  int FirstGhostCellIndices = int[NUMBER_OF_FACE_CELLS];
  FirstGhostCellIndices = this->GetFirstGhostCells(NUMBER_OF_FACE_CELLS);
  
  // Loop over first ghost cells
  for (icell = 0; icell < NUMBER_OF_FACE_CELLS; icell++) {
    
    index = FirstGhostCellIndices[icell];
    mctp = MonteCarloTracerParticles[index];
  
    // Loop over MC tracers in this cell
    while(mctp != NULL){
      nmc++;
      mctp = mctp->NextParticle;
    }
  }
  return nmc;
}