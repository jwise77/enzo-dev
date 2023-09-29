/***********************************************************************
/
/  RESET MONTE CARLO TRACER PARTICLES EXCHANGED THIS TIMESTEP
/
/  written by: Corey Brummel-Smith
/  date:       September, 2021
/
/  PURPOSE: Reset ExchangedThisTimestep for all Monte Carlo tracer
/           particles on all cells on this grid.
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h" 
#include "MonteCarloTracerParticle.h" 

int grid::ResetMonteCarloTracerParticlesExchangedThisTimestep()
{

  if (MonteCarloTracerParticles == NULL)
    return SUCCESS;
  
  /* Compute grid size. */
  int i, dim, size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];   

  for (i = 0; i < size; i++){
    MonteCarloTracerParticle *mc = MonteCarloTracerParticles[i];
    while(mc != NULL){
      mc->ExchangedThisTimestep = false;
      mc = mc->NextParticle;
    }
  }

  return SUCCESS;
}
