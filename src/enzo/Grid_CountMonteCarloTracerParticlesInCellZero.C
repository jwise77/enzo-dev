/***********************************************************************
/
/  GRID CLASS (Count Monte Carlo Tracer Particles in Cell Zero)
/
/  written by: Corey Brummel-Smith
/  date:       June, 2021
/
/  PURPOSE: Count the total number of Monte Carlo tracer particles 
/           in the first cell of this grid. The only time cell 0 has 
/           any MC tracer particles is when they've been inserted there
/           by CommunicationCollectMonteCarloTracerParticles(COPY_IN)
/           to be moved to their proper grid.
/
/  RETURNS:
/    N: Number of MC tracers
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

int grid::CountMonteCarloTracerParticlesInCellZero()
{
  int N = 0;
  MonteCarloTracerParticle *mc = MonteCarloTracerParticles[0];
  while(mc != NULL){
    N++;
    mc = mc->NextParticle;
  }
  return N;
}