/***********************************************************************
/
/  GRID CLASS (Count Monte Carlo Tracer Particles)
/
/  written by: Corey Brummel-Smith
/  date:       June, 2021
/
/  PURPOSE: Count the total number of Monte Carlo tracer particles 
/           in this grid and update the NumberOfMonteCarloTracerParticles
/           grid variable.
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

int grid::CountMonteCarloTracerParticles()
{
  /* Compute grid size. */
  int i, dim, size = 1, N = 0;
  MonteCarloTracerParticle *mc;

  if (MonteCarloTracerParticles == NULL) {
    //printf("\nproc%0d: CountMCTP: MonteCarloTracerParticles not allocated on grid %d, ProcessorNumber %d\n", MyProcessorNumber, this->ID, ProcessorNumber);
    return 0;
  }

  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];   

  for (i = 0; i < size; i++){
    mc = MonteCarloTracerParticles[i];
    while(mc != NULL){
      N++;
      mc = mc->NextParticle;
    }
  }
  this->NumberOfMonteCarloTracerParticles = N;
  return N;
}