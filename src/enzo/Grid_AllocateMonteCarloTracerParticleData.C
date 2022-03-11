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

int grid::AllocateMonteCarloTracerParticleData()
{

  if (this->MonteCarloTracerParticles != NULL) {
    printf("proc%d: AllocateMonteCarloTracerParticleData: Already allocated. Returning now.", MyProcessorNumber);
    return SUCCESS;
  }
  /* Compute grid size. */
  int size = 1;
  int i, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];   
  
  // allocate memory for array ("grid") of MC tracer lists pointers
  this->MonteCarloTracerParticles = new MonteCarloTracerParticle*[size];

  for (dim = 0; dim < GridRank; dim++)
 	 this->MassFlux[dim] = new float[size];

  for (int i = 0; i < size; i++){
    this->MonteCarloTracerParticles[i] = NULL;
    for (dim = 0; dim < GridRank; dim++)
      this->MassFlux[dim][i] = 0.0;
  }  
  return SUCCESS;
}
