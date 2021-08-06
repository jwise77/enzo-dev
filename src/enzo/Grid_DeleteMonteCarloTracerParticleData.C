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

void DeleteMonteCarloTracerParticleList(MonteCarloTracerParticle *&Node);

int grid::DeleteMonteCarloTracerParticleData()
{
  /* Compute grid size. */
  int size = 1;
  int i, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];   
  
  for (i = 0; i < size; i++) {
    if (this->MonteCarloTracerParticles[i] != NULL)
      DeleteMonteCarloTracerParticleList(this->MonteCarloTracerParticles[i]);
  }
  delete [] this->MonteCarloTracerParticles;
  for (dim = 0; dim < GridRank; dim++)    
    delete [] this->MassFlux[dim];
  return SUCCESS;
}