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

void grid::AllocateMonteCarloTracerParticleData():
{
  /* Compute grid size. */
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];   
  
  // allocate memory for array ("grid") of MC tracer lists pointers
  this->MonteCarloTracerParticles = new *MonteCarloTracerParticle[size];

  for (int i = 0; i < size; i++)
    this->MonteCarloTracerParticles[i] = NULL;

  this->ReducedCellMass = new float[size];

  for (int dim = 0; dim < GridRank; dim++)
 	 this->MassFlux[dim] = new float[size];

}