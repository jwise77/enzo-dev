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

void InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticle * &Node, MonteCarloTracerParticle * &NewNode);

int grid::CreateMonteCarloTracerParticles()
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int groupID = 0;
  int size = 1;
  int i, j, k, n, index;
  FLOAT pos[3];
  MonteCarloTracerParticle *newMC;

  this->AllocateMonteCarloTracerParticleData();

  // Loop over all cells
  for (k = 0; k < GridDimension[2]; k++) {
    pos[2] = (k + 0.5) * CellWidth[2][0];
    for (j = 0; j < GridDimension[1]; j++) {
      pos[1] = (j + 0.5) * CellWidth[1][0];
      for (i = 0; i < GridDimension[0]; i++) {
        pos[0] = (i + 0.5) * CellWidth[0][0];

        index = (k * GridDimension[1] + j) * GridDimension[0] + i;

        // Add MC tracers to this cell 
        for (n = 0; n < NumberOfMonteCarloTracerParticlesPerCell; n++) {
          newMC = new MonteCarloTracerParticle(this, index, groupID, GridLevel, Time, pos);
          InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticles[index], newMC);
          TotalNumberOfMonteCarloTracerParticles++;
        }
      }
    }
  }
  return SUCCESS;
}
