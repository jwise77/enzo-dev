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

  printf("\n%s\n", "CreateMonteCarloTracerParticles...");

  this->AllocateMonteCarloTracerParticleData();

  printf("\n%s\n", "AllocatedMonteCarloTracerParticles...");  

  // Loop over all active cells 
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    pos[2] = CellLeftEdge[2][k] + 0.5 * CellWidth[2][0];
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      pos[1] = CellLeftEdge[1][j] + 0.5 * CellWidth[1][0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
        pos[0] = CellLeftEdge[0][i] + 0.5 * CellWidth[0][0]; 

        /* DEBUG */
        if (i != GridStartIndex[0] || j != GridStartIndex[1] || k != GridStartIndex[2])
          continue;

        // if (i != GridEndIndex[0] || j != GridEndIndex[1] || k != GridEndIndex[2])
        //   continue;

        // if (i % 2 != 0 || j % 2 != 0 || k % 2 != 0)
        //   continue;

        /* END DEBUG */

        index = (k * GridDimension[1] + j) * GridDimension[0] + i;

        // Add MC tracers to this cell 
        for (n = 0; n < NumberOfMonteCarloTracerParticlesPerCell; n++) {
          //printf("\n(i,j,k,n) = (%d,%d,%d,%d) New Particle", i, j, k, n);  
          newMC = new MonteCarloTracerParticle(this, index, groupID, GridLevel, Time, pos);
          //printf("\nnewMC %p", newMC);            
          InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticles[index], newMC);
          //printf("\nInserted %p", newMC);            
          TotalNumberOfMonteCarloTracerParticles++;
        }
      }
    }
  }
  printf("\nTotalNumberOfMonteCarloTracerParticles %i\n", TotalNumberOfMonteCarloTracerParticles);
  printf("%s\n", "CreateMonteCarloTracerParticles Done.");

  return SUCCESS;
}
