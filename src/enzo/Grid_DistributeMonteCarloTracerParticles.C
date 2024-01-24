/***********************************************************************
/
/  GRID CLASS (Distribute Monte Carlo Tracer Particles)
/
/  written by: Corey Brummel-Smith
/  date:       March, 2022
/
/  PURPOSE: Move Monte Carlo Tracer particles from cell 0 to their
/           correct position (cell). This function is only called 
/           after the particles have been transfered to cell 0 from 
/           another grid during communication.
/
/  RETURNS:
/    SUCCESS or FAIL
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

void InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticle * &Node, MonteCarloTracerParticle * &NewNode);
MonteCarloTracerParticle *PopMonteCarloTracerParticle(MonteCarloTracerParticle * &Node); 

int grid::DistributeMonteCarloTracerParticles()
{
  if (MyProcessorNumber == ProcessorNumber) {  

    int dim, index, i, j, k;
    int index_ijk[MAX_DIMENSION];
    float DomainWidth[MAX_DIMENSION], DomainWidthInv[MAX_DIMENSION];  
    MonteCarloTracerParticle *mctp, *MoveMCTP;
    int ActiveDim[MAX_DIMENSION];

    for (int dim = 0; dim < 3; dim++)
      ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] +1;      

    if (this->MonteCarloTracerParticles == NULL)
      this->AllocateMonteCarloTracerParticleData();
      //ENZO_FAIL("\nDistributeMonteCarloTracerParticles: Attempting to distribute Monte Carlo Tracer particles in a grid that has not been initialized with Monte Carlo Tracer particles.\n");
    
    mctp = this->MonteCarloTracerParticles[0];
    this->MonteCarloTracerParticles[0] = NULL;

    if (mctp == NULL) {
      printf("\nproc%d: DistributeMonteCarloTracerParticles: No particles to distribute.\n", MyProcessorNumber);
      fflush(stdout);
      return SUCCESS;
    }

    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
      DomainWidthInv[dim] = 1.0/DomainWidth[dim];    
    }

    while (mctp != NULL) {

      MoveMCTP = PopMonteCarloTracerParticle(mctp);  // also advances to NextParticle

      /* Find which cell this particle belongs in */
      // for (dim = 0; dim < GridRank; dim++) {
      //   index_ijk[dim] = (int) (ActiveDim[dim] * 
      //                           (MoveMCTP->Position[dim] - GridLeftEdge[dim]) /
      //                           (GridRightEdge[dim] - GridLeftEdge[dim])) 
      //                         + NumberOfGhostZones;

      //   if (index_ijk[dim] < GridStartIndex[dim]){
      //     // FOR DEBUGGING
      //     printf("index[%d] = %d, GridStartIndex[%d] = %d, MP-P %d-%d", dim, index_ijk[dim], dim, GridStartIndex[dim], MyProcessorNumber, ProcessorNumber);
      //     fflush(stdout);
      //   }

      // }
      //index = GetIndex(index_ijk[0], index_ijk[1], index_ijk[2]);
      i = int((MoveMCTP->Position[0] - GridLeftEdge[0]) / CellWidth[0][0]);
      j = int((MoveMCTP->Position[1] - GridLeftEdge[1]) / CellWidth[1][0]);
      k = int((MoveMCTP->Position[2] - GridLeftEdge[2]) / CellWidth[2][0]);
      index = GRIDINDEX(i,j,k);      
      InsertMonteCarloTracerParticleAfter(this->MonteCarloTracerParticles[index], MoveMCTP);
    }

    printf("\nproc%d: DistributeMonteCarloTracerParticles: Success.\n", MyProcessorNumber);
    fflush(stdout);
    }
  return SUCCESS;
}