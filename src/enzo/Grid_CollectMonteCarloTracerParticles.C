/***********************************************************************
/
/  GRID CLASS (COLLECT PARTICLES INTO ONE PROCESSOR)
/
/  written by: John Wise
/  date:       May, 2009
/  modified1:  Corey Brummel-Smith, November, 2024 
/              ( Adapted for Monte Carlo Tracer Particles (MCTP) )
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

void DeleteMonteCarloTracerParticleList(MonteCarloTracerParticle * &Node);
MonteCarloTracerParticle* MonteCarloTracerParticleBufferToList(MonteCarloTracerParticleBuffer buffer);
void InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticle * &Node, MonteCarloTracerParticle * &NewNode);
 
int grid::CollectMonteCarloTracerParticles(int GridNum, int* &NumberToMove, 
		       int &StartIndex, int &EndIndex, 
		       mc_tracer_data* &List, int CopyDirection)
{
 
  /* Declarations. */

  int i, i0, j0, k0, index, dim, icell, nmc, grid, proc;
  int size = 1;
  MonteCarloTracerParticle *MoveMCTP;
  FLOAT pos[3];

  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];  

  /* ----------------------------------------------------------------- */
  /* Copy star out of grid. */

  if (CopyDirection == COPY_OUT) {

    /* If this is the correct processor, no copy-outs required. */

    if (MyProcessorNumber == ProcessorNumber)
      return SUCCESS;

    /* If there are no particles, we're done. */

    if (NumberOfMonteCarloTracerParticles == 0)
      return SUCCESS;      

    /* Add to the MC particle count to move */

    // NumberOfMonteCarloTracerParticles is still the number of local particles to COPY_OUT, not the
    // global total!
    NumberToMove[ProcessorNumber] += NumberOfMonteCarloTracerParticles;
 
    /* Move MC tracers to buffer then delete from the grid of linked lists */
    nmc = StartIndex;
    for (k0 = GridStartIndex[2]; k0 <= GridEndIndex[2]; k0++) {
      for (j0 = GridStartIndex[1]; j0 <= GridEndIndex[1]; j0++) {
        for (i0 = GridStartIndex[0]; i0 <= GridEndIndex[0]; i0++) {
 
          index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;

          // Compute particle position (cell-center)
          pos[2] = CellLeftEdge[2][k0] + 0.5 * CellWidth[2][0];
          pos[1] = CellLeftEdge[1][j0] + 0.5 * CellWidth[1][0];
          pos[0] = CellLeftEdge[0][i0] + 0.5 * CellWidth[0][0];    


          for (MoveMCTP = this->MonteCarloTracerParticles[index]; 
               MoveMCTP != NULL; 
               MoveMCTP = MoveMCTP->NextParticle,
               nmc++) {

            MoveMCTP->MonteCarloTracerParticleToBuffer(&List[nmc].data, pos);
            List[nmc].grid = GridNum;
            List[nmc].proc = ProcessorNumber;

          } // End loop over this cell's particles 

          DeleteMonteCarloTracerParticleList(MonteCarloTracerParticles[index]);      

        } // End for i0
      } // End for j0
    } // End for k0
    StartIndex = nmc; // I don't think this is needed because we're always looping over all particles

  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy particles back into grid. */
 
  else {

    if (MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    /* Count up total number. */

    printf("%sGrid_CollectMonteCarloTracerParticles: COPY_IN\nNumberOfMCTracers = %d (0?)", NumberOfMonteCarloTracerParticles);
 
    int TotalNumberOfMCTracers;
    int NumberOfNewMonteCarloTracerParticles = EndIndex - StartIndex;

    TotalNumberOfMCTracers = NumberOfMonteCarloTracerParticles + NumberOfNewMonteCarloTracerParticles;

    if (NumberOfNewMonteCarloTracerParticles > 0)
      for (i = StartIndex; i < EndIndex; i++) {
        
        /* Add particles to this grid. We don't yet know if 
           these particles belong in this grid so it doesn't matter where we
           put them. We arbitrarily put them in the first cell for now.
           NOTE: CurrentGrid may not be correct. This will need
           to be updated in CommunicationTransferParticles. I'M NOT SURE THIS IS TRUE. I THINK
           I SHOULD DISTRIBUTE BELOW*/
        
  	    MoveMCTP = MonteCarloTracerParticleBufferToList(List[i].data);
  	    MoveMCTP->CurrentGrid = this;
  	    InsertMonteCarloTracerParticleAfter(this->MonteCarloTracerParticles[0], MoveMCTP);
      }
      /* Move particles to their correct cell */
      DistributeMonteCarloTracerParticles();

      /* Set new number of stars in this grid. */
      NumberOfMonteCarloTracerParticles = TotalNumberOfMonteCarloTracerParticles;
 

  } // end: if (COPY_IN)
 
  return SUCCESS;
}
