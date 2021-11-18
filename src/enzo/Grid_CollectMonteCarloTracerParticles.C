/***********************************************************************
/
/  GRID CLASS (COLLECT PARTICLES INTO ONE PROCESSOR)
/
/  written by: John Wise
/  date:       May, 2009
/  modified1:  Corey Brummel-Smith, November, 2021 
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
		       int &MonteCarloTracerParticletIndex, int &EndIndex, 
		       mc_tracer_data* &List, int CopyDirection)
{
 
  /* Declarations. */

  int i, j, k, index, dim, icell, nmc, grid, proc;
  int size = 1, NumberOfMCTracersInGhostZones = 0;
  MonteCarloTracerParticle *MoveMCTP;
  NumberOfMCTracersInGhostZones = this->ReturnNumberOfMonteCarloTracerParticlesInFirstGhostCells();
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];  

  /* ----------------------------------------------------------------- */
  /* Copy star out of grid. */

  if (CopyDirection == COPY_OUT) {

    /* If this is the correct processor, no copy-outs required. */

    if (MyProcessorNumber == ProcessorNumber)
      return SUCCESS;

    /* If there are no particles in the ghost zones, we're done. */

    if (NumberOfMCTracersInGhostZones == 0)
      return SUCCESS;      

    /* Add to the MC particle count to move */

    // NumberOfMCTracersInGhostZones is still the number of local particles to COPY_OUT, not the
    // global total!
    NumberToMove[ProcessorNumber] += NumberOfMCTracersInGhostZones;
 
    /* Move MC tracers to buffer then delete from the grid of linked lists */

    // loop over all ghost cells neighboring the active region

    const int NUMBER_OF_FACE_CELLS = 2 * (GridDimension[0] * GridDimension[1] +
                                          GridDimension[0] * GridDimension[2] +
                                          GridDimension[1] * GridDimension[2]);

    int FirstGhostCellIndices = int[NUMBER_OF_FACE_CELLS];
    FirstGhostCellIndices = this->GetFirstGhostCells(NUMBER_OF_FACE_CELLS);

    for (icell = 0; icell < NUMBER_OF_FACE_CELLS; icell++) {
      index = FirstGhostCellIndices[icell];

      // Unravel index
      k = index / (GridDimension[0] * GridDimension[1]);
      j = (index / GridDimension[0]) % GridDimension[1];
      i = index % GridDimension[0];

      // Compute particle position (cell-center)
      pos[2] = (k + 0.5) * CellWidth[2][0];
      pos[1] = (j + 0.5) * CellWidth[1][0];
      pos[0] = (i + 0.5) * CellWidth[0][0];

      for (nmc = StartIndex, 
           MoveMCTP = this->MonteCarloTracerParticles[index]; 
           MoveMCTP != NULL; 
           MoveMCTP = MoveMCTP->NextParticle,
           nmc++;) {

        MoveMCTP->MonteCarloTracerParticleToBuffer(&List[nmc].data, pos);
        List[nmc].grid = GridNum;
        List[nmc].proc = ProcessorNumber;

      } // End loop over this cell's particles 

      StartIndex = nmc; // I don't think this is needed because we're always looping over all particles
      DeleteMonteCarloTracerParticleList(MonteCarloTracerParticles[index]);      

    } // End loop over ghost cells

  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy stars back into grid. */
 
  else {

    if (MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    /* Count up total number. */
    /* NumberOfMCTracersInGhostZones should be zero because if we're 
        using COPY_IN mode, All ghost zone particles should have been
        shared with other processors and deleted from this grid when 
        COPY_OUT mode was previously called. */
    if NumberOfMCTracersInGhostZones
 
    int TotalNumberOfMCTracersInGhostZones;
    int NumberOfNewMonteCarloTracerParticles = EndIndex - StartIndex;

    TotalNumberOfMCTracersInGhostZones = NumberOfMCTracersInGhostZones + NumberOfNewMonteCarloTracerParticles;

    // ***********************************************************
    // THIS NEEDS TO BE CHANGED BECAUSE WE NEED TO LOOP OVER CELLS AND PLACE PARTICLES
    // IN THE CORRECT CELL
    if (NumberOfNewMonteCarloTracerParticles > 0)
      for (i = StartIndex; i < EndIndex; i++) {
        
        /* Find if this particle belongs in this grid */
        


  	    MoveMCTP = MonteCarloTracerParticleBufferToList(List[i].data);
  	    MoveMCTP->CurrentGrid = this;
  	    InsertMonteCarloTracerParticleAfter(this->MonteCarloTracerParticles[index], MoveMCTP);
      }
    // ***********************************************************

 
    /* Set new number of stars in this grid. */
    // ****!!!!! THIS IS UNNECESSARY BECAUSE WE DON'T KEEP/NEED A RECORD OF THE NUMBER OF MC TRACERS *****/
    NumberOfMCTracersInGhostZones = TotalNumberOfMCTracersInGhostZones; 
 
  } // end: if (COPY_IN)
 
  return SUCCESS;
}
