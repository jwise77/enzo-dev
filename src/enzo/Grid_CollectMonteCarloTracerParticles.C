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
		       star_data* &List, int CopyDirection)
{
 
  /* Declarations. */

  int i, j, k, index, dim, n1, grid, proc, NumberOfMCTracersInGhostZones;
  int size = 1;
  MonteCarloTracerParticle *MoveMCTP;
  NumberOfMCTracersInGhostZones = this->CountMCTracersInGhostZones();
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];  

  /* ----------------------------------------------------------------- */
  /* Copy star out of grid. */

  if (CopyDirection == COPY_OUT) {

    /* If there are no particles in the ghost zones, we're done. */

    if (NumberOfMCTracersInGhostZones == 0)
      return SUCCESS;

    /* If this is the correct processor, no copy-outs required. */

    if (MyProcessorNumber == ProcessorNumber)
      return SUCCESS;

    /* Add to the MC particle count to move */

    // NumberOfMCTracersInGhostZones is still the number of local stars, not the
    // actual total!
    NumberToMove[ProcessorNumber] += NumberOfMCTracersInGhostZones;
 
    /* Move MC tracers to buffer then delete from the grid of linked lists */

    // loop over all ghost cells neighboring the active region

    // **** save code for later (DOESNT BELONG HERE) ****
    pos[2] = (k + 0.5) * CellWidth[2][0];
    pos[1] = (j + 0.5) * CellWidth[1][0];
    pos[0] = (i + 0.5) * CellWidth[0][0];
    index = (k * GridDimension[1] + j) * GridDimension[0] + i;
    // **************************************************
    int NumberOfFaceCells =  // count cell faces
    int GhostCellIndices = int[2
    /* x inner (left) face */
     
    i = StartIndex[0] - 1; // first ghost zone to the left of the active zone

    for (j = 0; j < GridDims[1]; j++)
    for (k = 0; k < GridDims[2]; k++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
    }

    /* x outer (right) face */
    i = EndIndex[0] + 1;
    for (j = 0; j < GridDims[1]; j++)
    for (k = 0; k < GridDims[2]; k++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
    }

    /* y inner (left) face */
     
    j = StartIndex[1] - 1;
    for (i = 0; i < GridDims[0]; i++)
    for (k = 0; k < GridDims[2]; k++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
    }

    /* y outer (right) face */

    j = EndIndex[1] + 1;
    for (i = 0; i < GridDims[0]; i++)
    for (k = 0; k < GridDims[2]; k++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
    }

    /* z inner (left) face */
     
    k = StartIndex[2] - 1;
    for (i = 0; i < GridDims[0]; i++)
    for (j = 0; j < GridDims[1]; j++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
    }

    /* z outer (right) face */

    k = EndIndex[2] + 1;
    for (i = 0; i < GridDims[0]; i++)
    for (j = 0; j < GridDims[1]; j++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
    }

    // END LOOP OVER GHOST CELLS  

        

    ////////////////////////////////////////////////////////////////////////////////
    // CLONE OF STAR PARTICLE CODE
    for (i = 0, n1 = MonteCarloTracerParticletIndex, MoveMCTP = MonteCarloTracerParticles; i < NumberOfMCTracersInGhostZones;
         i++, n1++, MoveMCTP = MoveMCTP->NextMonteCarloTracerParticle) {
      MoveMCTP->MonteCarloTracerParticleToBuffer(&List[n1].data);
      List[n1].grid = GridNum;
      List[n1].proc = ProcessorNumber;
    } // ENDFOR stars

    MonteCarloTracerParticletIndex = n1;
    DeleteMonteCarloTracerParticleList(MonteCarloTracerParticles);
    NumberOfMCTracersInGhostZones = 0;
    ////////////////////////////////////////////////////////////////////////////////

  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy stars back into grid. */
 
  else {

    if (MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    /* Count up total number. */
 
    int TotalNumberOfMCTracersInGhostZones;
    int NumberOfNewMonteCarloTracerParticles = EndIndex - MonteCarloTracerParticletIndex;

    TotalNumberOfMCTracersInGhostZones = NumberOfMCTracersInGhostZones + NumberOfNewMonteCarloTracerParticles;

    if (NumberOfNewMonteCarloTracerParticles > 0)
      for (i = MonteCarloTracerParticletIndex; i < EndIndex; i++) {
	MoveMCTP = MonteCarloTracerParticleBufferToList(List[i].data);
	MoveMCTP->CurrentGrid = this;
	MoveMCTP->GridID = this->ID;
	InsertMonteCarloTracerParticleAfter(this->MonteCarloTracerParticles, MoveMCTP);
      }
 
    /* Set new number of stars in this grid. */
 
    NumberOfMCTracersInGhostZones = TotalNumberOfMCTracersInGhostZones;
 
  } // end: if (COPY_IN)
 
  return SUCCESS;
}
