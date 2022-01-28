/***********************************************************************
/
/  GRID CLASS (COPY STARS INTO OR OUT OF GRID)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:  Robert Harkness
/  date:       April, 2006
/  modified2:  May, 2009 by John Wise: optimized version to transfer
/                particles in one sweep with collective calls.
/  modified3:  July, 2009 by John Wise: adapted for stars
/  modified4:  December, 2021 by Corey Brummel-Smith: adapted for 
                 Monte Carlo tracer particles
/
/  PURPOSE:
/
************************************************************************/
#ifdef USE_MPI
#include <mpi.h>
#endif
#include <stdio.h>
#include <string.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"

MonteCarloTracerParticle *PopMonteCarloTracerParticle(MonteCarloTracerParticle * &Node);
MonteCarloTracerParticle* MonteCarloTracerParticleBufferToList(MonteCarloTracerParticleBuffer buffer);
void InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticle * &Node, MonteCarloTracerParticle * &NewNode);
int search_lower_bound(int *arr, int value, int low, int high, 
		       int total);
 
int grid::CommunicationTransferMonteCarloTracerParticles(grid* Grids[], int NumberOfGrids,
				     int ThisGridNum, int TopGridDims[],
				     int *&NumberToMove, 
				     int StartIndex, int EndIndex, 
				     mc_tracer_data *&List, int *Layout, 
				     int *GStartIndex[], int *GridMap, 
				     int CopyDirection)
{
 
  /* Declarations. */
 
  int i, j, k, dim, grid, proc, grid_num, width, bin, CenterIndex, index,
      NumberOfMCTPInCellZero;
  int GridPosition[MAX_DIMENSION], index_ijk[MAX_DIMENSION];
  float DomainWidth[MAX_DIMENSION], DomainWidthInv[MAX_DIMENSION];  
  FLOAT r[MAX_DIMENSION];
  int *ToGrid, *pbin;
  MonteCarloTracerParticle *mctp, *MoveMCTP;

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    GridPosition[dim] = 0;
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    DomainWidthInv[dim] = 1.0/DomainWidth[dim];    
  }
 
  /* ----------------------------------------------------------------- */
  /* Copy Monte Carlo tracer particles out of grid. */
 
  if (CopyDirection == COPY_OUT) {

    /* Count the number of Monte Carlo Tracer particles to move. 
       All the MC tracers to move are now stored in the grid cell 0.
       This was done in 
       CommunicationCollectMonteCarloTracerParticles(COPY_IN) */

    NumberOfMCTPInCellZero = this->CountMonteCarloTracerParticlesInCellZero();

    /* If there are no Monte Carlo tracer particles to move, we're done. */

    if (NumberOfMCTPInCellZero == 0)
      return SUCCESS;
 
    /* Count the number of Monte Carlo tracer particles already moved */

    int PreviousTotalToMove = 0;
    for (i = 0; i < NumberOfProcessors; i++)
      PreviousTotalToMove += NumberToMove[i];

    /* Count Monte Carlo tracer particles to move.  
       Apply perioidic wrap to the Monte Carlo tracer particles. */
 
    ToGrid = new int[NumberOfMCTPInCellZero];

    // Periodic boundaries
    for (dim = 0; dim < GridRank; dim++) 
      for (mctp = MonteCarloTracerParticles[0]; mctp; mctp = mctp->NextParticle) {
      	if (mctp->Position[dim] > DomainRightEdge[dim])
      	  mctp->Position[dim] -= DomainWidth[dim];
      	else if (mctp->Position[dim] < DomainLeftEdge[dim])
      	  mctp->Position[dim] += DomainWidth[dim];
      }

    for (mctp = MonteCarloTracerParticles[0], i = 0; 
         mctp; 
         mctp = mctp->NextParticle, i++) {

      for (dim = 0; dim < GridRank; dim++) {

      	if (Layout[dim] == 1) {
      	  GridPosition[dim] = 0;
      	} 
        else {
      	  CenterIndex = 
      	  (int) (TopGridDims[dim] * 
      		 (mctp->Position[dim] - DomainLeftEdge[dim]) *
      		 DomainWidthInv[dim]);

      	  GridPosition[dim] = 
      	    search_lower_bound(GStartIndex[dim], CenterIndex, 0, Layout[dim],
      			       Layout[dim]);
      	  GridPosition[dim] = min(GridPosition[dim], Layout[dim]-1);

      	} // ENDELSE Layout

      } // ENDFOR dim

      grid_num = GridPosition[0] + Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
      grid = GridMap[grid_num];
      if (grid != ThisGridNum) {
      	proc = Grids[grid]->ReturnProcessorNumber();
      	NumberToMove[proc]++;
      }
      ToGrid[i] = grid;

    } // ENDFOR Monte Carlo tracer particles

    /* Allocate space. */
 
    int TotalToMove = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++)
      TotalToMove += NumberToMove[proc];

    if (TotalToMove > PreviousTotalToMove) {
 
      /* Move Monte Carlo tracer particles into list */

      // Increase the size of the list to include the Monte Carlo tracer particles from this grid

      mc_tracer_data *NewList = new mc_tracer_data[TotalToMove];
      memcpy(NewList, List, PreviousTotalToMove * sizeof(mc_tracer_data));
      delete [] List;
      List = NewList;

      int n1 = PreviousTotalToMove;

      mctp = MonteCarloTracerParticles[0];
      MonteCarloTracerParticles[0] = NULL;

      i = 0;
      while (mctp != NULL) {

	      MoveMCTP = PopMonteCarloTracerParticle(mctp);  // also advances to NextParticle
	      grid = ToGrid[i];

      	if (grid != ThisGridNum) {
          // MoveMCTP->Position was previously set in Grid_CollectMonteCarloTracerParticles
          MoveMCTP->MonteCarloTracerParticleToBuffer(&List[n1].data, MoveMCTP->Position);
      	  List[n1].grid = grid;
      	  List[n1].proc = MyProcessorNumber;
      	  n1++;
      	  delete MoveMCTP;
      	}  // ENDIF different processor

      	// Particle already in this grid (Only move from cell 0 to the correct cell)
      	else {
          for (dim = 0; dim < GridRank; dim++) {
              index_ijk[dim] = (int) (GridDimension[dim] * 
                                      (MoveMCTP->Position[dim] - DomainLeftEdge[dim]) *
                                      DomainWidthInv[dim]);
          }
          index = GetIndex(index_ijk[0], index_ijk[1], index_ijk[2]);
      	  InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticles[index], MoveMCTP);
      	}
      	i++;

      } // ENDWHILE Monte Carlo tracer particles
      
    } // ENDIF TotalToMove > PreviousTotalToMove

    delete [] ToGrid;
 
  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy Monte Carlo tracer particles back into grid. */
 
  else {
 
    int NumberOfNewMonteCarloTracerParticles = EndIndex - StartIndex;

    /* Compute which cell the particle needs to be deposited into. 
       Copy Monte Carlo tracer particles from buffer into the 
       the linked list for that cell */

    if (NumberOfNewMonteCarloTracerParticles > 0)

      for (i = StartIndex; i < EndIndex; i++) {

      	MoveMCTP = MonteCarloTracerParticleBufferToList(List[i].data);
      	MoveMCTP->CurrentGrid = this;

        for (dim = 0; dim < GridRank; dim++) {
            index_ijk[dim] = (int) (GridDimension[dim] * 
                                    (MoveMCTP->Position[dim] - DomainLeftEdge[dim]) *
                                    DomainWidthInv[dim]);
        }
        index = GetIndex(index_ijk[0], index_ijk[1], index_ijk[2]);
      	InsertMonteCarloTracerParticleAfter(this->MonteCarloTracerParticles[index], MoveMCTP);
      } // ENDFOR Monte Carlo tracer particles

  } // end: if (COPY_IN)
 
  return SUCCESS;
}
