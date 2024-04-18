/***********************************************************************
/
/  COMMUNICATION ROUTINE: TRANSFER MONTE CARLO TRACER PARTICLES
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified:   May, 2009 by John Wise: optimized version to transfer
/                particles in one sweep with collective calls.
/  modified:   July, 2009 by John Wise: adapted for stars
/  modified:   December, 2021 by Corey Brummel-Smith: adapted for 
                 Monte Carlo tracer particles
/
/  PURPOSE:
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
#include "SortCompareFunctions.h"
void my_exit(int status);
 
// function prototypes
 
int Enzo_Dims_create(int nnodes, int ndims, int *dims);
int search_lower_bound(int *arr, int value, int low, int high, 
		       int total);

// Remove define.  This method will always be used.
//#define KEEP_PARTICLES_LOCAL
 
int CommunicationTransferMonteCarloTracerParticles(grid *GridPointer[], int NumberOfGrids,
			       int TopGridDims[])
{

  if (NumberOfGrids == 1)
    return SUCCESS;

  int i, j, jstart, jend, dim, grid, proc, DisplacementCount, ThisCount;
  float ExactDims, ExactCount;

  /* Assume that the grid was split by Enzo_Dims_create, and create a
     map from grid number to an index that is determined from (i,j,k)
     of the grid partitions. */

  int *GridMap = new int[NumberOfGrids];
  int *StartIndex[MAX_DIMENSION];
  int Layout[MAX_DIMENSION], LayoutTemp[MAX_DIMENSION];
  int Rank, grid_num, bin, CenterIndex;
  int *pbin;
  int GridPosition[MAX_DIMENSION], Dims[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION];

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    Layout[dim] = 0;
    GridPosition[dim] = 0;
    StartIndex[dim] = NULL;
  }

  GridPointer[0]->ReturnGridInfo(&Rank, Dims, Left, Right); // need rank
  Enzo_Dims_create(NumberOfGrids, Rank, LayoutTemp);

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Layout[dim] = 0;
  for (dim = 0; dim < Rank; dim++)
    Layout[Rank-1-dim] = LayoutTemp[dim];

  /* For unequal splits of the topgrid, we need the start indices of
     the partitions.  It's easier to recalculate than to search for
     them. */

  for (dim = 0; dim < Rank; dim++) {

    StartIndex[dim] = new int[Layout[dim]+1];
    ExactDims = float(TopGridDims[dim]) / float(Layout[dim]);
    ExactCount = 0.0;
    DisplacementCount = 0;

    for (i = 0; i < Layout[dim]; i++) {
      ExactCount += ExactDims;
      if (dim == 0)
	ThisCount = nint(0.5*ExactCount)*2 - DisplacementCount;
      else
	ThisCount = nint(ExactCount) - DisplacementCount;
      StartIndex[dim][i] = DisplacementCount;
      DisplacementCount += ThisCount;
    } // ENDFOR i    

    StartIndex[dim][Layout[dim]] = TopGridDims[dim];

  } // ENDFOR dim

  for (grid = 0; grid < NumberOfGrids; grid++) {
    GridPointer[grid]->ReturnGridInfo(&Rank, Dims, Left, Right);
    for (dim = 0; dim < Rank; dim++) {

      if (Layout[dim] == 1) {
	GridPosition[dim] = 0;
      } else {

	CenterIndex = 
	  int(TopGridDims[dim] *
	      (0.5*(Right[dim]+Left[dim]) - DomainLeftEdge[dim]) /
	      (DomainRightEdge[dim] - DomainLeftEdge[dim]));

	GridPosition[dim] = 
	  search_lower_bound(StartIndex[dim], CenterIndex, 0, Layout[dim],
			     Layout[dim]);

      } // ENDELSE

    } // ENDFOR dim
    grid_num = GridPosition[0] + 
      Layout[0] * (GridPosition[1] + Layout[1]*GridPosition[2]);
    GridMap[grid_num] = grid;
  } // ENDFOR grids
 
  int *NumberToMove = new int[NumberOfProcessors];
  mc_tracer_data *SendList = NULL;
  mc_tracer_data *SharedList = NULL;
 
  for (i = 0; i < NumberOfProcessors; i++)
    NumberToMove[i] = 0;
 
  /* Generate the list of Monte Carlo tracer particle moves. */
  printf("\nComTrans");
  int Zero = 0;
  for (grid = 0; grid < NumberOfGrids; grid++)
    GridPointer[grid]->CommunicationTransferMonteCarloTracerParticles
      (GridPointer, NumberOfGrids, grid, TopGridDims, NumberToMove, Zero, Zero, 
       SendList, Layout, StartIndex, GridMap, COPY_OUT);

  int TotalNumberToMove = 0;
  for (proc = 0; proc < NumberOfProcessors; proc++)
    TotalNumberToMove += NumberToMove[proc];

  int NumberOfReceives = 0;

  /* Sort by grid number.  We no longer share with other processors
     because the Monte Carlo Tracer Particles are stored locally until
     CommunicationCollectParticles(SIBLINGS_ONLY). */

  SharedList = SendList;
  NumberOfReceives = TotalNumberToMove;
  int mc_tracer_data_size = sizeof(mc_tracer_data);
  std::sort(SharedList, SharedList+TotalNumberToMove, cmp_mc_tracer_grid());

  /* Copy Monte Carlo tracer particles back to grids */

  jstart = 0;
  jend = 0;

  // Copy shared Monte Carlo tracer particles to grids, if any
  if (NumberOfReceives > 0) {
    for (j = 0; j < NumberOfGrids && jend < NumberOfReceives; j++) {
      while (SharedList[jend].grid <= j) {
    	 jend++;
    	 if (jend == NumberOfReceives) break;
      }

      GridPointer[j]->CommunicationTransferMonteCarloTracerParticles
	      (GridPointer, NumberOfGrids, j, TopGridDims, NumberToMove, 
	      jstart, jend, SharedList, Layout, StartIndex, GridMap, COPY_IN);

      jstart = jend;
    } // ENDFOR grids
  } // ENDIF NumberOfRecieves > 0

  /* Cleanup. */

  if (SendList != SharedList)
    delete [] SendList;
  delete [] SharedList;
  delete [] NumberToMove;
  delete [] GridMap;
  for (dim = 0; dim < Rank; dim++)
    delete [] StartIndex[dim];

  CommunicationSumValues(&TotalNumberToMove, 1);
  if (debug)
    printf("CommunicationTransferMonteCarloTracerParticles: moved = %"ISYM"\n",
  	   TotalNumberToMove);

  return SUCCESS;
}

