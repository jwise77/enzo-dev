/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (SET A GRID'S BOUNDARY)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       November, 2005
/              Out-of-core handling for the boundary
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <hdf5.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

#define USE_PERIODIC

int Move_MonteCarloTracerParticles_From_CellA_to_CellB(MonteCarloTracerParticle *&headA, MonteCarloTracerParticle *&headB);

int ExternalBoundary::SetExternalBoundaryMonteCarloTracerParticles(int GridRank,
                int GridDims[], int GridOffset[], int StartIndex[], int EndIndex[],
                grid *GridData)
{
  /* declarations */
 
  int i, j, k, dim, index_ghost, bindex, index_periodic_active;
 
  /* error check: grid ranks */
 
  if (GridRank != BoundaryRank) {
    ENZO_VFAIL("GridRank(%"ISYM") != BoundaryRank(%"ISYM").\n",
            GridRank, BoundaryRank)
  }
 
  /* set Boundary conditions */

  // call is by Field - set all 6 faces
  // dim = 0, face = 0
  // dim = 0, face = 1
  // dim = 1, face = 0
  // dim = 1, face = 1
  // dim = 2, face = 0
  // dim = 2, face = 1
 
  if (BoundaryDimension[0] > 1 && GridOffset[0] == 0) {
 
    /* set x inner (left) face */
 
    i = StartIndex[0] - 1; // first ghost zone to the left of the active zone
    for (j = 0; j < GridDims[1]; j++)
    for (k = 0; k < GridDims[2]; k++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      bindex = j+GridOffset[1] + (k+GridOffset[2])*BoundaryDimension[1];

      switch (BoundaryType[0][0][0][bindex]) {
      case reflecting:
        ENZO_FAIL("Reflecting boundary conditions not implemented yet!\n");
        break;
      case outflow:
        ENZO_FAIL("Outflow boundary conditions not implemented yet!\n");
        break;
      case inflow:
        ENZO_FAIL("Inflow boundary conditions not implemented.\n");
        break;
      case periodic:
#ifdef USE_PERIODIC
        // TODO
        index_periodic_active = index_ghost + (EndIndex[0] - StartIndex[0] + 1);
        Move_MonteCarloTracerParticles_From_CellA_to_CellB(
          GridData->MonteCarloTracerParticles[index_ghost], 
          GridData->MonteCarloTracerParticles[index_periodic_active]);
#endif /* USE_PERIODIC */
        break;
      case shearing:
        ENZO_FAIL("Shearing boundary conditions not implemented yet.\n");
        break;
      case BoundaryUndefined:
            break;
      default:
        ENZO_VFAIL("BoundaryType %"ISYM" not recognized (x-left).\n",
            BoundaryType[0][0][0][bindex])
      }
    }
  }
 
  if (BoundaryDimension[0] > 1 && GridOffset[0]+GridDims[0] == BoundaryDimension[0]) {
   
    /* set x outer (right) face */
    i = EndIndex[0] + 1;
    for (j = 0; j < GridDims[1]; j++)
    for (k = 0; k < GridDims[2]; k++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      bindex = j+GridOffset[1] + (k+GridOffset[2])*BoundaryDimension[1];

      switch (BoundaryType[0][0][1][bindex]) {
      case reflecting:
        ENZO_FAIL("Reflecting boundary conditions not implemented yet!\n");
        break;
      case outflow:
        ENZO_FAIL("Outflow boundary conditions not implemented yet!\n");
        break;
      case inflow:
        ENZO_FAIL("Inflow boundary conditions not implemented.\n");
        break;
      case periodic:
#ifdef USE_PERIODIC
        // TODO
        index_periodic_active = index_ghost - (EndIndex[0] - StartIndex[0] + 1);
        Move_MonteCarloTracerParticles_From_CellA_to_CellB(
          GridData->MonteCarloTracerParticles[index_ghost], 
          GridData->MonteCarloTracerParticles[index_periodic_active]);

#endif /* USE_PERIODIC */
        break;
      case shearing:
        ENZO_FAIL("Shearing boundary conditions not implemented yet.\n");
        break;
      case BoundaryUndefined:
            break;
      default:
        ENZO_VFAIL("BoundaryType %"ISYM" not recognized (x-right).\n",
            BoundaryType[0][0][1][bindex])
      }
    }
  }
  /* set y inner (left) face */
 
  if (BoundaryDimension[1] > 1 && GridOffset[1] == 0) {

    j = StartIndex[1] - 1;
    for (i = 0; i < GridDims[0]; i++)
    for (k = 0; k < GridDims[2]; k++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      bindex = i+GridOffset[0] + (k+GridOffset[2])*BoundaryDimension[0];

      switch (BoundaryType[0][1][0][bindex]) {
      case reflecting:
        ENZO_FAIL("Reflecting boundary conditions not implemented yet!\n");
        break;
      case outflow:
        ENZO_FAIL("Outflow boundary conditions not implemented yet!\n");
        break;
      case inflow:
        ENZO_FAIL("Inflow boundary conditions not implemented.\n");
         break;
      case periodic:
#ifdef USE_PERIODIC
        // TODO
        index_periodic_active = index_ghost + (EndIndex[1] - StartIndex[1] + 1)*GridDims[0];
        Move_MonteCarloTracerParticles_From_CellA_to_CellB(
          GridData->MonteCarloTracerParticles[index_ghost], 
          GridData->MonteCarloTracerParticles[index_periodic_active]);

#endif /* USE_PERIODIC */
         break;
      case shearing:
        ENZO_FAIL("Shearing boundary conditions not implemented yet.\n");
        break;
      case BoundaryUndefined:
            break;
      default:
        ENZO_VFAIL("BoundaryType %"ISYM" not recognized (y-left).\n",
            BoundaryType[0][1][0][bindex])
      }
    }
  }
 
  if (BoundaryDimension[1] > 1 && GridOffset[1]+GridDims[1] == BoundaryDimension[1]) {
 
    /* set y outer (right) face */

    j = EndIndex[1] + 1;
    for (i = 0; i < GridDims[0]; i++)
    for (k = 0; k < GridDims[2]; k++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      bindex = i+GridOffset[0] + (k+GridOffset[2])*BoundaryDimension[0];

      switch (BoundaryType[0][1][1][bindex]) {
      case reflecting:
        ENZO_FAIL("Reflecting boundary conditions not implemented yet!\n");
        break;
      case outflow:
        ENZO_FAIL("Outflow boundary conditions not implemented yet!\n");
        break;
      case inflow:
        ENZO_FAIL("Inflow boundary conditions not implemented.\n");
        break;
      case periodic:
#ifdef USE_PERIODIC
        // TODO
        index_periodic_active = index_ghost - (EndIndex[1] - StartIndex[1] + 1)*GridDims[0];
        Move_MonteCarloTracerParticles_From_CellA_to_CellB(
          GridData->MonteCarloTracerParticles[index_ghost], 
          GridData->MonteCarloTracerParticles[index_periodic_active]);

#endif /* USE_PERIODIC */
        break;
      case shearing:
        ENZO_FAIL("Shearing boundary conditions not implemented yet.\n");
        break;
      case BoundaryUndefined:
            break;
      default:
        ENZO_VFAIL("BoundaryType %"ISYM" not recognized (y-right).\n",
            BoundaryType[0][1][1][bindex])
      }
    }
  }
 
  /* set z inner (left) face */
 
  if (BoundaryDimension[2] > 1 && GridOffset[2] == 0) {
 
    k = StartIndex[2] - 1;
    for (i = 0; i < GridDims[0]; i++)
    for (j = 0; j < GridDims[1]; j++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      bindex = i+GridOffset[0] + (j+GridOffset[1])*BoundaryDimension[0];

      switch (BoundaryType[0][2][0][bindex]) {
      case reflecting:
        ENZO_FAIL("Reflecting boundary conditions not implemented yet!\n");
        break;
      case outflow:
        ENZO_FAIL("Outflow boundary conditions not implemented yet!\n");
        break;
      case inflow:
        ENZO_FAIL("Inflow boundary conditions not implemented.\n");
        break;
      case periodic:
#ifdef USE_PERIODIC
        // TODO
        index_periodic_active = index_ghost + (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1];
        Move_MonteCarloTracerParticles_From_CellA_to_CellB(
          GridData->MonteCarloTracerParticles[index_ghost], 
          GridData->MonteCarloTracerParticles[index_periodic_active]);

#endif /* USE_PERIODIC */
        break;
      case shearing:
        ENZO_FAIL("Shearing boundary conditions not implemented yet.\n");
        break;
      case BoundaryUndefined:
            break;
      default:
        ENZO_VFAIL("BoundaryType %"ISYM" not recognized (z-left).\n",
            BoundaryType[0][2][0][bindex])
      }
    }
  }
 
  if (BoundaryDimension[2] > 1 && GridOffset[2]+GridDims[2] == BoundaryDimension[2]) {
 
    /* set z outer (right) face */
    k = EndIndex[2]+1;
    for (i = 0; i < GridDims[0]; i++)
    for (j = 0; j < GridDims[1]; j++) {

      index_ghost = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      bindex = i+GridOffset[0] + (j+GridOffset[1])*BoundaryDimension[0];

      switch (BoundaryType[0][2][1][bindex]) {
      case reflecting:
        ENZO_FAIL("Reflecting boundary conditions not implemented yet!\n");
        break;
      case outflow:
        ENZO_FAIL("Outflow boundary conditions not implemented yet!\n");
        break;
      case inflow:
        ENZO_FAIL("Inflow boundary conditions not implemented.\n");
        break;
      case periodic:
#ifdef USE_PERIODIC
        // TODO
        index_periodic_active = index_ghost - (EndIndex[2]-StartIndex[2]+1)*GridDims[0]*GridDims[1];
        Move_MonteCarloTracerParticles_From_CellA_to_CellB(
          GridData->MonteCarloTracerParticles[index_ghost], 
          GridData->MonteCarloTracerParticles[index_periodic_active]);

#endif /* USE_PERIODIC */
        break;
      case shearing:
        ENZO_FAIL("Shearing boundary conditions not implemented yet.\n");
        break;
      case BoundaryUndefined:
            break;
      default:
        fprintf(stderr, "BoundaryType %"ISYM" not recognized (z-right).\n",
            BoundaryType[0][2][1][bindex]);
        ENZO_FAIL("Unrecognized IO BoundaryType!\n");
      }
    }
  }

  return SUCCESS;
}