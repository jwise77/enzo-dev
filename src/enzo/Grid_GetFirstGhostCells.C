/***********************************************************************
/
/  GRID CLASS (GET INDICIES OF 1st GHOST CELLS)
/
/  written by: Corey Brummel-Smith
/  date:       Noveber, 2021
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

int* grid::GetFirstGhostCells(const int NUMBER_OF_FACE_CELLS)
{
    int GhostCellIndices = int[NUMBER_OF_FACE_CELLS];
    int n = 0;

    /* x inner (left) face */
     
    i = StartIndex[0] - 1; // first ghost zone to the left of the active zone

    for (j = 0; j < GridDims[1]; j++)
    for (k = 0; k < GridDims[2]; k++) {

      GhostCellIndices[n] = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      n++;
    }

    /* x outer (right) face */
    i = EndIndex[0] + 1;
    for (j = 0; j < GridDims[1]; j++)
    for (k = 0; k < GridDims[2]; k++) {

      GhostCellIndices[n] = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      n++;
    }

    /* y inner (left) face */
     
    j = StartIndex[1] - 1;
    for (i = 0; i < GridDims[0]; i++)
    for (k = 0; k < GridDims[2]; k++) {

      GhostCellIndices[n] = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      n++;
    }

    /* y outer (right) face */

    j = EndIndex[1] + 1;
    for (i = 0; i < GridDims[0]; i++)
    for (k = 0; k < GridDims[2]; k++) {

      GhostCellIndices[n] = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      n++;
    }

    /* z inner (left) face */
     
    k = StartIndex[2] - 1;
    for (i = 0; i < GridDims[0]; i++)
    for (j = 0; j < GridDims[1]; j++) {

      GhostCellIndices[n] = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      n++;
    }

    /* z outer (right) face */

    k = EndIndex[2] + 1;
    for (i = 0; i < GridDims[0]; i++)
    for (j = 0; j < GridDims[1]; j++) {

      GhostCellIndices[n] = i + j*GridDims[0] + k*GridDims[1]*GridDims[0];
      n++;
    }

    return GhostCellIndices;

}
