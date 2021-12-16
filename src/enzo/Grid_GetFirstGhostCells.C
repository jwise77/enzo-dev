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

void grid::GetFirstGhostCells(int GhostCellIndices[], const int NUMBER_OF_FACE_CELLS)
{
    int i, j, k, n = 0;

    /* x inner (left) face */
     
    i = GridStartIndex[0] - 1; // first ghost zone to the left of the active zone

    for (j = 0; j < GridDimension[1]; j++)
    for (k = 0; k < GridDimension[2]; k++) {

      GhostCellIndices[n] = i + j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];
      n++;
    }

    /* x outer (right) face */
    i = GridEndIndex[0] + 1;
    for (j = 0; j < GridDimension[1]; j++)
    for (k = 0; k < GridDimension[2]; k++) {

      GhostCellIndices[n] = i + j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];
      n++;
    }

    /* y inner (left) face */
     
    j = GridStartIndex[1] - 1;
    for (i = 0; i < GridDimension[0]; i++)
    for (k = 0; k < GridDimension[2]; k++) {

      GhostCellIndices[n] = i + j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];
      n++;
    }

    /* y outer (right) face */

    j = GridEndIndex[1] + 1;
    for (i = 0; i < GridDimension[0]; i++)
    for (k = 0; k < GridDimension[2]; k++) {

      GhostCellIndices[n] = i + j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];
      n++;
    }

    /* z inner (left) face */
     
    k = GridStartIndex[2] - 1;
    for (i = 0; i < GridDimension[0]; i++)
    for (j = 0; j < GridDimension[1]; j++) {

      GhostCellIndices[n] = i + j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];
      n++;
    }

    /* z outer (right) face */

    k = GridEndIndex[2] + 1;
    for (i = 0; i < GridDimension[0]; i++)
    for (j = 0; j < GridDimension[1]; j++) {

      GhostCellIndices[n] = i + j*GridDimension[0] + k*GridDimension[1]*GridDimension[0];
      n++;
    }

    return;

}
