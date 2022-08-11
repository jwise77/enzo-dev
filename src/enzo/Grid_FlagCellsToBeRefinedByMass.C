#define NO_DEBUG_MRP
/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY MASS)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
int grid::FlagCellsToBeRefinedByMass(int level, int method, int RestrictFlag)
{
 
  /* Return if this grid is not on this processor. */
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* error check */

  int ThisFlaggingMethod = CellFlaggingMethod[method];
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  if (MassFlaggingField == NULL && ThisFlaggingMethod == 2) {
    fprintf(stderr, "Mass Flagging Field is undefined.\n");
    return -1;
  }
 
  if (ParticleMassFlaggingField == NULL && ThisFlaggingMethod == 4) {
    fprintf(stderr, "Mass Flagging Field is undefined.\n");
    return -1;
  }
 
  /* compute size */
 
  int i, size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Compute the ModifiedMinimumMass */
 
  float ModifiedMinimumMassForRefinement =
    MinimumMassForRefinement[method]*POW(RefineBy,
		    level*MinimumMassForRefinementLevelExponent[method]);
  if (ProblemType == 28)
    ModifiedMinimumMassForRefinement = 0;
 
  /* Flag points */

  float *ffield;
  if (ThisFlaggingMethod == 2)
    ffield = MassFlaggingField;
  else if (ThisFlaggingMethod == 4)
    ffield = ParticleMassFlaggingField;
  else {
    ENZO_VFAIL("Unrecognized mass refinement flagging method (%"ISYM")\n", 
	    method)
  }

  /* If RestrictFlag is true, then cells will only be flagged if they
     are flagged by mass.  This happens when using must-refine
     particles and must happen after all of the other flagging
     methods. */

  int nflag = 0, nunflag = 0, nonzero = 0;
  if (RestrictFlag) {
    for (i = 0; i < size; i++) {
      if (FlaggingField[i] == 0 && ffield[i] > ModifiedMinimumMassForRefinement) {
        nunflag++;
      }
      if (FlaggingField[i] > 0 && ffield[i] > ModifiedMinimumMassForRefinement) {
        nflag++;
      }
      FlaggingField[i] = ((ffield[i] > ModifiedMinimumMassForRefinement) &&
			  (FlaggingField[i] > 0)) ? 1 : 0;
    }
  } else {
    for (i = 0; i < size; i++) {
      FlaggingField[i] += (ffield[i] > ModifiedMinimumMassForRefinement) ? 1 : 0;
      if (ffield[i] > tiny_number) {
        nonzero++;
      }
#ifdef DEBUG_MRP
      if (ffield[i] > ModifiedMinimumMassForRefinement && level == 3) {
        int ii, jj, kk;
        FLOAT cx, cy, cz;
        ii = i % GridDimension[0];
        kk = i / (GridDimension[0] * GridDimension[1]);
        jj = (i - kk * GridDimension[0] * GridDimension[1]) / GridDimension[0];
        cx = CellLeftEdge[0][ii] + 0.5 * CellWidth[0][ii];
        cy = CellLeftEdge[1][jj] + 0.5 * CellWidth[1][jj];
        cz = CellLeftEdge[2][kk] + 0.5 * CellWidth[2][kk];
        printf("FlagByMass[G%"ISYM"]: Cell %d flagged, position = %lf %lf %lf, rho = %g\n", this->ID, i, cx, cy, cz, ffield[i]);
      }
#endif
    }
  }

  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++)
    if (FlaggingField[i] > 0)
      NumberOfFlaggedCells++;

  // Print MRP bounding box
#ifdef DEBUG_MRP
  if (NumberOfFlaggedCells > 0 && level == 3) {
    FLOAT minr[3] = {1, 1, 1};
    FLOAT maxr[3] = {0, 0, 0};
    for (i = 0; i < NumberOfParticles; i++) {
      if (ParticleType[i] == PARTICLE_TYPE_MUST_REFINE) {
        for (int dim = 0; dim < GridRank; dim++) {
          minr[dim] = min(minr[dim], ParticlePosition[dim][i]);
          maxr[dim] = max(maxr[dim], ParticlePosition[dim][i]);
        }
      }
    }
    printf("FlagByMass[G%"ISYM"]: MRP bounding box = (%lf %lf %lf) -> (%lf %lf %lf)\n", this->ID, minr[0], minr[1], minr[2], maxr[0], maxr[1], maxr[2]);
  }

  printf("FlagByMass[G%"ISYM"][Restrict=%"ISYM"]: nflag = %"ISYM"/%"ISYM", nunflag = %"ISYM", nonzero = %"ISYM"/%"ISYM"\n", this->ID, RestrictFlag, nflag, NumberOfFlaggedCells, nunflag, nonzero, size);
#endif

  /* remove MassFlaggingField. */
 
  if (ThisFlaggingMethod == 2) {
    delete [] MassFlaggingField;
    MassFlaggingField = NULL;
  } else if (ThisFlaggingMethod == 4) {

    delete [] ParticleMassFlaggingField;
    ParticleMassFlaggingField = NULL;
  }
 
  return NumberOfFlaggedCells;
 
}
