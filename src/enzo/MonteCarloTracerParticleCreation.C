/***********************************************************************
/
/  CREATES MONTE CARL TRACER PARTICLES.
/
/  written by: Corey Brummel-Smith
/  date:       August, 2021
/
/  PURPOSE:
/
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
  
#include <string.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"
 
int MonteCarloTracerParticleCreation(FILE *fptr, HierarchyEntry &TopGrid,
			   TopGridData &MetaData)
{
  if (ParallelRootGridIO == FALSE) {
    // TODO: Extend to multiple grids
    if (TopGrid.GridData->CreateMonteCarloTracerParticles() == FAIL) {
      ENZO_FAIL("Error in grid->CreateMonteCarloTracerParticles.\n");
    }
  }
  return SUCCESS;
}
