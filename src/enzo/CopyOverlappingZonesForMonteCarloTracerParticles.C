/***********************************************************************
/
/  COPY OVERLAPPING GRIDS FUNCTION (Clone of CopyOverlappingZones but
/                                   for Monte Carlo Tracer Particles)
/
/  written by: Corey Brummel-Smith
/  date:       September, 2021
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
#include "TopGridData.h"
#include "LevelHierarchy.h"
 
 
int CopyOverlappingZonesForMonteCarloTracerParticles(grid* CurrentGrid, 
  TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], int level)
{
 
  LevelHierarchyEntry *Temp = LevelArray[level];
 
  while (Temp != NULL) {
    if (CurrentGrid->CheckForOverlap(Temp->GridData,
                     MetaData->LeftFaceBoundaryCondition,
                     MetaData->RightFaceBoundaryCondition,
                     &grid::CopyMonteCarloTracerParticlesFromGrid)
    == FAIL) {
      fprintf(stderr, "Error in grid->CopyMonteCarloTracerParticlesFromGrid.\n");
    }
    Temp = Temp->NextGridThisLevel;
  }
 
  return SUCCESS;
}
