/***********************************************************************
/
/  ACTIVE PARTICLE INITIALIZATION
/
/  written by: Corey Brummel-Smith
/  date:       September, 2021
/
/  PURPOSE: Reset ExchangedThisTimestep for all Monte Carlo tracer
/           particles on all grids on this level.
/
************************************************************************/
#ifdef USE_MPI
#endif
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "performance.h"
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
#include "CommunicationUtilities.h"

/* prototypes */

int MonteCarloTracerParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
               int NumberOfGrids, LevelHierarchyEntry *LevelArray[], int level)
{
  /* Reset ExchangedThisTimestep for all MC tracer particles on all grids */
  for (int grid1 = 0; grid1 < NumberOfGrids; grid1++) {
    Grids[grid1]->GridData->ResetMonteCarloTracerParticlesExchangedThisTimestep();
  }
  return SUCCESS;
}
