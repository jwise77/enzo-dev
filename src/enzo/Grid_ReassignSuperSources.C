#include <stdlib.h>
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
#include "LevelHierarchy.h"

int FindSuperSource(PhotonPackageEntry **PP, int &LeafID, int SearchNewTree);
int FindSuperSourceByPosition(PhotonPackageEntry **PP);

int grid::ReassignSuperSources(void)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (PhotonPackages.numPackages == 0)
    return SUCCESS;

  bool outside = false;
  int dim, LeafID;
  float radius2, dx;
  FLOAT OldPosition[MAX_DIMENSION];

  for (int idx = 0; idx < PhotonPackages.numPackages; idx++) {

    if (PhotonPackages.CurrentSource[idx] == NULL)
      continue;

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      OldPosition[dim] = PhotonPackages.CurrentSource[idx]->Position[dim];

    // Reassign super source by leaf ID
    LeafID = PhotonPackages.CurrentSource[idx]->LeafID;

    // Create a temporary bridge package on the stack for FindSuperSource
    PhotonPackageEntry tempPP;
    tempPP.CurrentSource = PhotonPackages.CurrentSource[idx];
    tempPP.SourcePosition[0] = PhotonPackages.SourceX[idx];
    tempPP.SourcePosition[1] = PhotonPackages.SourceY[idx];
    tempPP.SourcePosition[2] = PhotonPackages.SourceZ[idx];
    tempPP.Radius = PhotonPackages.Radius[idx];

    PhotonPackageEntry *tempPPPtr = &tempPP;
    FindSuperSource(&tempPPPtr, LeafID, TRUE);

    radius2 = 0;
    if (tempPP.CurrentSource != NULL) {
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	dx = tempPP.CurrentSource->Position[dim] - OldPosition[dim];
	radius2 += dx*dx;
      }
      outside = (radius2 > tempPP.CurrentSource->ClusteringRadius *
		tempPP.CurrentSource->ClusteringRadius);
    }

    /* In the case where the leaf ID has changed (check by change in
       position), find super source by position.  If LeafID is
       undefined, then FindSuperSource couldn't locate the leaf by its
       ID, so we search by position. */

    if (outside || LeafID == INT_UNDEFINED || tempPP.CurrentSource == NULL)
      if (FindSuperSourceByPosition(&tempPPPtr) == FAIL) {
	ENZO_FAIL("Error in FindSuperSourceByPosition.\n");
      }

    // Write back the updated CurrentSource
    PhotonPackages.CurrentSource[idx] = tempPP.CurrentSource;

  } // ENDFOR photons

  return SUCCESS;

}
