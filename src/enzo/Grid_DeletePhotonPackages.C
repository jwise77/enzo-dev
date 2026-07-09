/***********************************************************************
/
/  GRID CLASS (DELETES ALL PHOTONS)
/
/  written by: John H. Wise
/  date:       November, 2005
/  modified1:
/
/  PURPOSE:  
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);

int grid::DeletePhotonPackages(int DeleteHeadPointer) {

  PhotonPackages.free_arrays();
  FinishedPhotonPackages.free_arrays();
  PausedPhotonPackages.free_arrays();

  this->NumberOfPhotonPackages = 0;

  return SUCCESS;
}
