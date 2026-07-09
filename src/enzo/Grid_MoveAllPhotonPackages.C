#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (MOVE ALL PHOTONPACKAGES FROM SPECIFIED GRIDS TO THIS GRID)
/
/  written by: Tom Abel
/  date:       May, 2005
/  modified1:
/
/  PURPOSE: 
/
/    NOTE: We assume all the from grids are at the same level!
/    NOTE: Modelled after GBs Grid_MoveAllParticles.C
************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);

int grid::MoveAllPhotonPackages(int NumberOfGrids, grid* FromGrid[])
{

  if (NumberOfGrids < 1) {
    ENZO_VFAIL("grid::MoveAllPhotonPackages: NumberOfGrids(%"ISYM") must be > 0.\n", 
	    NumberOfGrids)
  }

  /* Determine total number of particles. */

  int TotalNumberOfPackages = NumberOfPhotonPackages;
  int i, j, gridcount, dim, *Number, *Type;

  for (gridcount = 0; gridcount < NumberOfGrids; gridcount++) 
    TotalNumberOfPackages += FromGrid[gridcount]->NumberOfPhotonPackages;
  if (TotalNumberOfPackages == 0)
    return SUCCESS;

  /* Debugging info. */

//  if (debug)
//    fprintf(stdout, "MoveAllPackages: %"ISYM" (before: ThisGrid = %"ISYM").\n",
//	    TotalNumberOfPackages, NumberOfPhotonPackages);

  // go to end of List
  /* Error check number of photons.  If a bad value, reset photons */

  if (NumberOfPhotonPackages < 0) {
    printf("MoveAllPackages: WARNING. Resetting photons. "
	   "NumberOfPhotons = %"ISYM"\n", NumberOfPhotonPackages);
    NumberOfPhotonPackages = 0;
    TotalNumberOfPackages = 0;
    PhotonPackages.free_arrays();
    printf("MoveAllPackages: deleted photons\n");
    return SUCCESS;
  }

  /* Connect lists of FromGrids' PhotonPackages to this list */

  int count = NumberOfPhotonPackages;

  for (gridcount = 0; gridcount < NumberOfGrids; gridcount++) {

   /* If on the same processor, just add. */

    if (MyProcessorNumber == ProcessorNumber &&
        MyProcessorNumber == FromGrid[gridcount]->ProcessorNumber) {
      PhotonPackageSoA *fromSoA = FromGrid[gridcount]->ReturnPhotonPackagePointer();
      for (int i = 0; i < fromSoA->numPackages; i++) {
        PhotonPackages.append(fromSoA->Flux[i], fromSoA->Type[i], fromSoA->Energy[i],
                              fromSoA->CrossSection[i], fromSoA->TimeInterval[i],
                              fromSoA->EmissionTime[i], fromSoA->CurrentTime[i],
                              fromSoA->Radius[i], fromSoA->ColumnDensity[i],
                              fromSoA->PixelNum[i], fromSoA->Level[i],
                              fromSoA->SourceX[i], fromSoA->SourceY[i],
                              fromSoA->SourceZ[i], fromSoA->SourcePositionDiff[i],
                              fromSoA->CurrentSource[i]);
      }
      fromSoA->free_arrays();
    }
    /* Otherwise, communicate. */
    
    else {
      if (MyProcessorNumber == ProcessorNumber ||
          MyProcessorNumber == FromGrid[gridcount]->ProcessorNumber) {
	if (DEBUG)
	  printf("MoveAllPackages(%"ISYM"): (COMM) moving %"ISYM" PhotonPackages. "
		 "grid #%"ISYM" of %"ISYM".\n", 
		 MyProcessorNumber, FromGrid[gridcount]->NumberOfPhotonPackages,
		 gridcount, NumberOfGrids);      
	if (FromGrid[gridcount]->CommunicationSendPhotonPackages(this, 
	       ProcessorNumber, NumberOfPhotonPackages, 
               FromGrid[gridcount]->NumberOfPhotonPackages, &PhotonPackages) == FAIL) {
	  ENZO_FAIL("Error in grid->CommunicationSendPhotonPackages.\n");
	}
      }
    }

  } // end: loop over grids.

  /* Set new number of particles in this grid. */

  NumberOfPhotonPackages = TotalNumberOfPackages;

  for (gridcount = 0; gridcount < NumberOfGrids; gridcount++)
    FromGrid[gridcount]->SetNumberOfPhotonPackages(0);

  return SUCCESS;
}
