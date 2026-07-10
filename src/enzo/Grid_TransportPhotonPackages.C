#define DEBUG 0
#define MYPROC MyProcessorNumber == ProcessorNumber
/***********************************************************************
/
/  GRID CLASS (TRANSPORT PHOTON PACKAGES)
/
/  written by: Tom Abel
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: This is the heart of the radiative transfer algorithm.
/    On each Grid we initialize photo and heating rates and then call
/    WalkPhotonPackage so all photon packages are transported along their
/    own directions and the photo-ionization and heating rates on 
/    on the grid are updated on the fly. 
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "phys_constants.h"

void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);
PhotonPackageEntry *PopPhoton(PhotonPackageEntry * &Node);
PhotonPackageEntry *DeletePhotonPackage(PhotonPackageEntry *PP);
int FindField(int field, int farray[], int numfields);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::TransportPhotonPackages(int level, int finest_level, 
				  ListOfPhotonsToMove **PhotonsToMove, 
				  int GridNum, grid **Grids0, int nGrids0, 
				  grid *ParentGrid, grid *CurrentGrid)
{

  int i,j,k, dim, index, count;
  grid *MoveToGrid;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0 || MultiSpecies < 1) 
    return SUCCESS;

  if (RadiativeTransfer < 1) 
    return SUCCESS;

  if (RadiativeTransfer > 0 && GridRank < 3) {
    ENZO_FAIL("Transfer in less than 3D is not implemented!\n");
  }

  if (PhotonPackages.numPackages == 0)
    return SUCCESS;

  /* Get units. */
  double MassUnits, RT_Units;
  float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits, 
    DensityUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, PhotonTime) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
  MassUnits = (double) DensityUnits * POW(LengthUnits, 3.0);
  RT_Units = (double) TimeUnits * POW(LengthUnits, -3.0);

  /* speed of light in code units. note this one is independent of
     a(t), and Modify the photon propagation speed by this
     parameter */

  float LightSpeed;
  LightSpeed = RadiativeTransferPropagationSpeedFraction * 
    (clight/VelocityUnits);

  float DomainWidth[MAX_DIMENSION];
  //double FinestCellVolume = pow(RefineBy, -3*(finest_level-level));
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    //FinestCellVolume *= CellWidth[dim][0];
  }

  // if (DEBUG) fprintf(stdout,"TransportPhotonPackage: initialize fields.\n");
  // if (DEBUG) fprintf(stdout,"TransportPhotonPackage: %"ISYM" %"ISYM" .\n",
  // 		     GridStartIndex[0], GridEndIndex[0]);

  if (DEBUG) {
    fprintf(stdout, "TransportPhotonPackage: done initializing.\n");
    fprintf(stdout, "[%d] counted %"ISYM" packages\n", this->ID, PhotonPackages.numPackages);
  }

  /* If requested, make vertex centered field (only when it doesn't
     exist ... see inside routine). */

  if (RadiativeTransferInterpolateField)
    for (i = 0; i < NumberOfBaryonFields; i++)
      if (FieldsToInterpolate[i] == TRUE)
	if (this->ComputeVertexCenteredField(i) == FAIL) {
	  ENZO_VFAIL("Error in grid->ComputeVertexCenteredField "
		  "(field %"ISYM").\n", i)
	}

  /* Calculate minimum photon flux before a ray is deleted */
  
  const double alphaB = 2.6e-13;
  float MinimumPhotonFlux, RecombinationTime;
  int gmethod = INT_UNDEFINED;
  for (i = 0; i < MAX_FLAGGING_METHODS; i++)
    if (CellFlaggingMethod[i] == 2) gmethod = i;
  if (gmethod == INT_UNDEFINED)
    MinimumPhotonFlux = POW(TopGridDx[0], GridRank) * POW(RefineBy, level); // estimate
  else
    MinimumPhotonFlux = MinimumMassForRefinement[gmethod] * 
      POW(RefineBy, level*MinimumMassForRefinementLevelExponent[gmethod]);
  if (ComovingCoordinates)
    MinimumPhotonFlux *= (float) ((RT_Units / TimeUnits) * (MassUnits / mh) * dtPhoton / 
				  (PhotonTime * RadiativeTransferHubbleTimeFraction));
  else {
    RecombinationTime = 1.0 / (alphaB * DensityUnits / mh) / TimeUnits;
    MinimumPhotonFlux *= (float) ((RT_Units / TimeUnits) * (MassUnits / mh) * dtPhoton / 
				  (10*RecombinationTime));
  }
  
  int dcount = 0;
  int tcount = 0;
  int pcount = 0;
  int trcount = 0;
  int DeleteMe, DeltaLevel, PauseMe;
  int prev_type = -1;
  float LightCrossingTime = RadiativeTransferRayMaximumLength * (VelocityUnits) /
    (clight * RadiativeTransferPropagationSpeedFraction); 
  FLOAT EndTime;
  if (MYPROC && DEBUG) {
    printf("RadiativeTransferRayMaximumLength = %g\t  RadiativeTransferPropagationSpeedFraction= %g\n",  RadiativeTransferRayMaximumLength, RadiativeTransferPropagationSpeedFraction);
    printf("LightCrossingTime = %f\n", LightCrossingTime);
  }
  if (RadiativeTransferAdaptiveTimestep)
    EndTime = PhotonTime+LightCrossingTime;
  else
    EndTime = PhotonTime+dtPhoton-PFLOAT_EPSILON;

  for (int idx = 0; idx < PhotonPackages.numPackages; idx++) {
    // Create temporary PhotonPackageEntry on the stack to bridge with the existing WalkPhotonPackage
    PhotonPackageEntry tempPP;
    tempPP.Photons = PhotonPackages.Flux[idx];
    tempPP.Type = PhotonPackages.Type[idx];
    tempPP.Energy = PhotonPackages.Energy[idx];
    tempPP.CrossSection = PhotonPackages.CrossSection[idx];
    tempPP.EmissionTimeInterval = PhotonPackages.TimeInterval[idx];
    tempPP.EmissionTime = PhotonPackages.EmissionTime[idx];
    tempPP.CurrentTime = PhotonPackages.CurrentTime[idx];
    tempPP.Radius = PhotonPackages.Radius[idx];
    tempPP.ColumnDensity = PhotonPackages.ColumnDensity[idx];
    tempPP.ipix = PhotonPackages.PixelNum[idx];
    tempPP.level = PhotonPackages.Level[idx];
    tempPP.SourcePosition[0] = PhotonPackages.SourceX[idx];
    tempPP.SourcePosition[1] = PhotonPackages.SourceY[idx];
    tempPP.SourcePosition[2] = PhotonPackages.SourceZ[idx];
    tempPP.SourcePositionDiff = PhotonPackages.SourcePositionDiff[idx];
    tempPP.CurrentSource = PhotonPackages.CurrentSource[idx];
    
    // Set dummy PreviousPackage/NextPackage to satisfy internal linked list validation in WalkPhotonPackage
    PhotonPackageEntry dummyPrev;
    PhotonPackageEntry dummyNext;
    tempPP.PreviousPackage = &dummyPrev;
    tempPP.NextPackage = &dummyNext;
    dummyPrev.NextPackage = &tempPP;
    dummyNext.PreviousPackage = &tempPP;

    PhotonPackageEntry *tempPPPtr = &tempPP;
    int retval = 0;
    DeleteMe = FALSE;
    PauseMe = FALSE;
    MoveToGrid = NULL;

    if (MYPROC && DEBUG) {
      if (prev_type != tempPP.Type) {
	fprintf(stdout, "%s: Radiation type = %d\n", __FUNCTION__, tempPP.Type);
	prev_type = tempPP.Type;
      }
    }

    if (tempPP.CurrentTime < EndTime) {
      retval = WalkPhotonPackage(&tempPPPtr,
				 &MoveToGrid, ParentGrid, CurrentGrid, Grids0, nGrids0,
				 DeleteMe, PauseMe, DeltaLevel, LightCrossingTime,
				 LightSpeed, level, MinimumPhotonFlux);
      tcount++;
    } else {
      /* If all work is finished, store in FinishedPhotonPackages and remove from active */
      FinishedPhotonPackages.append(tempPP);
      DeleteMe = TRUE;
    }

    if (PauseMe == TRUE) {
      if (DEBUG > 1) fprintf(stdout, "paused photon\n");
      this->RegridPausedPhotonPackage(&tempPPPtr, ParentGrid, &MoveToGrid, DeltaLevel,
				      DeleteMe, DomainWidth, LightSpeed);

      // Insert in paused photon list if it belongs in this grid.
      if (MoveToGrid == NULL && DeleteMe == FALSE) {
	PausedPhotonPackages.append(tempPP);
	DeleteMe = TRUE;
      }
      pcount++;
    }

    if (MoveToGrid != NULL) {
      if (DEBUG > 1) {
	fprintf(stdout, "moving photon from %p to %p\n", CurrentGrid, MoveToGrid);
      }
      ListOfPhotonsToMove *NewEntry = new ListOfPhotonsToMove;
      NewEntry->NextPackageToMove = (*PhotonsToMove)->NextPackageToMove;
      (*PhotonsToMove)->NextPackageToMove = NewEntry;
      
      // We must copy the ray to a standalone package to put in the move list
      PhotonPackageEntry *movedPP = new PhotonPackageEntry(tempPP);
      movedPP->PreviousPackage = NULL;
      movedPP->NextPackage = NULL;
      
      NewEntry->PhotonPackage = movedPP;
      NewEntry->FromGrid = CurrentGrid;
      NewEntry->ToGrid   = MoveToGrid;
      NewEntry->ToGridNum= MoveToGrid->GetGridID();
      NewEntry->ToLevel  = level + DeltaLevel;
      NewEntry->ToProcessor = MoveToGrid->ReturnProcessorNumber();
      NewEntry->PausedPhoton = PauseMe ? TRUE : FALSE;
      
      if (NewEntry->ToProcessor >= NumberOfProcessors ||
	  NewEntry->ToProcessor < 0) {
	tempPP.PrintInfo();
	ENZO_VFAIL("Grid %d, Invalid ToProcessor P%d", GridNum, 
		   NewEntry->ToProcessor)
      }
      trcount++;
      DeleteMe = TRUE;
    } // ENDIF MoveToGrid

    if (DeleteMe == TRUE) {
      if (DEBUG > 1) fprintf(stdout, "delete photon\n");
      dcount++;
      PhotonPackages.DeletePackage(idx);
      idx--; // Decrement to reprocess this index now occupied by the swapped element
    } else {
      // Write back modified fields
      PhotonPackages.Flux[idx] = tempPP.Photons;
      PhotonPackages.Type[idx] = tempPP.Type;
      PhotonPackages.Energy[idx] = tempPP.Energy;
      PhotonPackages.CrossSection[idx] = tempPP.CrossSection;
      PhotonPackages.TimeInterval[idx] = tempPP.EmissionTimeInterval;
      PhotonPackages.EmissionTime[idx] = tempPP.EmissionTime;
      PhotonPackages.CurrentTime[idx] = tempPP.CurrentTime;
      PhotonPackages.Radius[idx] = tempPP.Radius;
      PhotonPackages.ColumnDensity[idx] = tempPP.ColumnDensity;
      PhotonPackages.PixelNum[idx] = tempPP.ipix;
      PhotonPackages.Level[idx] = tempPP.level;
      PhotonPackages.SourceX[idx] = tempPP.SourcePosition[0];
      PhotonPackages.SourceY[idx] = tempPP.SourcePosition[1];
      PhotonPackages.SourceZ[idx] = tempPP.SourcePosition[2];
      PhotonPackages.SourcePositionDiff[idx] = tempPP.SourcePositionDiff;
      PhotonPackages.CurrentSource[idx] = tempPP.CurrentSource;
    }

    // Retrieve any child rays that were created by splitting
    if (tempPP.NextPackage != &dummyNext) {
      PhotonPackageEntry *currChild = tempPP.NextPackage;
      while (currChild != &dummyNext) {
        PhotonPackages.append(*currChild);
        PhotonPackageEntry *nextChild = currChild->NextPackage;
        delete currChild;
        currChild = nextChild;
      }
    }
  } // ENDFOR active packages

  if (DEBUG)
    fprintf(stdout, "grid::TransportPhotonPackage[%d]: "
	    "transported %"ISYM" deleted %"ISYM" paused %"ISYM" moved %"ISYM"\n",
	    this->ID, tcount, dcount, pcount, trcount);
  
  NumberOfPhotonPackages = PhotonPackages.numPackages + PausedPhotonPackages.numPackages + FinishedPhotonPackages.numPackages;

  return SUCCESS;
}

