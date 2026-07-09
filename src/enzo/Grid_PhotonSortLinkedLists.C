/***********************************************************************
/
/  GRID CLASS (SORT LINKED LISTS OF PHOTONS)
/
/  written by: Stephen Skory
/  date:       June, 2011
/  modified1:
/
/  PURPOSE:  Sorts the linked lists of arrays.
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

Eint32 compare_ss(const void *a, const void *b);
PhotonPackageEntry *LinkedListToArray(PhotonPackageEntry *Node, int n);
PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);
void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);

int grid::PhotonSortLinkedLists(void)
{

  if (MyProcessorNumber != ProcessorNumber) return SUCCESS;
  
  int nphotons, dim, count, i;
  PhotonPackageEntry *TempPP;
  
  // PhotonPackages
  nphotons = PhotonPackages.numPackages;
  if (nphotons > 0) {
    TempPP = new PhotonPackageEntry[nphotons];
    for (i = 0; i < nphotons; i++) {
      TempPP[i].Photons = PhotonPackages.Flux[i];
      TempPP[i].Type = PhotonPackages.Type[i];
      TempPP[i].Energy = PhotonPackages.Energy[i];
      TempPP[i].CrossSection = PhotonPackages.CrossSection[i];
      TempPP[i].EmissionTimeInterval = PhotonPackages.TimeInterval[i];
      TempPP[i].EmissionTime = PhotonPackages.EmissionTime[i];
      TempPP[i].CurrentTime = PhotonPackages.CurrentTime[i];
      TempPP[i].Radius = PhotonPackages.Radius[i];
      TempPP[i].ColumnDensity = PhotonPackages.ColumnDensity[i];
      TempPP[i].ipix = PhotonPackages.PixelNum[i];
      TempPP[i].level = PhotonPackages.Level[i];
      TempPP[i].SourcePositionDiff = PhotonPackages.SourcePositionDiff[i];
      TempPP[i].SourcePosition[0] = PhotonPackages.SourceX[i];
      TempPP[i].SourcePosition[1] = PhotonPackages.SourceY[i];
      TempPP[i].SourcePosition[2] = PhotonPackages.SourceZ[i];
      TempPP[i].CurrentSource = PhotonPackages.CurrentSource[i];
    }
    qsort(TempPP, nphotons, sizeof(PhotonPackageEntry), compare_ss);
    PhotonPackages.free_arrays();
    for (count = 0; count < nphotons; count++) {
      PhotonPackages.append(TempPP[count].Photons, TempPP[count].Type, TempPP[count].Energy,
                            TempPP[count].CrossSection, TempPP[count].EmissionTimeInterval,
                            TempPP[count].EmissionTime, TempPP[count].CurrentTime,
                            TempPP[count].Radius, TempPP[count].ColumnDensity,
                            TempPP[count].ipix, TempPP[count].level,
                            TempPP[count].SourcePosition[0], TempPP[count].SourcePosition[1],
                            TempPP[count].SourcePosition[2], TempPP[count].SourcePositionDiff,
                            TempPP[count].CurrentSource);
    }
    delete [] TempPP;
  }
  
  // FinishedPhotonPackages
  nphotons = FinishedPhotonPackages.numPackages;
  if (nphotons > 0) {
    TempPP = new PhotonPackageEntry[nphotons];
    for (i = 0; i < nphotons; i++) {
      TempPP[i].Photons = FinishedPhotonPackages.Flux[i];
      TempPP[i].Type = FinishedPhotonPackages.Type[i];
      TempPP[i].Energy = FinishedPhotonPackages.Energy[i];
      TempPP[i].CrossSection = FinishedPhotonPackages.CrossSection[i];
      TempPP[i].EmissionTimeInterval = FinishedPhotonPackages.TimeInterval[i];
      TempPP[i].EmissionTime = FinishedPhotonPackages.EmissionTime[i];
      TempPP[i].CurrentTime = FinishedPhotonPackages.CurrentTime[i];
      TempPP[i].Radius = FinishedPhotonPackages.Radius[i];
      TempPP[i].ColumnDensity = FinishedPhotonPackages.ColumnDensity[i];
      TempPP[i].ipix = FinishedPhotonPackages.PixelNum[i];
      TempPP[i].level = FinishedPhotonPackages.Level[i];
      TempPP[i].SourcePositionDiff = FinishedPhotonPackages.SourcePositionDiff[i];
      TempPP[i].SourcePosition[0] = FinishedPhotonPackages.SourceX[i];
      TempPP[i].SourcePosition[1] = FinishedPhotonPackages.SourceY[i];
      TempPP[i].SourcePosition[2] = FinishedPhotonPackages.SourceZ[i];
      TempPP[i].CurrentSource = FinishedPhotonPackages.CurrentSource[i];
    }
    qsort(TempPP, nphotons, sizeof(PhotonPackageEntry), compare_ss);
    FinishedPhotonPackages.free_arrays();
    for (count = 0; count < nphotons; count++) {
      FinishedPhotonPackages.append(TempPP[count].Photons, TempPP[count].Type, TempPP[count].Energy,
                                    TempPP[count].CrossSection, TempPP[count].EmissionTimeInterval,
                                    TempPP[count].EmissionTime, TempPP[count].CurrentTime,
                                    TempPP[count].Radius, TempPP[count].ColumnDensity,
                                    TempPP[count].ipix, TempPP[count].level,
                                    TempPP[count].SourcePosition[0], TempPP[count].SourcePosition[1],
                                    TempPP[count].SourcePosition[2], TempPP[count].SourcePositionDiff,
                                    TempPP[count].CurrentSource);
    }
    delete [] TempPP;
  }

  // PausedPhotonPackages
  nphotons = PausedPhotonPackages.numPackages;
  if (nphotons > 0) {
    TempPP = new PhotonPackageEntry[nphotons];
    for (i = 0; i < nphotons; i++) {
      TempPP[i].Photons = PausedPhotonPackages.Flux[i];
      TempPP[i].Type = PausedPhotonPackages.Type[i];
      TempPP[i].Energy = PausedPhotonPackages.Energy[i];
      TempPP[i].CrossSection = PausedPhotonPackages.CrossSection[i];
      TempPP[i].EmissionTimeInterval = PausedPhotonPackages.TimeInterval[i];
      TempPP[i].EmissionTime = PausedPhotonPackages.EmissionTime[i];
      TempPP[i].CurrentTime = PausedPhotonPackages.CurrentTime[i];
      TempPP[i].Radius = PausedPhotonPackages.Radius[i];
      TempPP[i].ColumnDensity = PausedPhotonPackages.ColumnDensity[i];
      TempPP[i].ipix = PausedPhotonPackages.PixelNum[i];
      TempPP[i].level = PausedPhotonPackages.Level[i];
      TempPP[i].SourcePositionDiff = PausedPhotonPackages.SourcePositionDiff[i];
      TempPP[i].SourcePosition[0] = PausedPhotonPackages.SourceX[i];
      TempPP[i].SourcePosition[1] = PausedPhotonPackages.SourceY[i];
      TempPP[i].SourcePosition[2] = PausedPhotonPackages.SourceZ[i];
      TempPP[i].CurrentSource = PausedPhotonPackages.CurrentSource[i];
    }
    qsort(TempPP, nphotons, sizeof(PhotonPackageEntry), compare_ss);
    PausedPhotonPackages.free_arrays();
    for (count = 0; count < nphotons; count++) {
      PausedPhotonPackages.append(TempPP[count].Photons, TempPP[count].Type, TempPP[count].Energy,
                                  TempPP[count].CrossSection, TempPP[count].EmissionTimeInterval,
                                  TempPP[count].EmissionTime, TempPP[count].CurrentTime,
                                  TempPP[count].Radius, TempPP[count].ColumnDensity,
                                  TempPP[count].ipix, TempPP[count].level,
                                  TempPP[count].SourcePosition[0], TempPP[count].SourcePosition[1],
                                  TempPP[count].SourcePosition[2], TempPP[count].SourcePositionDiff,
                                  TempPP[count].CurrentSource);
    }
    delete [] TempPP;
  }

  return SUCCESS;
}
