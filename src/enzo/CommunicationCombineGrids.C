/***********************************************************************
/
/  COMMUNICATION ROUTINE: COMBINE GRIDS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

#include <stdlib.h> 
 
/* function prototypes */
 
int CommunicationCombineGrids(HierarchyEntry *OldHierarchy,
			      HierarchyEntry **NewHierarchyPointer,
			      FLOAT WriteTime, int RestartDump = FALSE)
{
 
  /* If there is only one proc, then just point the new one at the old one. */
 
  if (NumberOfProcessors == 1 || ParallelRootGridIO == TRUE) {
    *NewHierarchyPointer = OldHierarchy;
    return SUCCESS;
  }
 
  /* Otherwise generate a new hierarchy entry and proceed. */
 
  HierarchyEntry *NewHierarchy = new HierarchyEntry;
  *NewHierarchyPointer = NewHierarchy;

  int NumberOfInterpolatedFields = 0;
  int Rank, dim, Dims[MAX_DIMENSION], NewDims[MAX_DIMENSION],
      SendOffset[MAX_DIMENSION], TempDims[MAX_DIMENSION],
      StartIndex[MAX_DIMENSION];
  FLOAT Left[MAX_DIMENSION], Right[MAX_DIMENSION], CellSize[MAX_DIMENSION];
 
  /* Compute dims, etc. for new grid assuming it fills entire domain. */
 
  OldHierarchy->GridData->ReturnGridInfo(&Rank, Dims, Left, Right);
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    if (dim < Rank)
      Dims[dim] -= 2*NumberOfGhostZones;
    CellSize[dim] = (Right[dim] - Left[dim])/FLOAT(Dims[dim]);
    NewDims[dim] = nint((DomainRightEdge[dim] - DomainLeftEdge[dim])/
			CellSize[dim])
                   + ((dim < Rank) ? 2*NumberOfGhostZones : 0);
  }
  if (debug)
    printf("CombineGrids: NewDims = %"ISYM" %"ISYM" %"ISYM"\n",
	   NewDims[0], NewDims[1], NewDims[2]);

  switch (OutputSmoothedDarkMatter) {
  case 1: NumberOfInterpolatedFields = 1; break;  // density
  case 2: NumberOfInterpolatedFields = 5; break;  // + rms velocity + 3-velocity
  }

  /* Generate a new grid. */
 
  NewHierarchy->GridData = new grid;
  NewHierarchy->GridData->InheritProperties(OldHierarchy->GridData);
  NewHierarchy->GridData->SetGravityParameters(
		       OldHierarchy->GridData->ReturnGravityBoundaryType());
  NewHierarchy->GridData->PrepareGrid(Rank, NewDims, DomainLeftEdge,
				      DomainRightEdge, 0);
 
  /* Loop over old grids and copy info. */
 
  HierarchyEntry *Temp = OldHierarchy;
  grid *NewGrid = NewHierarchy->GridData;

  if (MyProcessorNumber == ROOT_PROCESSOR) 
    NewGrid->WriteMCTP("ComCombine_NewGrid_A0"); //DEBUG

  int count = -1;//DEBUG

  while (Temp != NULL) {
      
      count++; // DEBUG
      printf("NEWGRIDA1 proc%d, count %i, NewGrid %p, Temp %p, TempNext %p\n", MyProcessorNumber, count, NewGrid, Temp, Temp->NextGridThisLevel);


    /* Compute region to send. */

    //if (MyProcessorNumber == ROOT_PROCESSOR){ // DEBUG
      char count_ch[1];
      char write_ch[10];
      //char rand_ch[5];
      //int rand_int = rand() % 100000;
      char *filename = new char[MAX_LINE_LENGTH];
      sprintf(count_ch, "%d", count);
      sprintf(write_ch, "%.0f", WriteTime);
      //sprintf(rand_ch, "%d", rand_int);
      strcpy(filename, "ComCombine_NewGrid_A1_c");
      strcat(filename, count_ch);
      strcat(filename, write_ch);
      //strcat(filename, rand_ch);
      NewGrid->WriteMCTP(filename); //DEBUG    
    //}

 
    grid *OldGrid = Temp->GridData;

    if (MyProcessorNumber == ROOT_PROCESSOR) 
      NewGrid->WriteMCTP("ComCombine_NewGrid_A2"); //DEBUG

    OldGrid->ReturnGridInfo(&Rank, TempDims, Left, Right);
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      SendOffset[dim] = (dim < Rank)? NumberOfGhostZones : 0;
      TempDims[dim] -= 2*SendOffset[dim];
      StartIndex[dim] = nint((Left[dim] - DomainLeftEdge[dim])/CellSize[dim])
	              + SendOffset[dim];
    }
 
    if (MyProcessorNumber == ROOT_PROCESSOR) 
      NewGrid->WriteMCTP("ComCombine_NewGrid_A3"); //DEBUG
    
    //Sometimes E isn't created by the time this code is run.Ensure it is.

    if(UseMHDCT == TRUE  && MyProcessorNumber == OldGrid->ReturnProcessorNumber() ){
      for(int field=0;field<3;field++)
        if( OldGrid->ElectricField[field] == NULL ){
          OldGrid->ElectricField[field] = new float[OldGrid->ElectricSize[field]];
          for(int i=0;i<OldGrid->ElectricSize[field];i++)
            OldGrid->ElectricField[field][i]=0.0;
        }
    }
    /* Copy grid region. */
 
    int RecvType = ((WriteTime < 0) && (RestartDump == FALSE)) ? 
                     NEW_ONLY : NEW_AND_OLD;
    int OldProc = OldGrid->ReturnProcessorNumber(),
        NewProc = NewGrid->ReturnProcessorNumber();
    CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
//  printf("(%"ISYM"): %"ISYM" --> %"ISYM"\n", MyProcessorNumber, OldProc, NewProc);
    if (MyProcessorNumber == NewProc || MyProcessorNumber == OldProc)
      if (NewGrid->CommunicationReceiveRegion(OldGrid, OldProc, ALL_FIELDS,
			      ((WriteTime < 0) ? NEW_ONLY : NEW_AND_OLD),
			      StartIndex, TempDims, FALSE) == FAIL) {
	ENZO_FAIL("Error in grid->CommunicationReceiveRegion.\n");
      }
    
    if (MyProcessorNumber == ROOT_PROCESSOR) 
      NewGrid->WriteMCTP("ComCombine_NewGrid_A4"); //DEBUG
  

    /* Copy interpolated (particle) regions */

    if ((MyProcessorNumber == NewProc || MyProcessorNumber == OldProc) &&
	NumberOfInterpolatedFields > 0) {
      OldGrid->ReturnGridInfo(&Rank, TempDims, Left, Right);
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	SendOffset[dim] = 0;
	TempDims[dim] -= 2*NumberOfGhostZones;
	StartIndex[dim] = nint((Left[dim] - DomainLeftEdge[dim])/CellSize[dim])
	  + SendOffset[dim];
      }
      if (NewGrid->CommunicationReceiveRegion(OldGrid, OldProc, INTERPOLATED_FIELDS,
			      NEW_ONLY, StartIndex, TempDims, FALSE) == FAIL) {
	ENZO_FAIL("Error in grid->CommunicationReceiveRegion.\n");
      }
    }

    if (MyProcessorNumber == ROOT_PROCESSOR) 
      NewGrid->WriteMCTP("ComCombine_NewGrid_A5"); //DEBUG    
 
    /* Copy particles. */
 
    if (OldGrid->CommunicationSendParticles(NewGrid, NewProc, 0,
	        OldGrid->ReturnNumberOfParticles(), -1) == FAIL) {
      ENZO_FAIL("Error in grid->CommunicationSendParticles.\n");
    }
    if (OldGrid->CommunicationSendActiveParticles(
            NewGrid, NewProc, false) == FAIL) {
      ENZO_FAIL("Error in grid->CommunicationSendActiveParticles.\n");
    }

    OldGrid->WriteMCTP("ComCombine_PreMoveCellZero"); //DEBUG

    if (OldGrid->MoveMonteCarloTracerParticlesToCellZero() == FAIL) {
      ENZO_FAIL("Error in grid->MoveMonteCarloTracerParticlesToCellZero.\n");
    }


    OldGrid->WriteMCTP("ComCombine_OldGrid_PreSend"); //DEBUG
    if (MyProcessorNumber == ROOT_PROCESSOR) 
      NewGrid->WriteMCTP("ComCombine_NewGrid_PreSend"); //DEBUG

    /* Send Monte Carlo Tracer Particles to new grid but keep particles in OldGrid */
    if (OldGrid->CommunicationSendMonteCarloTracerParticles(
            NewGrid, NewProc, 0) == FAIL) {
      ENZO_FAIL("Error in grid->CommunicationSendMonteCarloTracerParticles.\n");
    }

    OldGrid->WriteMCTP("ComCombine_OldGrid_PostSend"); //DEBUG
    if (MyProcessorNumber == ROOT_PROCESSOR) 
      NewGrid->WriteMCTP("ComCombine_NewGrid_PostSend"); //DEBUG

    if (OldGrid->DistributeMonteCarloTracerParticles() == FAIL) {
      ENZO_FAIL("Error in grid->DistributeMonteCarloTracerParticles.\n");
    } 

    // if (NewGrid->DistributeMonteCarloTracerParticles() == FAIL) {
    //   ENZO_FAIL("Error in grid->DistributeMonteCarloTracerParticles.\n");
    // }     
 
    OldGrid->WriteMCTP("ComCombine_OldGrid_PostDistribute"); //DEBUG
    if (MyProcessorNumber == ROOT_PROCESSOR) 
      NewGrid->WriteMCTP("ComCombine_NewGrid_PostDistribute"); //DEBUG
    /* Next Grid */
 
    Temp = Temp->NextGridThisLevel;
  }

//  printf("(%"ISYM"): done\n", MyProcessorNumber);
 
  /* Create a new first level of hierarchy entries that are all below the
     new one.  Below that, just point back into the old hierarchy. */
 
  NewHierarchy->ParentGrid = OldHierarchy->ParentGrid;
  NewHierarchy->NextGridNextLevel = NULL;
  NewHierarchy->NextGridThisLevel = NULL;
  Temp = OldHierarchy;
  HierarchyEntry *Previous = NULL;
  while (Temp != NULL) {
 
    HierarchyEntry *Temp2 = Temp->NextGridNextLevel;
 
    while (Temp2 != NULL) {
 
      HierarchyEntry *NewEntry = new HierarchyEntry;
      NewEntry->NextGridThisLevel = Previous;
      NewEntry->NextGridNextLevel = Temp2->NextGridNextLevel;
      NewEntry->ParentGrid        = NewHierarchy;
      NewEntry->GridData = Temp2->GridData;
      NewHierarchy->NextGridNextLevel = NewEntry;
      Previous = NewEntry;
 
      Temp2 = Temp2->NextGridThisLevel;
    }
 
    Temp = Temp->NextGridThisLevel;
  }
 
  return SUCCESS;
}
