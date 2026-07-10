/***********************************************************************
/
/  COMMUNICATION ROUTINE: RECEIVE PHOTONS (NON-BLOCKING ASYNC)
/
/  written by: John H. Wise
/  date:       November, 2005
/  modified1:  July, 2026 (Refactored to modern MPI-3 async consensus)
/
/  PURPOSE: Polls and receives count and data messages from other ranks
/           non-blockingly, and unpacks the photon packages.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#define DEBUG 0
#include <stdlib.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "GroupPhotonList.h"
#include "PhotonCommunication.h"

#ifdef USE_MPI
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int FindSuperSource(PhotonPackageEntry **PP, int &LeafID, 
		    int SearchNewTree = TRUE);
#endif /* USE_MPI */

int CommunicationReceiverPhotons(LevelHierarchyEntry *LevelArray[], bool block)
{
#ifdef USE_MPI
  if (NumberOfProcessors == 1)
    return SUCCESS;

  int active_count = 0;
  MPI_Request requests[MAX_PH_RECEIVE_BUFFERS + NumberOfProcessors + 1];
  int map_type[MAX_PH_RECEIVE_BUFFERS + NumberOfProcessors + 1];
  int map_index[MAX_PH_RECEIVE_BUFFERS + NumberOfProcessors + 1];
  
  // 1. Add consensus request
  if (PH_ConsensusRequest != MPI_REQUEST_NULL) {
    map_type[active_count] = 0;
    map_index[active_count] = 0;
    requests[active_count++] = PH_ConsensusRequest;
  }
  
  // 2. Add count receives
  for (int proc = 0; proc < NumberOfProcessors; proc++) {
    if (proc != MyProcessorNumber && PhotonMessageRequest[proc] != MPI_REQUEST_NULL) {
      map_type[active_count] = 1;
      map_index[active_count] = proc;
      requests[active_count++] = PhotonMessageRequest[proc];
    }
  }
  
  // 3. Add data receives
  for (int i = 0; i < MAX_PH_RECEIVE_BUFFERS; i++) {
    if (PH_CommunicationReceiveMPI_Request[i] != MPI_REQUEST_NULL) {
      map_type[active_count] = 2;
      map_index[active_count] = i;
      requests[active_count++] = PH_CommunicationReceiveMPI_Request[i];
    }
  }

  if (active_count == 0) {
    if (DEBUG) {
      printf("P%d: CommunicationReceiverPhotons: active_count == 0. Returning.\n", MyProcessorNumber);
      fflush(stdout);
    }
    return SUCCESS;
  }

  if (DEBUG) {
    printf("P%d: CommunicationReceiverPhotons(block=%d): active_count=%d, PH_ConsensusRequest=%s\n",
           MyProcessorNumber, block, active_count, (PH_ConsensusRequest != MPI_REQUEST_NULL) ? "ACTIVE" : "NULL");
    fflush(stdout);
  }

  Eint32 out_count = 0;
  Eint32 indices[MAX_PH_RECEIVE_BUFFERS + NumberOfProcessors + 1];
  MPI_Status statuses[MAX_PH_RECEIVE_BUFFERS + NumberOfProcessors + 1];

  if (block) {
    MPI_Waitsome(active_count, requests, &out_count, indices, statuses);
  } else {
    MPI_Testsome(active_count, requests, &out_count, indices, statuses);
  }

  if (DEBUG) {
    printf("P%d: CommunicationReceiverPhotons: Waitsome/Testsome returned out_count=%d\n",
           MyProcessorNumber, out_count);
    fflush(stdout);
  }

  if (out_count <= 0)
    return SUCCESS;

  // Generate grid array for unpacking if we have data receives completed
  int *nGrids[MAX_DEPTH_OF_HIERARCHY];
  HierarchyEntry **Grids[MAX_DEPTH_OF_HIERARCHY];
  bool grid_array_generated = false;

  for (int i = 0; i < out_count; i++) {
    int idx = indices[i];
    int type = map_type[idx];
    int index = map_index[idx];

    if (DEBUG) {
      printf("P%d: CommunicationReceiverPhotons completed event index=%d, type=%d, index=%d\n",
             MyProcessorNumber, idx, type, index);
      fflush(stdout);
    }

    if (type == 0) {
      // Consensus request completed
      PH_ConsensusRequest = requests[idx]; // Should be MPI_REQUEST_NULL now
      if (DEBUG) {
        printf("P%d: Consensus barrier request completed!\n", MyProcessorNumber);
        fflush(stdout);
      }
    }
    else if (type == 1) {
      // Count receive completed
      PhotonMessageRequest[index] = requests[idx]; // Should be MPI_REQUEST_NULL now
      int num_messages = PhotonMessageBuffer[index];
      if (DEBUG) {
        printf("P%d: Count receive completed from P%d: num_messages = %d\n", MyProcessorNumber, index, num_messages);
        fflush(stdout);
      }
      
      for (int msg = 0; msg < num_messages; msg++) {
        int slot = -1;
        for (int s = 0; s < MAX_PH_RECEIVE_BUFFERS; s++) {
          if (PH_CommunicationReceiveMPI_Request[s] == MPI_REQUEST_NULL) {
            slot = s;
            break;
          }
        }
        if (slot == -1) ENZO_FAIL("Exceeded MAX_PH_RECEIVE_BUFFERS!");

        GroupPhotonList *ReceiveBuffer = new GroupPhotonList[PHOTON_BUFFER_SIZE];
        PH_CommunicationReceiveBuffer[slot] = (char *) ReceiveBuffer;
        PH_CommunicationReceiveIndex++;

        if (DEBUG) {
          printf("P%d: Posting data Irecv from P%d in slot %d\n", MyProcessorNumber, index, slot);
          fflush(stdout);
        }

        MPI_Irecv(ReceiveBuffer, PHOTON_BUFFER_SIZE, MPI_PhotonList, index,
                  MPI_PHOTONGROUP_TAG, MPI_COMM_WORLD, &PH_CommunicationReceiveMPI_Request[slot]);
      }
      
      // Re-post count receive
      MPI_Irecv(&PhotonMessageBuffer[index], 1, MPI_INT, index, MPI_NPHOTON_TAG, MPI_COMM_WORLD, &PhotonMessageRequest[index]);
    }
    else if (type == 2) {
      // Data receive completed
      PH_CommunicationReceiveMPI_Request[index] = requests[idx]; // Should be MPI_REQUEST_NULL now
      PH_WorkReceived = 1;

      if (!grid_array_generated) {
        for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
          if (LevelArray[level] != NULL) {
            Grids[level] = NULL;
            nGrids[level] = new int;
            *nGrids[level] = GenerateGridArray(LevelArray, level, &Grids[level]);
          } else {
            nGrids[level] = NULL;
            Grids[level] = NULL;
          }
        }
        grid_array_generated = true;
      }

      // Unpack
      GroupPhotonList *RecvBuffer = (GroupPhotonList *) PH_CommunicationReceiveBuffer[index];
      int num_receives = 0;
      while (RecvBuffer[num_receives].ToLevel != BUFFER_END && num_receives < PHOTON_BUFFER_SIZE) {
        num_receives++;
      }

      if (DEBUG) {
        printf("P%d: Unpacking %d packages from slot %d\n", MyProcessorNumber, num_receives, index);
        fflush(stdout);
      }

      for (int k = 0; k < num_receives; k++) {
        int lvl = RecvBuffer[k].ToLevel;
        int gi = RecvBuffer[k].ToGrid;

        if (gi >= *nGrids[lvl]) continue;
        grid *ToGrid = Grids[lvl][gi]->GridData;
        if (ToGrid->ReturnProcessorNumber() != MyProcessorNumber) continue;

        PhotonPackageSoA *ToPP = RecvBuffer[k].PausedPhoton ? ToGrid->ReturnPausedPackagePointer() : ToGrid->ReturnPhotonPackagePointer();
        ToPP->append(RecvBuffer[k].buffer.Photons, RecvBuffer[k].buffer.Type,
                     RecvBuffer[k].buffer.Energy, RecvBuffer[k].buffer.CrossSection,
                     RecvBuffer[k].buffer.EmissionTimeInterval, RecvBuffer[k].buffer.EmissionTime,
                     RecvBuffer[k].buffer.CurrentTime, RecvBuffer[k].buffer.Radius,
                     RecvBuffer[k].buffer.ColumnDensity, RecvBuffer[k].buffer.ipix,
                     RecvBuffer[k].buffer.level, RecvBuffer[k].buffer.SourcePosition[0],
                     RecvBuffer[k].buffer.SourcePosition[1], RecvBuffer[k].buffer.SourcePosition[2],
                     RecvBuffer[k].buffer.SourcePositionDiff, NULL);

        int idx = ToPP->numPackages - 1;
        if (RadiativeTransferSourceClustering) {
          PhotonPackageEntry tempPP;
          int leafID = RecvBuffer[k].buffer.SuperSourceID;
          PhotonPackageEntry *tempPPPtr = &tempPP;
          FindSuperSource(&tempPPPtr, leafID);
          ToPP->CurrentSource[idx] = tempPP.CurrentSource;
        } else {
          ToPP->CurrentSource[idx] = NULL;
        }

        ToGrid->SetNumberOfPhotonPackages(ToGrid->ReturnNumberOfPhotonPackages() + 1);
      }

      delete [] RecvBuffer;
      PH_CommunicationReceiveBuffer[index] = NULL;
      PH_CommunicationReceiveIndex--;
    }
  }

  if (grid_array_generated) {
    for (int level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
      if (LevelArray[level] != NULL) {
        delete [] Grids[level];
        delete nGrids[level];
      }
    }
  }
#endif /* USE_MPI */
  return SUCCESS;
}
