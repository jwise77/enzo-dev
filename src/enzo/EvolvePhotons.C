/***********************************************************************
/
/  EVOLVE PHOTONS FUNCTION
/
/  written by: Tom Abel
/  date:       May 2004
/  modified1:  November 2005 by John Wise (parallelized it)
/                
/
/  PURPOSE:
/    This routine is the main photon evolution function. 
/    From here we call first the routines that control the emission:
/    grid::Shine
/  while (stil_photons_to_update)
/    transport all photon packages local to a grid
/    communicate all surviving photon packages to their new parent grids. 
/  endwhile  
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "preincludes.h"
#include "performance.h"
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
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
#include "CommunicationUtilities.h"

/* function prototypes */
void my_exit(int status);
int RadiationFieldCalculateRates(FLOAT Time);
int CommunicationReceiverPhotons(LevelHierarchyEntry *LevelArray[],
				 bool block = false);
int CommunicationTransferPhotons(LevelHierarchyEntry *LevelArray[], 
				 ListOfPhotonsToMove **AllPhotons, 
				 char *kt_global,
				 int &keep_transporting);
int RadiativeTransferLoadBalanceRevert(HierarchyEntry **Grids[], int *NumberOfGrids);
int CommunicationLoadBalancePhotonGrids(HierarchyEntry **Grids[], int *NumberOfGrids,
					int FirstTimeAfterRestart);
int RadiativeTransferMoveLocalPhotons(ListOfPhotonsToMove **AllPhotons,
				      int &keep_transporting);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int InitializePhotonCommunication(void);
int FinalizePhotonCommunication(void);
int CommunicationBufferedSendActiveCount(void);
int CommunicationBufferPurge(void);
RadiationSourceEntry* DeleteRadiationSource(RadiationSourceEntry *RS);
PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);
int CreateSourceClusteringTree(int nShine, SuperSourceData *SourceList,
			       LevelHierarchyEntry *LevelArray[]);
int CommunicationSyncNumberOfPhotons(LevelHierarchyEntry *LevelArray[]);
int RadiativeTransferComputeTimestep(LevelHierarchyEntry *LevelArray[],
				     TopGridData *MetaData, float dtLevelAbove,
				     int level);
int StarParticleFindAll(LevelHierarchyEntry *LevelArray[], Star *&AllStars);
int SetSubgridMarker(TopGridData &MetaData, 
		     LevelHierarchyEntry *LevelArray[], int level,
		     int UpdateReplicatedGridsOnly);
void PrintMemoryUsage(char *str);
void fpcol(Eflt64 *x, int n, int m, FILE *log_fptr);
double ReturnWallTime();



//#define NONBLOCKING_RT_OFF  // moved to a compile-time define
#define REPORT_PERF
#define MAX_ITERATIONS 1000000
#define DEBUG_RT 0

#ifdef REPORT_PERF
#define START_PERF() tt0 = ReturnWallTime();
#else
#define START_PERF() ;
#endif
#ifdef REPORT_PERF
#define END_PERF(A) \
  tt1 = ReturnWallTime(); \
  PerfCounter[A] += tt1-tt0;
#else
#define END_PERF(A) ;
#endif

/* EvolvePhotons function */
int EvolvePhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		  Star *&AllStars, FLOAT GridTime, int level, int LoopTime)
{

  struct ThinGridList
  {
    grid *ThisGrid;
    ThinGridList *NextGrid;
  };

  bool FirstTime = true;

#ifdef REPORT_PERF
  double ep0, tt0, tt1, PerfCounter[14];
  ep0 = ReturnWallTime();
  for (int i = 0; i < 14; i++)
    PerfCounter[i] = 0;
#endif

  if (!RadiativeTransfer)
    return SUCCESS;

  /* Only call on the finest level */

  if (LevelArray[level+1] != NULL && LoopTime)
    return SUCCESS;


  if (dtPhoton < 0)
    return SUCCESS;  

  if (GridTime <= PhotonTime)
    return SUCCESS;

  LCAPERF_START("EvolvePhotons");
  TIMER_START("EvolvePhotons");

  /* Declarations */



  /* For early termination with a background, calculate background
     intensities */

  RadiationFieldCalculateRates(PhotonTime+0.5*dtPhoton);

  TIMER_REGISTER("RayTracing");
  TIMER_REGISTER("RayCommunication");

  grid *Helper;
  int i, lvl, GridNum;
  LevelHierarchyEntry *Temp;
  RadiationSourceEntry *RS;

  /* Find the finest grid that hosts each radiation source, and store
     them in a GridList */

  ThinGridList *RS_GridList, *NewNode, *TempGridList, *TailNode, *LastNode, 
    *Destroyer;
  HierarchyEntry **Grids[MAX_DEPTH_OF_HIERARCHY];
  int nGrids[MAX_DEPTH_OF_HIERARCHY];
  for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    if (LevelArray[lvl] != NULL)
      nGrids[lvl] = GenerateGridArray(LevelArray, lvl, &Grids[lvl]);
    else
      nGrids[lvl] = 0;

  RS_GridList = new ThinGridList;
  RS_GridList->NextGrid = NULL;
  TailNode = RS_GridList;
  for (RS = GlobalRadiationSources->NextSource; RS; RS = RS->NextSource) {
    // Search for grid, if not defined. Skip sources that aren't born
    // for PhotonTest
    if (RS->GridID == INT_UNDEFINED && 
	!(RS->CreationTime > PhotonTime && ProblemType == 50)) {
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0; lvl--) {
	for (GridNum = 0; GridNum < nGrids[lvl]; GridNum++) {
	  if (MyProcessorNumber == Grids[lvl][GridNum]->GridData->
	      ReturnProcessorNumber())
	    if (Grids[lvl][GridNum]->GridData->PointInGrid(RS->Position)) {
	      RS->GridID = GridNum;
	      RS->GridLevel = lvl;
	    } // ENDIF PointInGrid
	} // ENDFOR grids
      } // ENDFOR level
    } // ENDIF undefined grid

    /* Now we know which grid, we can add to the GridList */

    NewNode = new ThinGridList;
    if (RS->GridLevel != INT_UNDEFINED)
      NewNode->ThisGrid = Grids[RS->GridLevel][RS->GridID]->GridData;
    else
      NewNode->ThisGrid = NULL;
    NewNode->NextGrid = NULL;
    TailNode->NextGrid = NewNode;
    TailNode = NewNode;

  } // ENDFOR RS

  /**********************************************************************
                       MAIN RADIATION TRANSPORT LOOP
   **********************************************************************/

  int NumberOfSources = 0;
  bool DeleteSources = FALSE;

  while (GridTime > PhotonTime) {

    /* Recalculate timestep if this isn't the first loop.  We already
       did this in RadiativeTransferPrepare */

    FLOAT dtLevelAbove;
    if (!FirstTime) {
      dtLevelAbove = LevelArray[level]->GridData->ReturnTimeStep();
      RadiativeTransferComputeTimestep(LevelArray, MetaData, dtLevelAbove, level);
    }

    if (debug && LoopTime == TRUE)
      printf("EvolvePhotons[%"ISYM"]: dt = %"GSYM", Time = %"FSYM", ", 
	     level, dtPhoton, PhotonTime);
      
    /* delete source if we are passed (or before) their lifetime (only
       if not restarting).  We must remove the host grid from the grid
       list as well. */

    RS = GlobalRadiationSources->NextSource;
    LastNode = RS_GridList;
    TempGridList = RS_GridList->NextGrid;
    NumberOfSources = 0;
    DeleteSources = (FirstTime ||
			  !(RadiativeTransferAdaptiveTimestep == FALSE &&
			    RadiativeTransferSourceClustering == TRUE));
			  
    while (RS != NULL) {
      if ( ((RS->CreationTime + RS->LifeTime) < PhotonTime) && LoopTime == TRUE &&
	   DeleteSources) {
	if (debug) {
	  fprintf(stdout, "\nEvolvePhotons: Deleted Source on lifetime limit \n");
	  fprintf(stdout, "EvolvePhotons:  %"GSYM" %"GSYM" %"GSYM" \n",
		  RS->CreationTime, RS->LifeTime, PhotonTime);
	}
	RS = DeleteRadiationSource(RS);
	
	// Remove the grid from the list
	LastNode->NextGrid = TempGridList->NextGrid;
	Destroyer = TempGridList;
	TempGridList = TempGridList->NextGrid;
	delete Destroyer;

      } else {
	if (RS->CreationTime < PhotonTime)
	  NumberOfSources++;                 // count sources
	RS = RS->NextSource;
	LastNode = TempGridList;
	TempGridList = TempGridList->NextGrid;
      }
    }

    if (debug) fprintf(stdout, "%"ISYM" SRC(s)\n", NumberOfSources);

    if (NumberOfSources == 0) {
      PhotonTime += dtPhoton;
      continue;
    }

  /* Temporarily load balance grids according to the number of ray
     segments.  We'll move the grids back at the end of this
     routine */

    if (RadiativeTransferLoadBalance && NumberOfSources > 0) {
      CommunicationLoadBalancePhotonGrids(Grids, nGrids, 
					  MetaData->FirstTimestepAfterRestart);
    }

    /* Initialize radiation fields */

    START_PERF();
    for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) 
	if (Temp->GridData->InitializeRadiativeTransferFields() == FAIL) {
	  ENZO_FAIL("Error in InitializeRadiativeTransferFields.\n");
	}
    END_PERF(0);

    /* create temperature fields for Compton heating */  

    if (RadiationXRayComptonHeating)  
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) 
	  if (Temp->GridData->InitializeTemperatureFieldForComptonHeating() == FAIL) {  
	    ENZO_FAIL("Error in InitializeTemperatureFieldForComptonHeating.\n");
	  }	
    /* Initialize Temperature Field for H2 shielding approximation */
    if(RadiativeTransferH2ShieldType == 1 || ProblemType == 50) {
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) 
	  if (Temp->GridData->InitializeTemperatureFieldForH2Shield() == FAIL) {
	    ENZO_FAIL("Error in InitializeTemperatureFieldForH2Shield.\n");
	  }
    }

#if USE_H2II_LOOKUP	
    /* Initialize Temperature Field for H2II cross section approximation */
    if(RadiativeTransferH2IIDiss == TRUE) {
      for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) 
	  if (Temp->GridData->InitializeTemperatureFieldForH2IICrossSection() == FAIL) {
	    ENZO_FAIL("Error in InitializeTemperatureFieldForH2II.\n");
	  }
    }
#endif
    
    for (i = 0; i < 4; i++)
      EscapedPhotonCount[i] = 0.0;

    for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++)
      FieldsToInterpolate[i] = FALSE;

    if (NumberOfSources == 0) {
      PhotonTime += dtPhoton;
      continue;
    }    

    /* Create tree that clusters the sources if requested.  While
       creating tree (type SuperSource), compute position of the super
       source in each leaf. */

    START_PERF();
    if (RadiativeTransferSourceClustering == TRUE) {
      CreateSourceClusteringTree(0, NULL, LevelArray);
    }
    END_PERF(1);

    // first identify sources and let them radiate 
    RS = GlobalRadiationSources->NextSource;
    TempGridList = RS_GridList->NextGrid;

    START_PERF(); 
    while (RS != NULL) {
      if (TempGridList->ThisGrid != NULL)
	TempGridList->ThisGrid->Shine(RS);
      TempGridList = TempGridList->NextGrid;
      RS = RS->NextSource;
    }    // while still sources 
    END_PERF(2);

#ifdef USE_MPI
    if (RadiativeTransferInterpolateField)
      CommunicationAllReduceValues(FieldsToInterpolate,
				   MAX_NUMBER_OF_BARYON_FIELDS, MPI_MAX);
#endif /* USE_MPI */  

    /* Initialize interpolated radiation fields */  

    //  if (RadiativeTransferInterpolateField)
    //    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    //      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
    //	if (Temp->GridData->AllocateInterpolatedRadiation() == FAIL) {
    //	  ENZO_FAIL("Error in grid->AllocateInterpolatedRadiation.\n");
    //	}


    /* Evolve all photons by fixed timestep. */
  
    ListOfPhotonsToMove *PhotonsToMove = new ListOfPhotonsToMove;
    PhotonsToMove->NextPackageToMove = NULL;

    int keep_transporting = 1;
    int iteration = 0;
    int local_work = 0;
    int globally_terminated = 0;
#ifdef USE_MPI
    int consensus_active = 0;  // 0=none, 1=barrier, 2=iallreduce pending
    int outstanding_receives = 0;
    int outstanding_sends = 0;
    Eint32 local_done_val = 0;   // MPI_Iallreduce send buffer (must persist)
    Eint32 global_done_val = 0;  // MPI_Iallreduce recv buffer (must persist)
    PH_WorkReceived = 0;
#endif

    HierarchyEntry **Temp0;
    int nGrids0 = GenerateGridArray(LevelArray, 0, &Temp0);
    grid **Grids0 = new grid*[nGrids0];
    for (i = 0; i < nGrids0; i++)
      Grids0[i] = Temp0[i]->GridData;

    /* Initialize nonblocking communication */

    START_PERF();
    InitializePhotonCommunication();
    END_PERF(3);

    /* Transport the rays! */

    PrintMemoryUsage("EvolvePhotons -- before loop");

    while (!globally_terminated && iteration++ < MAX_ITERATIONS) {
      if (DEBUG_RT) {
        printf("P%d: EvolvePhotons loop iteration %d starting\n", MyProcessorNumber, iteration);
        fflush(stdout);
      }

#ifdef USE_MPI
      // Step 1: Check completed Count & Data Receives (always non-blocking here)
      CommunicationReceiverPhotons(LevelArray, false);
      
      // Calculate local work and MPI status
      local_work = 0;
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++) {
        for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) {
          local_work += Temp->GridData->ReturnPhotonPackagePointer()->numPackages;
        }
      }

      outstanding_receives = PH_CommunicationReceiveIndex;
      outstanding_sends = CommunicationBufferedSendActiveCount();
      if (DEBUG_RT) {
        printf("P%d: Iteration %d status: local_work=%d, active_recv=%d, active_send=%d, consensus_active=%d\n",
               MyProcessorNumber, iteration, local_work, outstanding_receives, outstanding_sends, consensus_active);
        fflush(stdout);
      }
#else
      local_work = 0;
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++) {
        for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) {
          local_work += Temp->GridData->ReturnPhotonPackagePointer()->numPackages;
        }
      }
#endif

      // Step 2: Trace locally available rays in grids
      int pre_trace_photons = local_work;
      int has_photons_to_move = 0;
      if (local_work > 0) {
        START_PERF();
        TIMER_START("RayTracing");
        PhotonsToMove->NextPackageToMove = NULL;

        for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0 ; lvl--) {
          for (Temp = LevelArray[lvl], GridNum = 0; Temp; Temp = Temp->NextGridThisLevel, GridNum++) {
            if (Temp->GridHierarchyEntry->ParentGrid != NULL) 
              Helper = Temp->GridHierarchyEntry->ParentGrid->GridData;
            else
              Helper = NULL;

#ifdef BITWISE_IDENTICALITY
            Temp->GridData->PhotonSortLinkedLists();
#endif
            Temp->GridData->TransportPhotonPackages
              (lvl, level, &PhotonsToMove, GridNum, Grids0, nGrids0, Helper, 
               Temp->GridData);
          }
        }
        TIMER_STOP("RayTracing");
        END_PERF(4);
      }

      // Check if any photons crossed grid boundaries
      has_photons_to_move = (PhotonsToMove->NextPackageToMove != NULL) ? 1 : 0;

      // Step 3: Pack and send boundary rays non-blocking
      int temp_keep = 0;
      if (has_photons_to_move) {
        START_PERF();
        TIMER_START("RayCommunication");
        CommunicationTransferPhotons(LevelArray, &PhotonsToMove, NULL, temp_keep);
        TIMER_STOP("RayCommunication");
        END_PERF(5);
      }

      // Step 4: Check for completed sends
#ifdef USE_MPI
      CommunicationBufferPurge();
      outstanding_sends = CommunicationBufferedSendActiveCount();
#endif

      // Recount photons after tracing + sending
      int post_trace_photons = 0;
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
        for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
          post_trace_photons += Temp->GridData->ReturnPhotonPackagePointer()->numPackages;

      // Determine whether transport made progress this iteration.
      // "Done" photons that stay in grids (already traced for this
      // timestep) don't constitute active work — matching the original
      // blocking code's keep_transporting semantics.
      int transport_progress = has_photons_to_move ||
                               (post_trace_photons < pre_trace_photons);

      // When all photons have been traced, merge paused packages
      int nmerges = 0;
      if (RadiativeTransferSourceClustering && !transport_progress) {
        for (lvl = MAX_DEPTH_OF_HIERARCHY-1; lvl >= 0; lvl--) {
          for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel) {
            nmerges += Temp->GridData->MergePausedPhotonPackages();
          }
        }
        if (nmerges > 0) transport_progress = 1;
      }

      // Step 5: Distributed termination check & Consensus
      //
      // Three-phase non-blocking consensus:
      //   Phase 0 (consensus_active==0): Enter Ibarrier when locally idle
      //   Phase 1 (consensus_active==1): Barrier pending; continue processing
      //   Phase 1→2 transition: Barrier completed; drain eager messages, then
      //                         start MPI_Iallreduce (non-blocking!) with
      //                         current local state.  This transition happens
      //                         REGARDLESS of local idle status, avoiding the
      //                         deadlock where blocking MPI_Allreduce prevented
      //                         processing of pending rendezvous Irecvs.
      //   Phase 2 (consensus_active==2): Iallreduce pending; continue processing
      //   Phase 2 completion: If global_done==0, terminate.  Otherwise reset
      //                       to phase 0 for another round.
      //
#ifdef USE_MPI
      int current_local_idle = (!transport_progress && outstanding_receives == 0 && outstanding_sends == 0) ? 1 : 0;

      // --- Barrier→Iallreduce transition (runs even if NOT idle) ---
      if (consensus_active == 1 && PH_ConsensusRequest == MPI_REQUEST_NULL) {
        // Barrier completed (detected by CommunicationReceiverPhotons or MPI_Test).
        // Drain eager messages that may have been delivered during the barrier.
        CommunicationReceiverPhotons(LevelArray, false);  // drain pass 1: count→post data Irecv
        CommunicationReceiverPhotons(LevelArray, false);  // drain pass 2: data Irecv completion

        // Recompute state after drain
        outstanding_receives = PH_CommunicationReceiveIndex;
        outstanding_sends = CommunicationBufferedSendActiveCount();

        // Contribute current state to non-blocking allreduce.
        // "done" (0) only if truly idle AND no work received during the barrier.
        // Note: we don't check numPackages here — done photons sitting in
        // grids are not active work (they're already traced for this timestep).
        local_done_val = (outstanding_sends == 0 &&
                          outstanding_receives == 0 && PH_WorkReceived == 0) ? 0 : 1;
        global_done_val = 1; // safe default (don't terminate)
        MPI_Iallreduce(&local_done_val, &global_done_val, 1, MPI_INT, MPI_MAX,
                       MPI_COMM_WORLD, &PH_ConsensusRequest);
        consensus_active = 2;
        PH_WorkReceived = 0;

        if (DEBUG_RT) {
          printf("P%d: Barrier completed. Started Iallreduce with local_done=%d "
                 "(lw=%d os=%d or=%d)\n",
                 MyProcessorNumber, local_done_val, local_work,
                 outstanding_sends, outstanding_receives);
          fflush(stdout);
        }
      }

      // --- Iallreduce completion check (runs even if NOT idle) ---
      if (consensus_active == 2 && PH_ConsensusRequest == MPI_REQUEST_NULL) {
        // Iallreduce completed (detected by CommunicationReceiverPhotons
        // setting PH_ConsensusRequest = MPI_REQUEST_NULL, or by direct MPI_Test).
        if (DEBUG_RT) {
          printf("P%d: Iallreduce completed: global_done=%d\n",
                 MyProcessorNumber, global_done_val);
          fflush(stdout);
        }
        if (global_done_val == 0) {
          // All ranks reported done.  Re-verify local state to guard
          // against messages that arrived during the Iallreduce.
          outstanding_receives = PH_CommunicationReceiveIndex;
          outstanding_sends = CommunicationBufferedSendActiveCount();
          if (outstanding_sends == 0 &&
              outstanding_receives == 0 && PH_WorkReceived == 0) {
            globally_terminated = 1;
            if (DEBUG_RT) {
              printf("P%d: Globally terminated!\n", MyProcessorNumber);
              fflush(stdout);
            }
          } else {
            // Race: new work arrived during Iallreduce.  Reset.
            if (DEBUG_RT) {
              printf("P%d: Iallreduce said done but local state dirty "
                     "(lw=%d os=%d or=%d wr=%d). Resetting.\n",
                     MyProcessorNumber, local_work, outstanding_sends,
                     outstanding_receives, PH_WorkReceived);
              fflush(stdout);
            }
            consensus_active = 0;
            PH_ConsensusRequest = MPI_REQUEST_NULL;
          }
        } else {
          // Some rank had work — reset and try again when idle
          consensus_active = 0;
          PH_ConsensusRequest = MPI_REQUEST_NULL;
        }
      }

      // --- Phase 0: Enter new barrier when locally idle ---
      if (current_local_idle && consensus_active == 0 && !globally_terminated) {
        if (DEBUG_RT) {
          printf("P%d: Iteration %d is locally idle. Entering Ibarrier.\n",
                 MyProcessorNumber, iteration);
          fflush(stdout);
        }
        PH_WorkReceived = 0;
        MPI_Ibarrier(MPI_COMM_WORLD, &PH_ConsensusRequest);
        consensus_active = 1;
      }

      // --- Phase 1 or 2: Explicit MPI_Test to drive progress ---
      if (consensus_active >= 1 && PH_ConsensusRequest != MPI_REQUEST_NULL) {
        Eint32 test_flag = 0;
        MPI_Test(&PH_ConsensusRequest, &test_flag, MPI_STATUS_IGNORE);
        // If completed, PH_ConsensusRequest is now MPI_REQUEST_NULL.
        // The transition logic above will handle it on the next iteration.
      }

      // Step 6: CPU Yielding / Waitsome if idle and waiting.
      // Only block when there IS a pending consensus request to wait on;
      // blocking with PH_ConsensusRequest==NULL deadlocks because
      // Waitsome would have only count Irecvs (nothing to wake it).
      if (!globally_terminated && !transport_progress && outstanding_sends == 0 &&
          nmerges == 0 && consensus_active >= 1 &&
          PH_ConsensusRequest != MPI_REQUEST_NULL) {
        CommunicationReceiverPhotons(LevelArray, true);
      }
#else
      if (!transport_progress && nmerges == 0) {
        globally_terminated = 1;
      }
#endif
    }

    if (DEBUG_RT) {
      printf("P%d: EvolvePhotons loop finished. Calling FinalizePhotonCommunication()\n", MyProcessorNumber);
      fflush(stdout);
    }
    FinalizePhotonCommunication();
    if (DEBUG_RT) {
      printf("P%d: EvolvePhotons finished FinalizePhotonCommunication()\n", MyProcessorNumber);
      fflush(stdout);
    }

    /* Move all finished photon packages back to their original place,
       PhotonPackages.  For the adaptive timestep, we don't carryover
       any photons to the next timestep. */

    START_PERF();
    if (RadiativeTransferAdaptiveTimestep)  
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  Temp->GridData->DeletePhotonPackages();  
    else
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  Temp->GridData->MoveFinishedPhotonsBack();
    END_PERF(7);
    PrintMemoryUsage("EvolvePhotons -- deleted photons");

    /* If we're keeping track of photon escape fractions on multiple
       processors, collect photon counts from all processors */

    FILE *fptr;

    if (RadiativeTransferPhotonEscapeRadius > 0) {
      CommunicationSumValues(EscapedPhotonCount, 4);
      if (MyProcessorNumber == ROOT_PROCESSOR) {

	/* Open f_esc file for writing */

	if (TotalEscapedPhotonCount[0] <= 0) {
	  if ((fptr = fopen(PhotonEscapeFilename, "w")) == NULL) {
	    ENZO_VFAIL("Error opening file %s\n", PhotonEscapeFilename)
	  }
	  fprintf(fptr, 
		  "# Time TotalPhotons fesc(0.5rvir) fesc(rvir) fesc(2rvir)\n");
	} else {
	  if ((fptr = fopen(PhotonEscapeFilename, "a")) == NULL) {
	    ENZO_VFAIL("Error opening file %s\n", PhotonEscapeFilename)
	  }
	}

	for (i = 0; i < 4; i++)
	  TotalEscapedPhotonCount[i] += EscapedPhotonCount[i];

	fprintf(fptr, "%"GOUTSYM" %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", 
		PhotonTime, TotalEscapedPhotonCount[0],
		TotalEscapedPhotonCount[1] / TotalEscapedPhotonCount[0], 
		TotalEscapedPhotonCount[2] / TotalEscapedPhotonCount[0], 
		TotalEscapedPhotonCount[3] / TotalEscapedPhotonCount[0]);

	fclose(fptr);

      } // ENDIF ROOT_PROCESSOR

    } // ENDIF RTPhotonEscapeRadius

    PhotonTime += dtPhoton;

    delete PhotonsToMove;
    delete [] Grids0;
    delete [] Temp0;

    /* Sync photon counts after ray tracing */

    START_PERF();
    CommunicationSyncNumberOfPhotons(LevelArray);
    END_PERF(11);

    /* Delete baryon fields on temporary "fake" grid on
       ProcessorNumber and revert it back to OriginalProcessorNumber,
       which was saved in CommunicationLoadBalancePhotonGrids.  Photon
       packages must be moved back, too. */

    if (RadiativeTransferLoadBalance && NumberOfSources > 0)
      RadiativeTransferLoadBalanceRevert(Grids, nGrids);

    /************************************************************************/
    /********************* Coupled rate & energy solver *********************/
    /************************************************************************/

    int debug_store = debug;
    debug = FALSE;

    // Divide the photo-ionization and photo-heating rates by the
    // number of particles (rho * dx^3)
    START_PERF();
    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	if (Temp->GridData->RadiationPresent() == TRUE)
	  Temp->GridData->FinalizeRadiationFields();
    END_PERF(8);

    /* Set the optically-thin H2 dissociation rates */

    START_PERF();
    if (RadiativeTransferOpticallyThinH2)
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  Temp->GridData->AddH2Dissociation(AllStars, NumberOfSources);
    END_PERF(10);

    START_PERF();
    if (RadiativeTransferCoupledRateSolver)
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  if (Temp->GridData->RadiationPresent() == TRUE) {

	    int RTCoupledSolverIntermediateStep = TRUE;

#ifdef USE_GRACKLE
            if (grackle_data->use_grackle == TRUE){
              grackle_data->radiative_transfer_intermediate_step = (Eint32) RTCoupledSolverIntermediateStep;
              if (Temp->GridData->GrackleWrapper() == FAIL){
                ENZO_FAIL("Error in GrackleWrapper.\n");
              }
            } else
#endif // USE_GRACKLE
	    {
	      Temp->GridData->SolveRateAndCoolEquations(RTCoupledSolverIntermediateStep);
	    }

	  } /* ENDIF radiation */
    END_PERF(9);

    /* Clean up temperature field */

    if (RadiationXRayComptonHeating)
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  if (Temp->GridData->FinalizeTemperatureFieldForComptonHeating() == FAIL) {  
	    ENZO_FAIL("Error in FinalizeTemperatureFieldForComptonHeating.\n");
	  }	
    if (RadiativeTransferH2ShieldType == 1 || ProblemType == 50)
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  if (Temp->GridData->FinalizeTemperatureFieldForH2Shield() == FAIL) {  
	    ENZO_FAIL("Error in FinalizeTemperatureFieldForH2Shield.\n");
	  }
#if USE_H2II_LOOKUP	
    if (RadiativeTransferH2IIDiss == TRUE) {
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  if (Temp->GridData->FinalizeTemperatureFieldForH2IICrossSection() == FAIL) {  
	    ENZO_FAIL("Error in FinalizeTemperatureFieldForH2Shield.\n");
	  }	
    }
#endif
    debug = debug_store;

    /* If we're using the HII restricted timestep, get the global
       maximum kph in I-fronts. */

    START_PERF();
    if (RadiativeTransferHIIRestrictedTimestep) {
      float LocalMaximumkph = -1e20;
      for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY-1; lvl++)
	for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	  LocalMaximumkph = max(LocalMaximumkph,
				Temp->GridData->ReturnMaximumkphIfront());
      LocalMaximumkph = CommunicationMaxValue(LocalMaximumkph);
      MetaData->GlobalMaximumkphIfront = LocalMaximumkph;
    }
    END_PERF(12);

#ifdef DEBUG
    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->ErrorCheckPhotonNumber(lvl);
#endif

    if (!LoopTime)
      break;
    
    FirstTime = false;

  } // ENDWHILE GridTime >= PhotonTime

  /* Cleanup photon memory pool if we're deleting all photons between
     timesteps, i.e. no need to save photons */

  START_PERF();
#ifdef MEMORY_POOL
  const int PhotonMemorySize = MEMORY_POOL_SIZE;
  int PhotonSize = sizeof(PhotonPackageEntry);
  if (RadiativeTransferAdaptiveTimestep) {
    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->DeletePhotonPackages(TRUE);
    delete PhotonMemoryPool;
    PhotonMemoryPool = new MPool::MemoryPool(PhotonMemorySize*PhotonSize,
					     PhotonSize,
					     PhotonMemorySize*PhotonSize/4);
    for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
      for (Temp = LevelArray[lvl]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->InitializePhotonPackages();
  }
#endif
  END_PERF(13);

#ifdef REPORT_PERF
  if (debug) printf("EvolvePhotons: total time = %g\n", ReturnWallTime()-ep0);
  for (int i = 0; i < 1; i++) {  // Only report on ROOT, not all
    CommunicationBarrier();
    if (MyProcessorNumber == i) {

      printf("P%d:", MyProcessorNumber);
      fpcol(PerfCounter, 14, 14, stdout);
      fflush(stdout);
    }
  }
#endif

  /* Double check and clean up photon sources that have expired */

  if (NumberOfSources > 0){
    RS = GlobalRadiationSources->NextSource;

    /* Clean up photon sources if time is up */
    while (RS != NULL){
      if ( ((RS->CreationTime + RS->LifeTime) < PhotonTime) && LoopTime == TRUE &&
             DeleteSources) {

        if (debug) {
          fprintf(stdout, "\nEvolvePhotons: Deleted Source on lifetime limit \n");
          fprintf(stdout, "EvolvePhotons:  %"GSYM" %"GSYM" %"GSYM" \n",
                                             RS->CreationTime, RS->LifeTime, PhotonTime);
        }

        RS = DeleteRadiationSource(RS);
        NumberOfSources--;

      }  else {
        RS = RS->NextSource;
      }

    } // end while RS loop
  } // end source clean up

  /* Delete GridList and Reset Grid IDs in static sources */

  ThinGridList *Orphan;
  TempGridList = RS_GridList;
  while (TempGridList != NULL) {
    Orphan = TempGridList;
    TempGridList = TempGridList->NextGrid;
    if (Orphan != NULL) delete Orphan;
  }

  if (ProblemType == 50) {
    for (RS = GlobalRadiationSources->NextSource; RS; RS = RS->NextSource) {
      RS->GridID = INT_UNDEFINED;
      RS->GridLevel = INT_UNDEFINED;
    }
  }

  /* Delete grid lists */

  for (lvl = 0; lvl < MAX_DEPTH_OF_HIERARCHY; lvl++)
    if (nGrids[lvl] > 0) delete [] Grids[lvl];

  LCAPERF_STOP("EvolvePhotons");
  TIMER_STOP("EvolvePhotons");
  return SUCCESS;

}
