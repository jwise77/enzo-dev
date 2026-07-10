#define DEBUG 0
/***********************************************************************
/
/  GRID CLASS (SEND PHOTONS FROM REAL GRID TO 'FAKE' (REPLICATED) GRID)
/
/  written by: John H. Wise
/  date:       November, 2005
/  modified1:
/
/  PURPOSE: 
/
/    NOTE: We assume all the from grids are at the same level!
/    NOTE: Modelled after GBs Grid_CommunicationSendParticles.C
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

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
#include "communication.h"
#include "CommunicationUtilities.h"
#include "GroupPhotonList.h"

#ifdef USE_MPI
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, 
                              int Target, int Tag, MPI_Comm CommWorld, 
			      int BufferSize);
static int FirstTimeCalled = TRUE;
static MPI_Datatype PhotonBufferType;
#endif /* USE_MPI */

void my_exit(int status);
int FindSuperSource(PhotonPackageEntry **PP, int &LeafID, 
		    int SearchNewTree = TRUE);
PhotonPackageEntry* DeletePhotonPackage(PhotonPackageEntry *PP);
void InsertPhotonAfter(PhotonPackageEntry * &Node, PhotonPackageEntry * &NewNode);

int grid::CommunicationSendPhotonPackages(grid *ToGrid, int ToProcessor,
					  int ToNumber, int FromNumber, 
					  PhotonPackageSoA *ToPP)
{

  int index, dim, temp_int;

  if (CommunicationShouldExit(ProcessorNumber, ToProcessor))
    return SUCCESS;

  if (FromNumber == 0) 
    return SUCCESS;

  /* Allocate memory */

  PhotonBuffer *buffer = NULL;
#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = (PhotonBuffer *) CommunicationReceiveBuffer[CommunicationReceiveIndex];
  else
#endif /* USE_MPI */
    buffer = new PhotonBuffer[FromNumber];

  if (DEBUG)
    printf("SendPhotonPackages(%"ISYM"): Sending %"ISYM" photon packages "
	   "from P%"ISYM"->P%"ISYM".\n",
	   MyProcessorNumber, FromNumber, ProcessorNumber, ToProcessor);

  /* If this is from processor, pack photons */

  if (MyProcessorNumber == ProcessorNumber) {
    for (index = 0; index < PhotonPackages.numPackages; index++) {
      buffer[index].Photons		 = PhotonPackages.Flux[index];
      buffer[index].Type		 = PhotonPackages.Type[index];
      buffer[index].Energy		 = PhotonPackages.Energy[index];
      buffer[index].EmissionTimeInterval = PhotonPackages.TimeInterval[index];
      buffer[index].EmissionTime	 = PhotonPackages.EmissionTime[index];
      buffer[index].CurrentTime          = PhotonPackages.CurrentTime[index];
      buffer[index].ColumnDensity        = PhotonPackages.ColumnDensity[index];
      buffer[index].CrossSection         = PhotonPackages.CrossSection[index];
      buffer[index].Radius		 = PhotonPackages.Radius[index];
      buffer[index].ipix		 = PhotonPackages.PixelNum[index];
      buffer[index].level		 = PhotonPackages.Level[index];
      for (dim = 0; dim < GridRank; dim++)
	buffer[index].SourcePosition[dim] = (dim == 0) ? PhotonPackages.SourceX[index] :
                                            ((dim == 1) ? PhotonPackages.SourceY[index] :
                                             PhotonPackages.SourceZ[index]);
      buffer[index].SourcePositionDiff   = PhotonPackages.SourcePositionDiff[index];

      if (PhotonPackages.CurrentSource[index] != NULL)
	buffer[index].SuperSourceID = PhotonPackages.CurrentSource[index]->LeafID;
      else
	buffer[index].SuperSourceID = -1;
      
      if (PhotonPackages.CurrentTime[index] < 0 || PhotonPackages.CurrentTime[index] > 1e10) {
	ENZO_VFAIL("CTPhotons[0][P%"ISYM"->P%"ISYM"]: "
		"(%"ISYM" of %"ISYM") Bad photon time %"GSYM"\n",
		ProcessorNumber, ToProcessor, index, NumberOfPhotonPackages, 
		PhotonPackages.CurrentTime[index])
      }
    }

    if (DEBUG)
      printf("CommSendPhotons(P%"ISYM"): Counted %"ISYM" photons.\n", MyProcessorNumber,
	     index);

    /* Now that we're done packing the photons, delete them */
    PhotonPackages.free_arrays();

    /* Check if we packed all of the photons */

    if (index != FromNumber) {
      fprintf(stdout, "CommSendPhotons WARNING: Counted %"ISYM" photon packages, but"
	      " FromNumber = %"ISYM"\n", index, FromNumber);
      FromNumber = min(index, FromNumber);
      fprintf(stdout, "CommSendPhotons: Correcting FromNumber to %"ISYM"\n", 
	      FromNumber);
      //ENZO_FAIL("Photon package mismatch!\n");
    }

  } /* ENDIF PackPhotons */

#ifdef USE_MPI

  /* Send buffer if the processor numbers aren't identical */

  if (ProcessorNumber != ToProcessor) {

    MPI_Status status;
    MPI_Datatype DataTypeByte = MPI_BYTE;
    MPI_Arg PhotonBufferSize;
    MPI_Arg Count = FromNumber;
    MPI_Arg Source = ProcessorNumber;
    MPI_Arg Dest = ToProcessor;
    MPI_Arg stat;

    if (FirstTimeCalled) {
      PhotonBufferSize = sizeof(PhotonBuffer);
      //  fprintf(stderr, "Size of ParticleMoveList %"ISYM"\n", Count);
      stat = MPI_Type_contiguous(PhotonBufferSize, DataTypeByte, &PhotonBufferType);
      stat |= MPI_Type_commit(&PhotonBufferType);
      if (stat != MPI_SUCCESS) my_exit(EXIT_FAILURE);
      FirstTimeCalled = FALSE;
    }

    if (MyProcessorNumber == ProcessorNumber) {
      if (DEBUG)
	printf("PhotonSend(P%"ISYM"): Sending %"ISYM" photons to processor %"ISYM".\n",
	       MyProcessorNumber, FromNumber, ToProcessor);
      CommunicationBufferedSend(buffer, Count, PhotonBufferType, Dest,
				MPI_PHOTON_TAG, MPI_COMM_WORLD, BUFFER_IN_PLACE);
    }

    if (MyProcessorNumber == ToProcessor) {

      if (DEBUG) 
	printf("PhotonSend(P%"ISYM"): Receiving %"ISYM" photons from processor %"ISYM".\n",
	       MyProcessorNumber, FromNumber, ProcessorNumber);

      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
	MPI_Irecv(buffer, Count, PhotonBufferType, Source, MPI_PHOTON_TAG,
		  MPI_COMM_WORLD,
		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);

	CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
	CommunicationReceiveGridTwo[CommunicationReceiveIndex] = ToGrid;
	CommunicationReceiveCallType[CommunicationReceiveIndex] = 15;
	CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] = ToNumber;
	CommunicationReceiveArgumentInt[1][CommunicationReceiveIndex] = FromNumber;
	CommunicationReceiveBuffer[CommunicationReceiveIndex] = (float *) buffer;
	CommunicationReceiveDependsOn[CommunicationReceiveIndex] = 
	  CommunicationReceiveCurrentDependsOn;
	CommunicationReceiveIndex++;


      } // ENDIF post receive

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
	if (MPI_Recv(buffer, Count, PhotonBufferType, Source,
		     MPI_PHOTON_TAG, MPI_COMM_WORLD, &status) != MPI_SUCCESS) {
	  fprintf(stderr, "P(%"ISYM"): MPI_Recv error %"ISYM"\n", MyProcessorNumber,
		  status.MPI_ERROR);
	  fprintf(stderr, "P(%"ISYM"): TransferSize = %"ISYM" ProcessorNumber = %"ISYM"\n", 
		  MyProcessorNumber, Count*sizeof(PhotonBuffer), ProcessorNumber);
	  char errstr[MPI_MAX_ERROR_STRING];
	  Eint32 errlen;
	  MPI_Error_string(status.MPI_ERROR, errstr, &errlen);
	  ENZO_VFAIL("MPI Error %s\n",errstr)
	}

    } /* ENDIF (MyProcessorNumber == ToProcessor) */

  } /* ENDIF (ProcessorNumber != ToProcessor) */
#endif /* USE_MPI */

  /* If this is the to processor, unpack fields */

  if (MyProcessorNumber == ToProcessor && 
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {

    for (index = 0; index < FromNumber; index++) {

      ToPP->append(buffer[index].Photons, buffer[index].Type,
                   buffer[index].Energy, buffer[index].CrossSection,
                   buffer[index].EmissionTimeInterval, buffer[index].EmissionTime,
                   buffer[index].CurrentTime, buffer[index].Radius,
                   buffer[index].ColumnDensity, buffer[index].ipix,
                   buffer[index].level, buffer[index].SourcePosition[0],
                   buffer[index].SourcePosition[1], buffer[index].SourcePosition[2],
                   buffer[index].SourcePositionDiff, NULL);

      int idx = ToPP->numPackages - 1;

      if (buffer[index].CurrentTime < 0 || buffer[index].CurrentTime > 1e10) {
	ENZO_VFAIL("CTPhotons[1][P%"ISYM"->P%"ISYM"]: "
		"(%"ISYM" of %"ISYM") Bad photon time %"GSYM"\n",
		ProcessorNumber, ToProcessor, index, FromNumber, 
		buffer[index].CurrentTime)
      }

      if (RadiativeTransferSourceClustering) {
	PhotonPackageEntry tempPP;
	int leafID = buffer[index].SuperSourceID;
	PhotonPackageEntry *tempPPPtr = &tempPP;
	if (FindSuperSource(&tempPPPtr, leafID) == FAIL) {
	  ENZO_FAIL("Error in FindSuperSource.\n");
	}
	ToPP->CurrentSource[idx] = tempPP.CurrentSource;
      } else
	ToPP->CurrentSource[idx] = NULL;

    } /* ENDFOR index */

    /* Only delete the buffer if we're in receive mode (in send mode
       it will be deleted by CommunicationBufferedSend and if we're in
       post-receive mode then it will be deleted when we get to
       receive-mode). */

    delete [] buffer;

  }  /* ENDIF (MyProcessorNumber == ToProcessor) */

  return SUCCESS;

}
