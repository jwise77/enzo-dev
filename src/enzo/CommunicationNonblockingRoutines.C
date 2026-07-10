/***********************************************************************
/
/  NON-BLOCKING COMMUNICATION ROUTINES FOR PHOTONS
/
/  written by: John H. Wise
/  date:       December, 2007
/  modified1:  July, 2026 (Refactored to modern MPI-3 async consensus)
/
/  PURPOSE: Implements photon communication initialization, finalization,
/           and count message sends.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "GroupPhotonList.h"
#include "PhotonCommunication.h"

#ifdef USE_MPI
int CommunicationBufferPurge(void);
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, int Target,
			      int Tag, MPI_Comm CommWorld, int BufferSize);
#endif /* USE_MPI */

/******************************************************************/

int InitializePhotonCommunication(void)
{
#ifdef USE_MPI
  int proc, i;

  // Initialize MPI datatype on first call
  static int FirstTimeCalled = TRUE;
  if (FirstTimeCalled) {
    MPI_Type_contiguous(sizeof(GroupPhotonList), MPI_BYTE, &MPI_PhotonList);
    MPI_Type_commit(&MPI_PhotonList);
    FirstTimeCalled = FALSE;
  }

  // Initialize data receive slots
  for (i = 0; i < MAX_PH_RECEIVE_BUFFERS; i++) {
    PH_CommunicationReceiveMPI_Request[i] = MPI_REQUEST_NULL;
    PH_CommunicationReceiveBuffer[i] = NULL;
  }
  PH_CommunicationReceiveIndex = 0;
  PH_ConsensusRequest = MPI_REQUEST_NULL;

  // Initialize count receives
  PhotonMessageIndex = 0;
  PhotonMessageMaxIndex = 0;
  for (proc = 0; proc < NumberOfProcessors; proc++) {
    PhotonMessageRequest[proc] = MPI_REQUEST_NULL;
    PhotonMessageBuffer[proc] = 0;
  }

  for (proc = 0; proc < NumberOfProcessors; proc++) {
    if (proc != MyProcessorNumber) {
      MPI_Irecv(&PhotonMessageBuffer[proc], 1, MPI_INT, proc, MPI_NPHOTON_TAG, MPI_COMM_WORLD, &PhotonMessageRequest[proc]);
    }
  }
#endif /* USE_MPI */
  return SUCCESS;
}

/**********************************************************************/

int FinalizePhotonCommunication(void)
{
#ifdef USE_MPI
  int proc, i;

  // Cancel and wait for all outstanding count receives
  for (proc = 0; proc < NumberOfProcessors; proc++) {
    if (PhotonMessageRequest[proc] != MPI_REQUEST_NULL) {
      MPI_Cancel(&PhotonMessageRequest[proc]);
      MPI_Wait(&PhotonMessageRequest[proc], MPI_STATUS_IGNORE);
    }
  }

  // Cancel and wait for all outstanding data receives
  for (i = 0; i < MAX_PH_RECEIVE_BUFFERS; i++) {
    if (PH_CommunicationReceiveMPI_Request[i] != MPI_REQUEST_NULL) {
      MPI_Cancel(&PH_CommunicationReceiveMPI_Request[i]);
      MPI_Wait(&PH_CommunicationReceiveMPI_Request[i], MPI_STATUS_IGNORE);
      GroupPhotonList *RecvBuffer = (GroupPhotonList *) PH_CommunicationReceiveBuffer[i];
      delete [] RecvBuffer;
      PH_CommunicationReceiveBuffer[i] = NULL;
    }
  }
  PH_CommunicationReceiveIndex = 0;

  if (PH_ConsensusRequest != MPI_REQUEST_NULL) {
    MPI_Cancel(&PH_ConsensusRequest);
    MPI_Wait(&PH_ConsensusRequest, MPI_STATUS_IGNORE);
  }

  CommunicationBufferPurge();
#endif /* USE_MPI */
  return SUCCESS;
}

/**********************************************************************/

int CommunicationNumberOfPhotonSends(int *nPhoton, int size)
{
#ifdef USE_MPI
  Eint32 NumberOfMessages, proc;

  for (proc = 0; proc < NumberOfProcessors; proc++) {
    if (proc != MyProcessorNumber) {
      NumberOfMessages = nPhoton[proc] / size;
      if (nPhoton[proc] % size > 0) NumberOfMessages++;
      
      CommunicationBufferedSend(&NumberOfMessages, 1, MPI_INT, proc,
				MPI_NPHOTON_TAG, MPI_COMM_WORLD, sizeof(Eint32));
    }
  }
#endif /* USE_MPI */
  return SUCCESS;
}
