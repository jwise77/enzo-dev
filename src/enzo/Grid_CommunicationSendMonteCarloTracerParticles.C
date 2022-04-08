/***********************************************************************
/
/  GRID CLASS (SEND MONTE CARLO TRACER PARTICLES FROM REAL GRID TO 'FAKE'
/             (REPLICATED) GRID)
/
/  written by: Corey Brummel-Smith
/  date:       January, 2022
/  modified1:  
/
/  NOTES:  Adapted from grid::CommunicationSendMonteCarloTracerParticles().
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
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
#include "communication.h"
#include "CommunicationUtilities.h"

/* function prototypes */

#ifdef USE_MPI
static int FirstTimeCalled = TRUE;
static MPI_Datatype MPI_MCTP; 
int CommunicationBufferedSend(void *buffer, int size, MPI_Datatype Type, 
                              int Target, int Tag, MPI_Comm CommWorld, 
			      int BufferSize);
#endif /* USE_MPI */
MonteCarloTracerParticle* MonteCarloTracerParticleBufferToList(MonteCarloTracerParticleBuffer buffer);
void InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticle * &Node, MonteCarloTracerParticle * &NewNode);
void DeleteMonteCarloTracerParticleList(MonteCarloTracerParticle * &Node);

/* Send particle from this grid to ToGrid on processor ToProcessor. */

int grid::CommunicationSendMonteCarloTracerParticles(grid *ToGrid, int ToProcessor)
{

  int i, j, k, n,dim, index, TransferSize;
  int index_ijk[MAX_DIMENSION];
  float DomainWidth[MAX_DIMENSION], DomainWidthInv[MAX_DIMENSION];  
  MonteCarloTracerParticleBuffer *buffer = NULL;
  MonteCarloTracerParticle *mctp;
  grid *DestinationGrid;

  if (CommunicationShouldExit(ProcessorNumber, ToProcessor))
    return SUCCESS;

  TransferSize = NumberOfMonteCarloTracerParticles;
  printf("\nproc%d: CommSendMCTP: TransferSize (NumberOfMonteCarloTracerParticles) %d", MyProcessorNumber, TransferSize);

  if (TransferSize == 0)
    return SUCCESS;

  /* Allocate buffer in ToProcessor.  This is automatically done in
     MonteCarloTracerParticleListToBuffer in the local processor. */

#ifdef USE_MPI
  if (CommunicationDirection == COMMUNICATION_RECEIVE)
    buffer = (MonteCarloTracerParticleBuffer*) CommunicationReceiveBuffer[CommunicationReceiveIndex];
  else
#endif
    buffer = new MonteCarloTracerParticleBuffer[TransferSize];

  /* If this is the from processor, fill the buffer and delete mc tracer particles from this grid. */

  if (MyProcessorNumber == ProcessorNumber) {
    printf("\nproc%d: CommSendMC: package buffer and delete", MyProcessorNumber);      
    this->MonteCarloTracerParticles[0]->MonteCarloTracerParticleListToBuffer(buffer, NumberOfMonteCarloTracerParticles);
    DeleteMonteCarloTracerParticleList(this->MonteCarloTracerParticles[0]);
    // *** Add flag to not delete particles when needed (e.g. combine grids) ****
  }
    
  /* Send buffer. */

#ifdef USE_MPI

  /* only send if processor numbers are not identical */

  if (ProcessorNumber != ToProcessor) {

    if (FirstTimeCalled) {
      MPI_Type_contiguous(sizeof(MonteCarloTracerParticleBuffer), MPI_BYTE, &MPI_MCTP);
      MPI_Type_commit(&MPI_MCTP);
      FirstTimeCalled = FALSE;
    }

    MPI_Status status;
    MPI_Arg PCount, Count = TransferSize;
    MPI_Arg Source = ProcessorNumber;
    MPI_Arg Dest = ToProcessor;
    MPI_Arg stat;

#ifdef MPI_INSTRUMENTATION
    starttime = MPI_Wtime();
#endif
    if (MyProcessorNumber == ProcessorNumber){
      printf("\nproc%d: CommSendMCTP: Sending %d MCTPs from proc%d to proc%d\n", MyProcessorNumber, TransferSize, ProcessorNumber, ToProcessor);
      CommunicationBufferedSend(buffer, Count, MPI_MCTP, 
				Dest, MPI_SENDMCTP_TAG, MPI_COMM_WORLD, 
				BUFFER_IN_PLACE);
    }

    if (MyProcessorNumber == ToProcessor) {
      printf("\nproc%d: CommSendMCTP: Receiving %d MCTPs on proc%d from proc%d\n", MyProcessorNumber, TransferSize, ToProcessor, ProcessorNumber);
      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
      	MPI_Irecv(buffer, Count, MPI_MCTP, Source,
      		  MPI_SENDMCTP_TAG, MPI_COMM_WORLD,
      		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);

      	CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
      	CommunicationReceiveGridTwo[CommunicationReceiveIndex] = ToGrid;
      	CommunicationReceiveCallType[CommunicationReceiveIndex] = 23;
      	CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] = NumberOfMonteCarloTracerParticles;

      	CommunicationReceiveBuffer[CommunicationReceiveIndex] = (float *) buffer;
      	CommunicationReceiveDependsOn[CommunicationReceiveIndex] = CommunicationReceiveCurrentDependsOn;
      	CommunicationReceiveIndex++;
      }

      if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
	MPI_Recv(buffer, Count, MPI_MCTP, Source,
		 MPI_SENDMCTP_TAG, MPI_COMM_WORLD, &status);

    } // ENDIF MyProcessorNumber == ToProcessor

#ifdef MPI_INSTRUMENTATION
    /* Zhiling Lan's instrumented part */
    endtime = MPI_Wtime();
    timer[7] += endtime-starttime;
    counter[7] ++;
    timer[8] += double(TransferSize);
    timer[28] += double(TransferSize*TransferSize);
    timer[27] += (endtime-starttime)*(endtime-starttime);
#endif /* MPI_INSTRUMENTATION */
  
  } // end: if (ProcessorNumber != ToProcessor)

#endif /* USE_MPI */

  /* If this is the to processor, unpack fields. */

  if (MyProcessorNumber == ToProcessor &&
      (CommunicationDirection == COMMUNICATION_SEND_RECEIVE ||
       CommunicationDirection == COMMUNICATION_RECEIVE)) {

    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
      DomainWidthInv[dim] = 1.0/DomainWidth[dim];    
    }

    if (ToGrid->MonteCarloTracerParticles == NULL) {
      ToGrid->AllocateMonteCarloTracerParticleData();
    }

    for (n = 0; n < TransferSize; n++) {
      mctp = MonteCarloTracerParticleBufferToList(buffer[n]);
      mctp->CurrentGrid = ToGrid;
      
      /* Find which cell this particle belongs in */
      for (dim = 0; dim < GridRank; dim++) {
          index_ijk[dim] = (int) ((NumberOfGhostZones + ToGrid->GridDimension[dim]) * 
                                  (mctp->Position[dim] - DomainLeftEdge[dim]) *
                                  DomainWidthInv[dim]);
      }
      index = GetIndex(index_ijk[0], index_ijk[1], index_ijk[2]); 
      InsertMonteCarloTracerParticleAfter(ToGrid->MonteCarloTracerParticles[index], mctp);
    }
    printf("\nproc%d: CommSendMCTP: Inserted %d MCTPs from buffer into ToGrid", MyProcessorNumber, TransferSize);
    delete [] buffer;
			  
  } // end: if (CommunicationDirection...)

  return SUCCESS;
}

