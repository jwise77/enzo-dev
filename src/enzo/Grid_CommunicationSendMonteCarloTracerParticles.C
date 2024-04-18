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

int grid::CommunicationSendMonteCarloTracerParticles(grid *ToGrid, int ToProcessor, int NumberOfParticlesToTransfer, bool DeleteParticles, bool COM_COMBINE)
{

  int i, j, k, n,dim, index, TransferSize;
  int index_ijk[MAX_DIMENSION], ToGridActiveDim[MAX_DIMENSION];
  MonteCarloTracerParticleBuffer *buffer = NULL;
  MonteCarloTracerParticle *mctp;
  grid *DestinationGrid;

  if (CommunicationShouldExit(ProcessorNumber, ToProcessor)){
      printf("\nCOMSEND_EXIT: proc%d, thisGridProc%d, ToGridProc%d, this %p, ToGrid %p",
               MyProcessorNumber, ProcessorNumber, ToProcessor, this, ToGrid); //DEBUG 
    return SUCCESS;
  }

  TransferSize = NumberOfParticlesToTransfer;
  if (COM_COMBINE){
    printf("\nCOMSEND_TRANSFER: proc%d, ToProc%d: CommSendMCTP: TransferSize (NumberOfMonteCarloTracerParticles) %d", MyProcessorNumber, ToProcessor, TransferSize);
    fflush(stdout);   
  }

  if (TransferSize == 0){
    //printf("\nproc%d, ToProc%d: TransferSize == 0 return", MyProcessorNumber, ToProcessor, TransferSize);
    return SUCCESS;
  }

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
    //printf("\nproc%d: CommSendMC: package buffer and delete", MyProcessorNumber);      
    //this->WriteMCTP("ComSend_this_B0");
    this->MonteCarloTracerParticles[0]->MonteCarloTracerParticleListToBuffer(buffer, TransferSize);
    //this->WriteMCTP("ComSend_this_B1");
    if (DeleteParticles)
      DeleteMonteCarloTracerParticleList(this->MonteCarloTracerParticles[0]);

    if (COM_COMBINE) //DEBUG
      printf("\nCOMSEND_CALLED_1: Pack buffer proc%d, thisGridProc%d, ToGridProc%d, this %p, ToGrid %p",
               MyProcessorNumber, ProcessorNumber, ToProcessor, this, ToGrid);
    fflush(stdout);
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
      fflush(stdout);
      CommunicationBufferedSend(buffer, Count, MPI_MCTP, 
				Dest, MPI_SENDMCTP_TAG, MPI_COMM_WORLD, 
				BUFFER_IN_PLACE);
    }

    if (MyProcessorNumber == ToProcessor) {
      printf("\nproc%d: CommSendMCTP: Receiving %d MCTPs on proc%d from proc%d\n", MyProcessorNumber, TransferSize, ToProcessor, ProcessorNumber);
      fflush(stdout);
      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
      	MPI_Irecv(buffer, Count, MPI_MCTP, Source,
      		  MPI_SENDMCTP_TAG, MPI_COMM_WORLD,
      		  CommunicationReceiveMPI_Request+CommunicationReceiveIndex);

      	CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
      	CommunicationReceiveGridTwo[CommunicationReceiveIndex] = ToGrid;
      	CommunicationReceiveCallType[CommunicationReceiveIndex] = 23;
      	CommunicationReceiveArgumentInt[0][CommunicationReceiveIndex] = TransferSize;

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

    // for (dim = 0; dim < MAX_DIMENSION; dim++)
    //   ToGridActiveDim[dim] = ToGrid->GridEndIndex[dim] - ToGrid->GridStartIndex[dim] +1;           

    if (ToGrid->MonteCarloTracerParticles == NULL) {
      // if (COM_COMBINE) 
      //   printf("%s\n", "COMSEND_NULL");
      ToGrid->AllocateMonteCarloTracerParticleData();
    }

    if (COM_COMBINE) //DEBUG
      printf("\nCOMSEND_CALLED_2: proc%d, thisGridProc%d, ToGridProc%d, this %p, ToGrid %p",
               MyProcessorNumber, ProcessorNumber, ToProcessor, this, ToGrid);
    fflush(stdout);

    //ToGrid->WriteMCTP("ComSend_ToGrid_A0"); //DEBUG    

    // printf("COMSEND_INFO: proc%d, thisGridProc%d, ToGridProc%d, this %p, ToGrid %p, thisDims %d %d %d, ToGridDism %d %d %d, thisLE %.2f %.2f %.2f, ToGridLE %.2f %.2f %.2f, thisRE %.2f %.2f %.2f, ToGridRE %.2f %.2f %.2f \n",
    // MyProcessorNumber, ProcessorNumber, ToProcessor, this, ToGrid, 
    // this->GridDimension[0], this->GridDimension[1], this->GridDimension[2],
    // ToGrid->GridDimension[0], ToGrid->GridDimension[1], ToGrid->GridDimension[2],
    // this->GridLeftEdge[0], this->GridLeftEdge[1], this->GridLeftEdge[2],
    // ToGrid->GridLeftEdge[0], ToGrid->GridLeftEdge[1], ToGrid->GridLeftEdge[2],
    // this->GridRightEdge[0], this->GridRightEdge[1], this->GridRightEdge[2],
    // ToGrid->GridRightEdge[0], ToGrid->GridRightEdge[1], ToGrid->GridRightEdge[2]); // DEBUG

    for (n = 0; n < TransferSize; n++) {
      mctp = MonteCarloTracerParticleBufferToList(buffer[n]);
      mctp->CurrentGrid = ToGrid;
      
      i = int((mctp->Position[0] - ToGrid->GridLeftEdge[0]) / ToGrid->CellWidth[0][0]);
      j = int((mctp->Position[1] - ToGrid->GridLeftEdge[1]) / ToGrid->CellWidth[1][0]);
      k = int((mctp->Position[2] - ToGrid->GridLeftEdge[2]) / ToGrid->CellWidth[2][0]);
      
      index = ((k + ToGrid->GridStartIndex[2]) * ToGrid->GridDimension[1] +
               (j + ToGrid->GridStartIndex[1])) * ToGrid->GridDimension[0] +
               (i + ToGrid->GridStartIndex[0]);
      //index = GRIDINDEX(i,j,k);

      // if (COM_COMBINE) // DEBUG
      //   printf("\nCOMSEND_INDEX, proc%d, thisProc%d, ToProc%d: (%d, %d, %d), pos (%.4f, %.4f, %.4f)", MyProcessorNumber, ProcessorNumber, ToGrid->ProcessorNumber, i, j, k, mctp->Position[0], mctp->Position[1], mctp->Position[2]);

      InsertMonteCarloTracerParticleAfter(ToGrid->MonteCarloTracerParticles[index], mctp);        
    }

    //ToGrid->WriteMCTP("ComSend_ToGrid_A1"); //DEBUG      


    printf("\nproc%d: CommSendMCTP: Inserted %d MCTPs from buffer into ToGrid. ToProcessor %d", MyProcessorNumber, TransferSize, ToProcessor);
    fflush(stdout);
    delete [] buffer;
			  
  } // end: if (CommunicationDirection...)

  return SUCCESS;
}

