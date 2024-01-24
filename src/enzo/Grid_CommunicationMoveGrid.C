/***********************************************************************
/
/  GRID CLASS (MOVE A GRID FROM ONE PROCESSOR TO ANOTHER)
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
/  INPUTS:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif

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
 
/* function prototypes */
 
 
 
int grid::CommunicationMoveGrid(int ToProcessor, int MoveParticles, 
				int DeleteAllFields, int MoveSubgridMarker)
{
  printf("\nStart CommMoveGrid");

  int dim;
  int Zero[] = {0, 0, 0};
  FLOAT FZero[] = {0.0, 0.0, 0.0};

  //CommunicationDirection = COMMUNICATION_SEND_RECEIVE;
    //printf("\nALL CommMoveGrid: MyProcessorNumber %d, ProcessorNumber %d, ToProcessor %d\n", MyProcessorNumber, ProcessorNumber, ToProcessor);
  
  if ((MyProcessorNumber == ProcessorNumber ||
       MyProcessorNumber == ToProcessor) &&
      ProcessorNumber != ToProcessor) {
    printf("\nSOME CommMoveGrid: MyProcessorNumber %d, ProcessorNumber %d, ToProcessor %d\n", MyProcessorNumber, ProcessorNumber, ToProcessor);
    fflush(stdout);
    /* Copy baryons. */
 
    if (NumberOfBaryonFields > 0) {
#ifdef USE_MPI
      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE) {
	CommunicationReceiveGridOne[CommunicationReceiveIndex] = this;
	CommunicationReceiveGridTwo[CommunicationReceiveIndex] = this;
	CommunicationReceiveCallType[CommunicationReceiveIndex] = 16;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  CommunicationReceiveArgumentInt[dim][CommunicationReceiveIndex] =
	    GridDimension[dim];
      }
#endif
      this->CommunicationSendRegion(this, ToProcessor, ALL_FIELDS,
				    NEW_ONLY, Zero, GridDimension);
    }
 
    /* Copy particles. */

    if (NumberOfParticles > 0 && MoveParticles == TRUE)
      this->CommunicationSendParticles(this, ToProcessor, 0,
				       NumberOfParticles, 0);
    /* Copy active particles */

    if (MoveParticles == TRUE)
      this->CommunicationSendActiveParticles(this, ToProcessor);
    
    /* Copy stars */

    if (NumberOfStars > 0 && MoveParticles == TRUE)
      this->CommunicationSendStars(this, ToProcessor);

    /* Copy Monte Carlo tracer particles */
    //** DEBUG **
    char ToProcessorStr[1];
    char *filename = new char[MAX_LINE_LENGTH];
    sprintf(ToProcessorStr, "%d", ToProcessor);  
    strcpy(filename, "Pre_ComSend_in_ComMov");
    strcat(filename, ToProcessorStr);
    //** END DEBUG **
    this->WriteMCTP(filename);
    if (MoveParticles == TRUE){
      printf("\nTHREE proc%d, ToProc%d: ComMovGrid: this->NMCTP %d", MyProcessorNumber, ToProcessor, this->GetNumberOfMonteCarloTracerParticles());
      this->CommunicationSendMonteCarloTracerParticles(this, ToProcessor);   
    }
    strcpy(filename, "Post_ComSend_in_ComMov");
    strcat(filename, ToProcessorStr);       
    this->WriteMCTP(filename);
    delete filename;

    /* Copy photon packages */

#ifdef TRANSFER
    if (NumberOfPhotonPackages > 0)
      this->CommunicationSendPhotonPackages(this, ToProcessor, 
					    NumberOfPhotonPackages, 
					    NumberOfPhotonPackages, 
					    &PhotonPackages);
    if (MoveSubgridMarker == TRUE)
      this->CommunicationSendSubgridMarker(this, ToProcessor);
#endif /* TRANSFER */    

    /* Delete fields on old grid. */
 
    if (MyProcessorNumber == ProcessorNumber && ProcessorNumber != ToProcessor &&
	(CommunicationDirection == COMMUNICATION_SEND ||
	 CommunicationDirection == COMMUNICATION_SEND_RECEIVE)) {
      if (DeleteAllFields == TRUE) {
	if (MoveParticles == TRUE)
	  this->DeleteAllFields();
	else
	  this->DeleteAllButParticles();
#ifdef UNUSED
	if (MoveSubgridMarker == TRUE)
	  delete [] SubgridMarker;
#endif
      } else {
	this->DeleteBaryonFields();
      }
    }
    
  } // ENDIF right processor
 
  /* Update processor number. */
  
  if (CommunicationDirection == COMMUNICATION_SEND_RECEIVE)
    ProcessorNumber = ToProcessor;
 
  return SUCCESS;
}
 
