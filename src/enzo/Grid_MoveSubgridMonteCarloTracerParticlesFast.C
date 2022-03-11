/***********************************************************************
/
/  GRID CLASS (MOVE MONTE CARLO TRACER PARTICLES FROM SPECIFIED GRID TO
/              THIS GRID)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Corey Brummel-Smith (modified MoveSubgridParticlesFast to
/               work with Monte Carlo tracer particles)
/  date:       Jan. 2022 
/
/  PURPOSE:
/
************************************************************************/
 
//
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"

 
int CommunicationBroadcastValue(int *Value, int BroadcastProcessor);
void InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticle * &Node, MonteCarloTracerParticle * &NewNode);
MonteCarloTracerParticle *PopMonteCarloTracerParticle(MonteCarloTracerParticle * &Node); 
 
int grid::MoveSubgridMonteCarloTracerParticlesFast(int NumberOfSubgrids, grid* ToGrids[],
				   int AllLocal)
{

  /* If there are no subgrids or particles to move, we're done. */
  if (NumberOfSubgrids == 0)
    return SUCCESS;
 
  int i, j, k, i0, j0, k0, dim, index, subgrid, n;
  FLOAT pos[3];
  MonteCarloTracerParticle *mctp, *MoveMCTP;
 
  /* Initialize. */
 
  int *ParticlesToMove = new int[NumberOfSubgrids];
  mc_tracer_data *SendList = NULL;

  for (i = 0; i < NumberOfSubgrids; i++)
    ParticlesToMove[i] = 0;
 
  /* Error check. */
 
  if (BaryonField[NumberOfBaryonFields] == NULL &&
      MyProcessorNumber == ProcessorNumber) {
    ENZO_FAIL("Subgrid field not present.\n");
  }


  /* Loop over particles and count the number in each subgrid. */
  
  if (MyProcessorNumber == ProcessorNumber) {

    int NumberOfMonteCarloTracerParticles = this->CountMonteCarloTracerParticles();
    //printf("\nproc%d: MoveSubgridMonteCarloTracerParticlesFast: Number of MCTPs = %"ISYM"", MyProcessorNumber, NumberOfMonteCarloTracerParticles);

    if (NumberOfMonteCarloTracerParticles == 0)
      return SUCCESS; 
         
    for (k0 = GridStartIndex[2]; k0 <= GridEndIndex[2]; k0++) {
      for (j0 = GridStartIndex[1]; j0 <= GridEndIndex[1]; j0++) {
        for (i0 = GridStartIndex[0]; i0 <= GridEndIndex[0]; i0++) {
 
          /* Compute the cell index for particles in this cell */
     
          index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;   

          /* Loop over all particles in this cell */
          for (mctp = this->MonteCarloTracerParticles[index]; mctp != NULL; 
               mctp = mctp->NextParticle) { 

            /* Find subgrid number of this particle, and add to count. */
       
            subgrid = nint(BaryonField[NumberOfBaryonFields][index])-1;

            if (subgrid >= 0)
              ParticlesToMove[subgrid]++;
            if (subgrid < -1 || subgrid > NumberOfSubgrids-1) {
              ENZO_VFAIL("particle subgrid (%"ISYM"/%"ISYM") out of range\n", subgrid,
               NumberOfSubgrids)
            }            
          } // end: loop over particles 
        } // end: loop over i
      } // end: loop over j
    } // end: loop over k
  } // end: if (MyProcessorNumber)
 
  /* Communicate number of send particles to subgrids */
  if (AllLocal == FALSE)
    for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
      if (CommunicationBroadcastValue(&ParticlesToMove[subgrid],
				      ProcessorNumber) == FAIL) {
	       ENZO_FAIL("Error in CommunicationBroadcastValue.\n");
      }

  /* Allocate space on all the subgrids with particles.
     *** THIS SHOULD PROBABLY BE DONE REGARDLESS OF WHETHER THERE ARE PARTICLES THERE OR NOT? *** */
 
  if (MyProcessorNumber == ProcessorNumber)
    for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++)
 
      if (ParticlesToMove[subgrid] > 0) {
 
        /* Check if MonteCarloTracerParticles have already been allocated */
      	if (ToGrids[subgrid]->MonteCarloTracerParticles != NULL) {
      	  ENZO_VFAIL("\nproc%d: Monte Carlo tracer particles already in subgrid %"ISYM" (n=%"ISYM", nm=%"ISYM")\n",
      		  MyProcessorNumber, subgrid, ToGrids[subgrid]->MonteCarloTracerParticles, ParticlesToMove[subgrid])
      	}
        
      	ToGrids[subgrid]->AllocateMonteCarloTracerParticleData();
       
      	//printf("\nproc%d: MoveSubgridMonteCarloTracerParticles: subgrid[%"ISYM"] = %"ISYM"", MyProcessorNumber, subgrid, ParticlesToMove[subgrid]);
 
      } // end: if (ParticlesToMove > 0)

  //printf("\nproc%d:\nthis(OldGrid) %p\nToGrids[0] %p\nToGrids[1] %p\n", MyProcessorNumber, this, ToGrids[0], ToGrids[1]); 
 
  if (MyProcessorNumber == ProcessorNumber) {  
 
    /* Loop over particles and move them to the appropriate ToGrid, depending
       on their position. */
    int count = 0;
    for (k0 = GridStartIndex[2]; k0 <= GridEndIndex[2]; k0++) {
      for (j0 = GridStartIndex[1]; j0 <= GridEndIndex[1]; j0++) {
        for (i0 = GridStartIndex[0]; i0 <= GridEndIndex[0]; i0++) {
          count++;
 
          index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;

          // Compute particle position (cell-center)
          pos[2] = CellLeftEdge[2][k0] + 0.5 * CellWidth[2][0];
          pos[1] = CellLeftEdge[1][j0] + 0.5 * CellWidth[1][0];
          pos[0] = CellLeftEdge[0][i0] + 0.5 * CellWidth[0][0];

          /* Find subgrid number of this particle, and move it. */
     
          subgrid = nint(BaryonField[NumberOfBaryonFields][index])-1;
          

            /* Loop over particles in this cell */
                        
          mctp = MonteCarloTracerParticles[index];
          MonteCarloTracerParticles[index] = NULL;
          while (mctp != NULL) {
            if (subgrid >= 0) {

              /* Pop particle from this grid/cell and insert it into the subgrid (aka ToGrid) */

              MoveMCTP = PopMonteCarloTracerParticle(mctp);  // also advances mctp to NextParticle

              /* set particle position */
              for (dim = 0; dim < MAX_DIMENSION; dim++)
                MoveMCTP->Position[dim] = pos[dim];

              /* Insert this particle into the first cell of the subgrid. This is not where it belongs
                 but at this point, no communication has been done so these particles live on a "fake"
                 grid on this processor,  the "real" grid may be on a different processor. These 
                 particles are temporarily stored in cell 0 and then put into the correct cell on the
                 "real" grid when they are received in CommunicationSendParticles(). */
              InsertMonteCarloTracerParticleAfter(ToGrids[subgrid]->MonteCarloTracerParticles[0], MoveMCTP);
              this->NumberOfMonteCarloTracerParticles--;
              ToGrids[subgrid]->NumberOfMonteCarloTracerParticles++;
              //printf("\nproc%d: inserted particle %d (%p) from this grid (%p), into ToGrids[%d] (%p)", MyProcessorNumber, count, MoveMCTP, this, subgrid, ToGrids[subgrid]);
            } // end: if (subgrid >= 0)
            else
              InsertMonteCarloTracerParticleAfter(this->MonteCarloTracerParticles[index], MoveMCTP);
          } // end: while (mctp != NULL)          
        } // end: loop over i
      } // end: loop over j
    } // end: loop over k
    
    delete [] BaryonField[NumberOfBaryonFields];
    BaryonField[NumberOfBaryonFields] = NULL;
 
  } // end: if (MyProcessorNumber)
 
  /* Transfer particles from fake to real grids (and clean up).
     **** MAYBE JUST DELETE THIS. MIGHT BE NEEDED WHEN CALLED IN REBUILD HIERARCHY. 
     NEED TO CHECK THIS **** */
  
  for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++) {
    if ((MyProcessorNumber == ProcessorNumber ||
         MyProcessorNumber == ToGrids[subgrid]->ProcessorNumber) &&
	       ProcessorNumber != ToGrids[subgrid]->ProcessorNumber)
      if (ParticlesToMove[subgrid] != 0) {
        printf("\nthis->ProcessorNumber: %d\nToGrids->ProcessorNumber: %d\n", this->ProcessorNumber, ToGrids[subgrid]->ProcessorNumber);
      	if (this->CommunicationSendMonteCarloTracerParticles(ToGrids[subgrid], ToGrids[subgrid]->ProcessorNumber)
      	    == FAIL) {
      	  ENZO_FAIL("Error in grid->CommunicationSendParticles.\n");
      	}
      	if (MyProcessorNumber == ProcessorNumber)
      	  ToGrids[subgrid]->DeleteAllFields();
      }
  }
 
  delete [] ParticlesToMove;
 

  // DEBUG
  // for (k0 = 0; k0 <= GridDimension[2]; k0++) {
  //   for (j0 = 0; j0 <= GridDimension[1]; j0++) {
  //     for (i0 = 0; i0 <= GridDimension[0]; i0++) {
  //       index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;
  //       if ((k0 < GridStartIndex[2] or k0 > GridEndIndex[2]) or
  //           (j0 < GridStartIndex[1] or j0 > GridEndIndex[1]) or
  //           (i0 < GridStartIndex[0] or i0 > GridEndIndex[0]))
  //           printf("\n* Ghost Cell * ");
  //       else
  //           printf("\nActive Cell "); 
  //         mctp = this->MonteCarloTracerParticles[index];
  //         printf(" MCTP[%d](0) %p", index, mctp);  
  //         int count = 1;   
  //         while (mctp != NULL) {
  //           mctp = mctp->NextParticle;
  //           printf("\nMCTP[%d](%d) %p", index, count, mctp);  
  //           count++;   
  //         }
  //     }
  //   }
  // }
  // // END DEBUG

  return SUCCESS;
}
