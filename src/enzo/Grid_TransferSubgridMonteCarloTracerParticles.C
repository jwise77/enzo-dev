/***********************************************************************
/
/  GRID CLASS (COPY SUBGRID MONTE CARLO TRACER PARTICLES INTO OR OUT OF GRID)
/
/  written by: Greg Bryan
/  date:       January, 1999
/  modified1:  Robert Harkness
/  date:       April, 2006
/  modified2:  May, 2009 by John Wise: modified to move subgrid particles
/  modified2:  July, 2009 by John Wise: modified to move Monte Carlo tracer particles
/  modified3:  July, 2024 by Corey Brummel-Smith: modified to move Monte
/               Carlo tracer particles
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
#include "Hierarchy.h"

MonteCarloTracerParticle *PopMonteCarloTracerParticle(MonteCarloTracerParticle * &Node);
MonteCarloTracerParticle* MonteCarloTracerParticleBufferToList(MonteCarloTracerParticleBuffer buffer);
void InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticle * &Node, MonteCarloTracerParticle * &NewNode);
 
int grid::TransferSubgridMonteCarloTracerParticles(grid* Subgrids[], int NumberOfSubgrids, 
                   int* &NumberToMove, int iStart, 
                   int iEnd, mc_tracer_data* &List, 
                   bool KeepLocal, bool ParticlesAreLocal,
                   int CopyDirection, int IncludeGhostZones,
                               int CountOnly)
{
 
  /* Declarations. */

  int i, j, index, dim, n1, grid, proc;
  int i0, j0, k0;
  FLOAT pos[3];
  MonteCarloTracerParticle *mctp, *MoveMCTP;
  
  //IncludeGhostZones = 1;//****DEBUG****

  /* ----------------------------------------------------------------- */
  /* Copy Monte Carlo tracer particles out of grid. */

  if (CopyDirection == COPY_OUT) {

    /* If Monte Carlo tracer particles aren't distributed over several processors, exit
       if this isn't the host processor. */

    if (ParticlesAreLocal && MyProcessorNumber != ProcessorNumber)
      return SUCCESS;

    /* If there are no Monte Carlo tracer particles to move, we're done. */

    if (NumberOfMonteCarloTracerParticles == 0)
      return SUCCESS;

    /* Set boundaries (with and without ghost zones) */

    int StartIndex[] = {1,1,1}, EndIndex[] = {1,1,1};
    if (IncludeGhostZones)
      for (dim = 0; dim < GridRank; dim++) {
        StartIndex[dim] = 0;
        EndIndex[dim] = GridDimension[dim]-1;
      }
    else
      for (dim = 0; dim < GridRank; dim++) {
        StartIndex[dim] = GridStartIndex[dim];
        EndIndex[dim] = GridEndIndex[dim];
      }
 
    /* Count the number of Monte Carlo tracer particles already moved */

    int PreviousTotalToMove = 0;
    for (i = 0; i < NumberOfProcessors; i++)
      PreviousTotalToMove += NumberToMove[i];
 
    /* Count Monte Carlo tracer particles to move */

    
    int *subgrid = NULL;
    subgrid = new int[NumberOfMonteCarloTracerParticles];

    /* Loop over all cells and all mctps in cell. 
       Get subgrid number and assign subgrid processor*/


    i = 0;
    for (k0 = StartIndex[2]; k0 <= EndIndex[2]; k0++) {
      for (j0 = StartIndex[1]; j0 <= EndIndex[1]; j0++) {
        for (i0 = StartIndex[0]; i0 <= EndIndex[0]; i0++) {
 
          index = (k0*GridDimension[1] + j0)*GridDimension[0] + i0;

          // Compute particle position (cell-center)
          pos[2] = CellLeftEdge[2][k0] + 0.5 * CellWidth[2][0];
          pos[1] = CellLeftEdge[1][j0] + 0.5 * CellWidth[1][0];
          pos[0] = CellLeftEdge[0][i0] + 0.5 * CellWidth[0][0];

          /* Loop over particles in this cell */
                        
          mctp = MonteCarloTracerParticles[index];

          while (mctp != NULL) {

            for (dim = 0; dim < MAX_DIMENSION; dim++)
              mctp->Position[dim] = pos[dim];

            /* Find and store subgrid number of this mctp, and add to
           count. */
       
            subgrid[i] = nint(BaryonField[NumberOfBaryonFields][index])-1;
            if (subgrid[i] >= 0) {
              if (KeepLocal)
                proc = MyProcessorNumber;
              else
                proc = Subgrids[subgrid[i]]->ReturnProcessorNumber();
              NumberToMove[proc]++;
            }
            if (subgrid[i] < -1 || subgrid[i] > NumberOfSubgrids-1) {
              ENZO_VFAIL("mctp subgrid (%"ISYM"/%"ISYM") out of range\n", 
              subgrid[i], NumberOfSubgrids)
            }

            i++; // Increment particle counter
          
          } // End while mctp
        } // End for i0
      } // End for j0
    } // End for k0
      

    if (CountOnly == TRUE) {
      delete [] subgrid;
      return SUCCESS;
    }

 
    int TotalToMove = 0;
    for (proc = 0; proc < NumberOfProcessors; proc++)
      TotalToMove += NumberToMove[proc];
 
    if (TotalToMove > PreviousTotalToMove) {

      /* Move Monte Carlo tracer particles */

      MoveMonteCarloTracerParticlesToCellZero();

      n1 = PreviousTotalToMove;
      NumberOfMonteCarloTracerParticles = 0;
      mctp = MonteCarloTracerParticles[0];
      MonteCarloTracerParticles = NULL;
      i = 0;
      while (mctp != NULL) {

        MoveMCTP = PopMonteCarloTracerParticle(mctp); // also advances to NextParticle

        if (subgrid[i] >= 0) {
              MoveMCTP->MonteCarloTracerParticleToBuffer(&List[n1].data, MoveMCTP->Position);
          delete MoveMCTP;

          List[n1].grid = subgrid[i];
          List[n1].proc = (KeepLocal) ? MyProcessorNumber :
            Subgrids[subgrid[i]]->ReturnProcessorNumber();
          n1++;

        } // ENDIF move to subgrid

        else {
          InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticles[0], MoveMCTP);
          NumberOfMonteCarloTracerParticles++;
        }

        i++;

      } // ENDWHILE Monte Carlo tracer particles

    } // ENDIF Monte Carlo tracer particles to move

    DistributeMonteCarloTracerParticles();

    delete [] subgrid;
 
  } // end: if (COPY_OUT)
 
  /* ----------------------------------------------------------------- */
  /* Copy Monte Carlo tracer particles back into grid. */
 
  else {

    /* Count up total number. */
 
    int TotalNumberOfMonteCarloTracerParticles;
    int NumberOfNewMonteCarloTracerParticles = iEnd - iStart;

    TotalNumberOfMonteCarloTracerParticles = NumberOfMonteCarloTracerParticles + NumberOfNewMonteCarloTracerParticles;
 
    /* Copy Monte Carlo tracer particles from buffer into linked list */
    
    if (NumberOfNewMonteCarloTracerParticles > 0)
      for (i = iStart; i < iEnd; i++) {
        MoveMCTP = MonteCarloTracerParticleBufferToList(List[i].data);
        MoveMCTP->GridID = List[i].grid;
        MoveMCTP->CurrentGrid = this;
        // FALSE if going to a subgrid
        if (IncludeGhostZones == FALSE)
          MoveMCTP->IncreaseLevel();

        // Insert into cell zero. Will be distributed to correct cell later.
        InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticles[0], MoveMCTP);
      } // ENDFOR Monte Carlo tracer particles
 
    /* Set new number of Monte Carlo tracer particles in this grid. */
 
    NumberOfMonteCarloTracerParticles = TotalNumberOfMonteCarloTracerParticles;
  
  } // end: if (COPY_IN)

  DistributeMonteCarloTracerParticles();

 
  return SUCCESS;
}
