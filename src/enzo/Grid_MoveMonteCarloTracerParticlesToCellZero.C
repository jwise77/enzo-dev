/***********************************************************************
/
/  GRID CLASS (Move Monte Carlo Tracer Particles To Cell Zero)
/
/  written by: Corey Brummel-Smith
/  date:       June, 2021
/
/  PURPOSE: Move Monte Carlo Tracer particles from their
/           correct position (cell) to cell 0. This is necessary because 
/           CommunicationSendMonteCarloTracerParticles assumes the
/           particles have already been gathered in cell zero before
/           sending them to another processor.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h" 

void InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticle * &Node, MonteCarloTracerParticle * &NewNode);
MonteCarloTracerParticle *PopMonteCarloTracerParticle(MonteCarloTracerParticle * &Node); 

int grid::MoveMonteCarloTracerParticlesToCellZero()
{

  int i0, j0, k0, dim, index, subgrid, n;
  FLOAT pos[3];
  MonteCarloTracerParticle *mctp, *MoveMCTP;

  if (MyProcessorNumber == ProcessorNumber) {  

    if (this->MonteCarloTracerParticles == NULL)
      ENZO_FAIL("\nMoveMonteCarloTracerParticlesToCellZero: Attempting to Move Monte Carlo Tracer particles in a grid that has not been initialized with Monte Carlo Tracer particles.\n");
 
    /* Loop over all cells and move particles in each cell to cell 0 */
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

          /* Loop over particles in this cell */
                        
          mctp = MonteCarloTracerParticles[index];
          MonteCarloTracerParticles[index] = NULL;
          while (mctp != NULL) {

            /* Pop particle from this cell and insert it into cell 0 */

            MoveMCTP = PopMonteCarloTracerParticle(mctp);  // also advances mctp to NextParticle

            /* set particle position */
            for (dim = 0; dim < MAX_DIMENSION; dim++)
              MoveMCTP->Position[dim] = pos[dim];
            /* Insert this particle into the first cell of this grid. This is not where it belongs
               but it must be placed in cell zero before sending it to another processor. These 
               particles are put into the correct cell on the proper grid when they are received 
               in CommunicationSendParticles(). */
            InsertMonteCarloTracerParticleAfter(this->MonteCarloTracerParticles[0], MoveMCTP);
          } // end: while (mctp != NULL)          
        } // end: loop over i
      } // end: loop over j
    } // end: loop over k
  }
  return SUCCESS;
}