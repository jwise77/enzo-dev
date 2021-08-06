/***********************************************************************
/
/  GRID CLASS (Advect Monte Carlo Tracer Particles)
/
/  written by: Corey Brummel-Smith
/  date:       June, 2021
/
/  PURPOSE: Advect all Monte Carlo tracer particles on this grid. 
/           Particles are exchanged from one cell to its neighbors 
/           probabilistically based on the outgoing mass flux from 
/           one cell to another. Just like SolvePPM_DE, this is done
/           one dimension at a time. Each cycle, the order (e.g. y, z, x)
/           is cycled to avoid a preference for one direction.
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

int grid::AdvectMonteCarloTracerParticles(int CycleNumber)
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber || !UseHydro || !MonteCarloTracerParticlesOn)
    return SUCCESS; 

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
    
  // only need DensNum
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum);

  MonteCarloTracerParticle *tracersCellC, *tracersCellR; 
  float DeltaMR, probability;    
  int nxz, nyz, nzz, ixyz, iLR, face;
  int n, nLR, dim1, dim2, dim3, indexC, indexL, indexR;

  // Initialize reduced mass with the old mass (density) field
  float ReducedMass[size];
  for (indexC = 0; indexC < size; indexC++)
    ReducedMass[indexC] = OldBaryonField[DensNum][indexC];

  nxz = GridEndIndex[0] - GridStartIndex[0] + 1;
  nyz = GridEndIndex[1] - GridStartIndex[1] + 1;
  nzz = GridEndIndex[2] - GridStartIndex[2] + 1;  

  ixyz = CycleNumber % GridRank; // to permute which dimension we sweep through first
  iLR  = CycleNumber % 2;        // to permute which cell (left, right) we transfer particles to first

  /* Start and end indicies for each dimension of the face loop*/
  // outermost loop
  int dim3start[3] = {0, 0, 0};
  int dim3end[3]   = {GridDimension[2], GridDimension[0], GridDimension[1]};
  // mid loop
  int dim2start[3] = {0, 0, 0};
  int dim2end[3]   = {GridDimension[1], GridDimension[2], GridDimension[0]};
  // innermost loop
  int dim1start[3] = {GridStartIndex[0],     GridStartIndex[1],     GridStartIndex[2]};
  int dim1end[3]   = {GridEndIndex[0] + 1,   GridEndIndex[1] + 1,   GridEndIndex[2] + 1};  

  // Loop over faces. Start with either x, y, or z based on CycleNumber.
  for (n = ixyz; n < ixyz+GridRank; n++) {

    if ((n % GridRank == 0) && nxz > 1)
      face = 0 // x-sweep
    if ((n % GridRank == 1) && nyz > 1)
      face = 1 // y-sweep
    if ((n % GridRank == 2) && nzz > 1)
      face = 2 // z-sweep

    for (dim3 = dim3start[face]; dim3 < dim3end[face]; dim3++) {
      for (dim2 = dim2start[face]; dim2 < dim2end[face]; dim2++) {
        for (dim1 = dim1start[face]; dim1 < dim1end[face]; dim1++) {  
            
          if (face == 0) { 
            indexC = (dim3 * GridDimension[1] + dim2) * GridDimension[0] + dim1;
            indexL = (dim3 * GridDimension[1] + dim2) * GridDimension[0] + dim1-1;
            indexR = (dim3 * GridDimension[1] + dim2) * GridDimension[0] + dim1+1;
          }
          else if (face == 1) { 
            indexC = (dim2 * GridDimension[1] + dim1)   * GridDimension[0] + dim3;
            indexL = (dim2 * GridDimension[1] + dim1-1) * GridDimension[0] + dim3;
            indexR = (dim2 * GridDimension[1] + dim1+1) * GridDimension[0] + dim3;
          }
          else if (face == 2) { 
            indexC = ((dim1)   * GridDimension[1] + dim3) * GridDimension[0] + dim2;
            indexL = ((dim1-1) * GridDimension[1] + dim3) * GridDimension[0] + dim2;
            indexR = ((dim1+1) * GridDimension[1] + dim3) * GridDimension[0] + dim2;
          }
          else{
            printf("%s\n", "Something very bad and unexpected happened");
          }

          DeltaML = MassFlux[face][indexC]; // mass transfered  into  this cell from the left
          DeltaMR = MassFlux[face][indexR]; // mass transfered out of this cell  to  the right
          // dnu(i) = dslice(i,j) + (df(i,j) - df(i+1,j)) 
          //                        ^------^   ^--------^
          //                           ∆ML        ∆MR             

          for (nLR = iLR; nLR < iLR+2; nLR++) {
            if (nLR % 2 == 0) {              
              /* Exchange particles for outgoing left mass fluxe (DeltaML < 0) */
              if (DeltaML < 0) {
                probability = DeltaML / ReducedMass[indexC];
                this->Transfer_MonteCarloTracerParticles_From_CellA_to_CellB(
                  MonteCarloTracerParticles[indexC],
                  MonteCarloTracerParticles[indexL],
                  probability);

                // Update reduced mass
                ReducedMass[indexC] += DeltaML;
              }   
            }              
            if (nLR % 2 == 1) {
              /* Exchange particles for outgoing right mass fluxe (DeltaMR > 0) */
              if (DeltaMR > 0) {
                probability = DeltaMR / ReducedMass[indexC];
                this->Transfer_MonteCarloTracerParticles_From_CellA_to_CellB(
                  MonteCarloTracerParticles[indexC],
                  MonteCarloTracerParticles[indexR],
                  probability);

                // Update reduced mass
                ReducedMass[indexC] -= DeltaMR;
              }
            }
          } // END loop over iLR (left and right cells)  
        } // END loop over dim1
      } // END loop over dim2
    } // END loop over dim3
  } // END n 
}
