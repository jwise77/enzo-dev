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

  MonteCarloTracerParticle *tracersCellC, *tracersCellL, *tracersCellR; 
  float DeltaML, DeltaMR, probability;    
  int nxz, nyz, nzz, ixyz, iLR;
  int i, j, k, n, nLR, indexC, indexL, indexR;
  int istart, iend, jstart, jend, kstart, kend;

  // Initialize reduced mass with the old mass (density) field
  float ReducedMass[size];
  for (indexC = 0; indexC < size; indexC++)
    ReducedMass[indexC] = OldBaryonField[DensNum][indexC];

  nxz = GridEndIndex[0] - GridStartIndex[0] + 1;
  nyz = GridEndIndex[1] - GridStartIndex[1] + 1;
  nzz = GridEndIndex[2] - GridStartIndex[2] + 1;  
  ixyz = CycleNumber % GridRank; // to permute which dimension we sweep through first
  iLR  = CycleNumber % 2;        // to permute which cell (left, right) we transfer particles to first

  for (n = ixyz; n < ixyz+GridRank; n++) {

    // Update in x-direction
    if ((n % GridRank == 0) && nxz > 1) {

      istart = GridStartIndex[0];
      iend   = GridEndIndex[0] + 1;
      jstart = 0;
      jend   = GridDimension[1];
      kstart = 0;
      kend   = GridDimension[2];      

      for (k = kstart; k < kend; k++) {
        for (j = jstart; j < jend; j++) {
          for (i = istart; i < iend; i++) {  
              
    	      indexC = (k * GridDimension[1] + j) * GridDimension[0] + i;
            indexL = (k * GridDimension[1] + j) * GridDimension[0] + i - 1;
            indexR = (k * GridDimension[1] + j) * GridDimension[0] + i + 1;

            DeltaML = MassFlux[0][indexC];  // mass transfered into   this cell from the left
            DeltaMR = MassFlux[0][indexR];  // mass transfered out of this cell from the right
            // dnu(i) = dslice(i,j) + (df(i,j) - df(i+1,j)) 
            //                        ^------^   ^--------^
            //                           ∆ML        ∆MR   

            for (nLR = iLR; nLR < iLR+2; nLR++) {
              if (nLR % 2 == 0) {              
                /* Exchange particles for outgoing left mass fluxes (DeltaML < 0) */
                if (DeltaML < 0) {
                  probability  =  DeltaML / ReducedMass[indexC];
                  tracersCellC = MonteCarloTracerParticles[indexC];
                  tracersCellL = MonteCarloTracerParticles[indexL];
                  this->Transfer_MonteCarloTracerParticles_From_CellA_to_CellB(tracersCellC, 
                    tracersCellL, probability);

                  // Update reduced mass
                  ReducedMass[indexC] += DeltaML;
                }   
              }              
              if (nLR % 2 == 1) {
                /* Exchange particles for outgoing right mass fluxes (DeltaMR > 0) */
                if (DeltaMR > 0) {
        	        probability  =  DeltaMR / ReducedMass[indexC];
                  tracersCellC = MonteCarloTracerParticles[indexC];
                  tracersCellR = MonteCarloTracerParticles[indexR];
                  this->Transfer_MonteCarloTracerParticles_From_CellA_to_CellB(tracersCellC, 
                    tracersCellR, probability);

                  // Update reduced mass
                  ReducedMass[indexC] -= DeltaMR;
                }
              }
            } // END loop over iLR (left and right cells)         
    	    } // END loop over i
    	  } // END loop over j
      } // END loop over k
    } // ENDIF x-direction


    // Update in y-direction
    if ((n % GridRank == 1) && nyz > 1) {

      jstart = GridStartIndex[1];
      jend   = GridEndIndex[1] + 1;
      kstart = 0;
      kend   = GridDimension[2];
      istart = 0;
      iend   = GridDimension[0];

      for (i = istart; i < iend; i++) {
        for (k = kstart; k < kend; k++) {
          for (j = jstart; j < jend; j++) {        
    	      indexC = (k * GridDimension[1] + j  ) * GridDimension[0] + i;
            indexL = (k * GridDimension[1] + j-1) * GridDimension[0] + i;
            indexR = (k * GridDimension[1] + j+1) * GridDimension[0] + i;

            DeltaML = MassFlux[1][indexC];  // mass transfered into   this cell from the left
            DeltaMR = MassFlux[1][indexR];  // mass transfered out of this cell from the right

            for (nLR = iLR; nLR < iLR+2; nLR++) {
              if (nLR % 2 == 0) {              
                /* Exchange particles for outgoing left mass fluxes (DeltaML < 0) */
                if (DeltaML < 0) {
                  probability  =  DeltaML / ReducedMass[indexC];
                  tracersCellC = MonteCarloTracerParticles[indexC];
                  tracersCellL = MonteCarloTracerParticles[indexL];
                  this->Transfer_MonteCarloTracerParticles_From_CellA_to_CellB(tracersCellC, 
                    tracersCellL, probability);

                  // Update reduced mass
                  ReducedMass[indexC] += DeltaML;
                }   
              }              
              if (nLR % 2 == 1) {
                /* Exchange particles for outgoing right mass fluxes (DeltaMR > 0) */
                if (DeltaMR > 0) {
                  probability  =  DeltaMR / ReducedMass[indexC];
                  tracersCellC = MonteCarloTracerParticles[indexC];
                  tracersCellR = MonteCarloTracerParticles[indexR];
                  this->Transfer_MonteCarloTracerParticles_From_CellA_to_CellB(tracersCellC, 
                    tracersCellR, probability);

                  // Update reduced mass
                  ReducedMass[indexC] -= DeltaMR;
                }
              }
            } // END loop over iLR (left and right cells)                   
    	    } // END j
    	  } // END k
      } // END i
    } // ENDIF y-direction

      
    // Update in z-direction
    if ((n % GridRank == 2) && nzz > 1) {

      kstart = GridStartIndex[2];
      kend   = GridEndIndex[2] + 1;
      istart = 0;
      iend   = GridDimension[0];
      jstart = 0;
      jend   = GridDimension[1];      

      for (j = 0; j < GridDimension[1]; j++) {
        for (i = istart; i < iend; i++) {
          for (k = kstart; k < kend; k++) {        

    	      indexC = (  k   * GridDimension[1] + j) * GridDimension[0] + i;
            indexL = ((k-1) * GridDimension[1] + j) * GridDimension[0] + i;
            indexR = ((k+1) * GridDimension[1] + j) * GridDimension[0] + i;

            DeltaML = MassFlux[2][indexC]; // mass transfered into   this cell from the left
            DeltaMR = MassFlux[2][indexR]; // mass transfered out of this cell from the right

            for (nLR = iLR; nLR < iLR+2; nLR++) {
              if (nLR % 2 == 0) {              
                /* Exchange particles for outgoing left mass fluxes (DeltaML < 0) */
                if (DeltaML < 0) {
                  probability  =  DeltaML / ReducedMass[indexC];
                  tracersCellC = MonteCarloTracerParticles[indexC];
                  tracersCellL = MonteCarloTracerParticles[indexL];
                  this->Transfer_MonteCarloTracerParticles_From_CellA_to_CellB(tracersCellC, 
                    tracersCellL, probability);

                  // Update reduced mass
                  ReducedMass[indexC] += DeltaML;
                }   
              }              
              if (nLR % 2 == 1) {
                /* Exchange particles for outgoing right mass fluxes (DeltaMR > 0) */
                if (DeltaMR > 0) {
                  probability  =  DeltaMR / ReducedMass[indexC];
                  tracersCellC = MonteCarloTracerParticles[indexC];
                  tracersCellR = MonteCarloTracerParticles[indexR];
                  this->Transfer_MonteCarloTracerParticles_From_CellA_to_CellB(tracersCellC, 
                    tracersCellR, probability);

                  // Update reduced mass
                  ReducedMass[indexC] -= DeltaMR;
                }
              }
            } // END loop over iLR (left and right cells)                  
    	    } // END k
    	  } // END i
      } // END j
	  } // ENDIF y-direction	
  } // END n 
  return SUCCESS;
}