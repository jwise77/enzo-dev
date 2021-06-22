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

int grid::AdvectMonteCarloTracerParticles()
{

  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber || !UseHydro || !UseMonteCarloTracerParticles)
    return SUCCESS;	

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  int nxz, nyz, nzz, ixyz;

  nxz = GridEndIndex[0] - GridStartIndex[0] + 1;
  nyz = GridEndIndex[1] - GridStartIndex[1] + 1;
  nzz = GridEndIndex[2] - GridStartIndex[2] + 1;

  ixyz = CycleNumber % GridRank;

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
	
  // only need DensNum
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum);

  MonteCarloTracerParticle *headA, *headB; 
  float probability;
  int i, j, k, n, index2, index3;

  // Initialize reduced mass with the old mass (density) field
  float *reduced_mass = OldBaryonField[DensNum];

  for (n = ixyz; n < ixyz+GridRank; n++) {

    // Update in x-direction
    if ((n % GridRank == 0) && nxz > 1) {

	  for (j = 0; j < GridDimension[1]; j++) {
	    index2 = j * GridDimension[0];
	    for (i = 0; i < GridDimension[0]; i++) {
	      index3 = (k*GridDimension[1] + j) * GridDimension[0] + i;
	      probability =  // *** TODO ***
	    }
	  }
    } // ENDIF x-direction


    // Update in y-direction
    if ((n % GridRank == 1) && nyz > 1) {

	  for (k = 0; k < GridDimension[2]; k++) {
	    index2 = k * GridDimension[1];
	    for (j = 0; j < GridDimension[1]; j++) {
	      index3 = (k*GridDimension[1] + j) * GridDimension[0] + i;
	      probability =  // *** TODO ***
	    }
	  }   	
    } // ENDIF y-direction

      
      // Update in z-direction
    if ((n % GridRank == 2) && nzz > 1) {

	  for (i = 0; i < GridDimension[0]; i++) {
	    index2 = i * GridDimension[2];
	    for (k = 0; k < GridDimension[2]; k++) {
	      index3 = (k*GridDimension[1] + j) * GridDimension[0] + i;
	      probability = // *** TODO ***
	    }
	  }    	
	} // ENDIF y-direction	


}