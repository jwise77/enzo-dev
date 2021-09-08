/***********************************************************************
/
/  GRID CLASS (Print Monte Carlo Tracer Particles)
/
/  written by: Corey Brummel-Smith
/  date:       September, 2021
/
/  PURPOSE: Print Monte Carlo Tracer Particle data to a file using
/           python dictionary syntax so that the output can be read
/           in by running it as python code.
/           
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ostream>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h" 

void grid::PrintMonteCarloTracerParticlePythonDictionary(std::ostream& os)
{
  FLOAT pos[3];
  int i, j, k, index, counter;
  MonteCarloTracerParticle *mc;
  MonteCarloTracerParticle *mc0;
  
  os << "\nMCTracers = {}" << std::endl;
  
  for (k = 0; k < GridDimension[2]; k++) {
    pos[2] = (k + 0.5) * CellWidth[2][0];
    for (j = 0; j < GridDimension[1]; j++) {
      pos[1] = (j + 0.5) * CellWidth[1][0];
      for (i = 0; i < GridDimension[0]; i++) {
        pos[0] = (i + 0.5) * CellWidth[0][0];
        counter = 0;
        os << "\nMCTracers[" << index << "] = {}";
  
        index = i + GridDimension[0]*(j + GridDimension[1]*k);
        mc = MonteCarloTracerParticles[index];
        mc0 = MonteCarloTracerParticles[index];
  
        if (mc0 != NULL)
          if (mc0->PrevParticle != NULL)
              os << "\nMCTracers[" << index << "] is NOT a MC Tracer HEAD";
  
        while (mc != NULL){
          
          os << "\nMCTracers[" << index << "]" << "[" << counter << "][initial_position] = (" << mc->InitialPosition[1] << ", " << mc->InitialPosition[2] << ", " << mc->InitialPosition[3] << ")"; 
          os << "\nMCTracers[" << index << "]" << "[" << counter << "][position] = (" << pos[1] << ", " << pos[2] << ", " << pos[3] << ")"; 
          os << "\nMCTracers[" << index << "]" << "[" << counter << "][UniqueID] = " << mc->UniqueID;
          //os << "\nMCTracers[" << index << "]" << "[" << counter << "][GroupID] = " << mc->GroupID;
          os << "\nMCTracers[" << index << "]" << "[" << counter << "][ExchangeCount] = " << mc->ExchangeCount;
          counter++;
          mc = mc->NextParticle;
        }
        os << "\nMCTracers[" << index << "]" << "[" <<  "N"  << "] = " << counter << std::endl;
      }
    }
  }
  return;
}
