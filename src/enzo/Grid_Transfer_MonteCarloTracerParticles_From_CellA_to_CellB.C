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


int grid::Transfer_MonteCarloTracerParticles_From_CellA_to_CellB(MonteCarloTracerParticle *&headA, MonteCarloTracerParticle *&headB, double probability)
{
  // No particles to exchange
  if (headA == NULL)
      return 0;

  // Ensure headA and headB are actually heads

  int a,b,c,d,e,f;
  
  // Pointers for tracking the position in the particle list
  MonteCarloTracerParticle *next = NULL;
  MonteCarloTracerParticle *transfer = NULL;
  MonteCarloTracerParticle *current = headA;

  double random_unit_interval;
  
  while (current != NULL)
  {
    // Store the next particle in cellA's particle list becuase the
    // current particle may get cut out of cellA's particle list.
    next = current->NextParticle; // unneccesary b/c current is updated in PopMCT

    random_unit_interval = (float) rand() / (float) RAND_MAX;
    //random_unit_interval = 0;

    // Probabalistically transfer particle from cellA to cellB if it 
    // was not transfered into cellA durring the current timestep.
    if (random_unit_interval < probability && !current->ExchangedThisTimestep)
    {
      // Exchange particle from cellA to cellB  

      // Reassign this cell's head (headA) before transfering the particle to cellB
      //  if the current particle is at the head of this cell
      if (current->PrevParticle == NULL)
        headA = current->NextParticle;

      transfer = PopMonteCarloTracerParticle(current);
      InsertMonteCarloTracerParticleAfter(headB, transfer);
      
      // Update this particle's exchange count and set exchanged flag
      transfer->ExchangeCount++;
      transfer->ExchangedThisTimestep = true;
      current = next;      
    }
    else // we skipped the current particle
    {
      current = next;
    }
  }
  return SUCCESS;
}