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

int grid::Transfer_MonteCarloTracerParticles_From_CellA_to_CellB(MonteCarloTracerParticle *&headA, MonteCarloTracerParticle *&headB, double probability)
{
  // No particles to exchange
  if (headA == NULL)
      return 0;
  
  // Pointers for tracking the position in the particle list
  MonteCarloTracerParticle *previous, *next = NULL;
  MonteCarloTracerParticle *current = headA;

  double random_unit_interval;
  
  while (current != NULL)
  {
    // Store the next particle in cellA's particle list becuase the
    // current particle may get cut out of cellA's particle list.
    next = current->NextParticle;

    random_unit_interval = rand()/RAND_MAX;

    // Probabalistically transfer particle from cellA to cellB if it 
    // was not transfered into cellA durring the current timestep.
    if (random_unit_interval < probability && !current->EchangedThisTimestep)
    {
      // Exchange particle from cellA to cellB
      // cut current particle out of cellA particle list and
      // reassign headA if current == headA
      if (previous == NULL)
        headA = current->NextParticle;
      else
        previous->NextParticle = current->NextParticle;

      // Insert cellA_current at cellB_head and redirect cellB_head
      if headB != NULL
          current->next = headB->next;
      headB = current;
      
      // Update this particle's exchange count and set exchanged flag
      current->ExchangeCount++;
      current->EchangedThisTimestep = true;
      current = next;      
    }
    else // we skipped the current particle
    {
      previous = current;
      current = next;
    }
  }
  return 0;
}