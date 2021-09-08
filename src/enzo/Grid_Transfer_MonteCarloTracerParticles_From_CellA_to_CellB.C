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
  printf("\n%s", "Transfer_MonteCarloTracerParticles_From_CellA_to_CellB...");

  // No particles to exchange
  if (headA == NULL)
      return 0;

  // Ensure headA and headB are actually heads

  int a,b,c,d,e,f;

  // a = (headA->PrevParticle ? headA->PrevParticle->UniqueID : -1);
  // b = headB ? (headB->PrevParticle ? headB->PrevParticle->UniqueID : -1) : -2;
  // c = (headA->NextParticle ? headA->NextParticle->UniqueID : -1);
  // d = headB ? (headB->NextParticle ? headB->NextParticle->UniqueID : -1) : -2;

  // printf("\nheadA %p, headB %p\n", headA, headB);        
  // printf("headA_ID       = %i, headB_ID       = %i\n", headA->UniqueID, headB->UniqueID);
  // printf("headA->Prev_ID = %i, headB->Prev_ID = %i\n", a, b);
  // printf("headA->Next_ID = %i, headB->Next_ID = %i\n", c, d);    

  // if (headA->PrevParticle != NULL || headB->PrevParticle != NULL) {
  //   printf("%s\n", "headA and headB must MCTracer linked list heads");
  //   return -1;
  // }
  
  // Pointers for tracking the position in the particle list
  MonteCarloTracerParticle *next = NULL;
  MonteCarloTracerParticle *transfer = NULL;
  MonteCarloTracerParticle *current = headA;


  double random_unit_interval;
  
  while (current != NULL)
  {
    //printf("\n\t%s", " current != NULL");

    // Store the next particle in cellA's particle list becuase the
    // current particle may get cut out of cellA's particle list.
    next = current->NextParticle; // unneccesary b/c current is updated in PopMCT

    random_unit_interval = (float) rand() / (float) RAND_MAX;

    // Probabalistically transfer particle from cellA to cellB if it 
    // was not transfered into cellA durring the current timestep.
    if (random_unit_interval < probability && !current->EchangedThisTimestep)
    {
      // a = current->PrevParticle ? current->PrevParticle->UniqueID : -1;
      // b = headB->PrevParticle ? headB->PrevParticle->UniqueID : -1;
      // c = current->NextParticle ? current->NextParticle->UniqueID : -1;
      // d = headB->NextParticle ? headB->NextParticle->UniqueID : -1;

      // printf("\ncurrent_ID       = %i, headB_ID       = %i\n", current->UniqueID, headB->UniqueID);
      // printf("current->Prev_ID = %i, headB->Prev_ID = %i\n", a, b);
      // printf("current->Next_ID = %i, headB->Next_ID = %i\n", c, d);  
      // printf("\t%s", "Exchange...");

      // Exchange particle from cellA to cellB  

      // Reassign this cell's head (headA) before transfering the particle to cellB
      //  if the current particle is at the head of this cell
      if (current->PrevParticle == NULL)
        headA = current->NextParticle;

      transfer = PopMonteCarloTracerParticle(current);
      InsertMonteCarloTracerParticleAfter(headB, transfer);

      // if (headA){
      //   a = headA->PrevParticle ? headA->PrevParticle->UniqueID : -1;
      //   c = headA->NextParticle ? headA->NextParticle->UniqueID : -1;
      // }
      // else{
      //   a = b = -2;
      // }
      // b = headB->PrevParticle ? headB->PrevParticle->UniqueID : -1;
      // d = headB->NextParticle ? headB->NextParticle->UniqueID : -1;
      // e = headA ? headA->UniqueID : -2;
      // f = headB ? headB->UniqueID : -2;

      // printf("\nheadA_ID       = %i, headB_ID       = %i\n", e, f);
      // printf("headA->Prev_ID = %i, headB->Prev_ID = %i\n", a, b);
      // printf("headA->Next_ID = %i, headB->Next_ID = %i\n", c, d);
      // printf("headA %p, current %p, transfer %p, headB %p", headA, current, transfer, headB);        
      
      // Update this particle's exchange count and set exchanged flag
      transfer->ExchangeCount++;
      transfer->EchangedThisTimestep = true;
      current = next;      
    }
    else // we skipped the current particle
    {
      current = next;
    }
  }
  printf("\n%s", " Done.");
  return SUCCESS;
}