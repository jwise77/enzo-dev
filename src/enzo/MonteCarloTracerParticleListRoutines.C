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

void InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticle * &Node, MonteCarloTracerParticle * &NewNode)
{
  //printf("\nInserting mc");
  if (Node == NULL)
    Node = NewNode;
  else {
    if (NewNode == Node->NextParticle) {
      //printf("\nNode already in list?!\n");
      exit(1);
    }
    //printf("\nNewNode->PrevParticle = Node; %p -> %p = %p", NewNode, NewNode->PrevParticle, Node);
    NewNode->PrevParticle = Node;
    //printf("\nNewNode->NextParticle = Node->NextParticle;; %p -> %p = %p -> ", NewNode, NewNode->NextParticle, Node, Node->NextParticle);
    NewNode->NextParticle = Node->NextParticle;
    if (Node->NextParticle != NULL)
      Node->NextParticle->PrevParticle = NewNode;
    Node->NextParticle = NewNode;
  }
  return;
}

MonteCarloTracerParticle *PopMonteCarloTracerParticle(MonteCarloTracerParticle * &Node)
{
  MonteCarloTracerParticle *result = Node;
  if (Node->PrevParticle != NULL)
    Node->PrevParticle->NextParticle = Node->NextParticle;
  if (Node->NextParticle != NULL)
    Node->NextParticle->PrevParticle = Node->PrevParticle;
  Node = Node->NextParticle;
  result->NextParticle = NULL;
  result->PrevParticle = NULL;
  return result;
}

void DeleteMonteCarloTracerParticle(MonteCarloTracerParticle * &Node)
{
  MonteCarloTracerParticle *Orphan = PopMonteCarloTracerParticle(Node);
  //printf("\nOrphan = %p", Orphan);
  if (Orphan != NULL) delete Orphan;
  return;
}

void DeleteMonteCarloTracerParticleList(MonteCarloTracerParticle * &Node)
{
  MonteCarloTracerParticle *tmp = Node;
  while (tmp)  // delete all linked MC tracers
    DeleteMonteCarloTracerParticle(tmp);
  Node = NULL;
  return;
}

// MonteCarloTracerParticle *MonteCarloTracerParticleListToArray(MonteCarloTracerParticle *Node, int n)
// {
//   int dim, count = 0;
//   MonteCarloTracerParticle *result = new MonteCarloTracerParticle[n];
//   MonteCarloTracerParticle *tmp = Node;
//   while (tmp != NULL) {
//     result[count++] = *tmp;
//     tmp = tmp->NextParticle;
//   }
//   return result;
// }

// /* Since InsertMonteCarloTracerParticleAfter puts the node after the head node.  We insert
//    the nodes in a fashion to preserve the order of the array. */

MonteCarloTracerParticle* MonteCarloTracerParticleBufferToList(MonteCarloTracerParticleBuffer buffer)
{
  MonteCarloTracerParticle *result = NULL;
  result = new MonteCarloTracerParticle(buffer);
  return result;
}

// MonteCarloTracerParticle* MonteCarloTracerParticleBufferToList(MonteCarloTracerParticleBuffer *buffer, int n) 
// {
//   int i;
//   MonteCarloTracerParticle *result = NULL, *NewNode = NULL;
//   if (n > 0) {
//     NewNode = new MonteCarloTracerParticle(buffer, 0);
//     InsertMonteCarloTracerParticleAfter(result, NewNode);
//   }
//   for (i = n-1; i > 0; i--) {
//     NewNode = new MonteCarloTracerParticle(buffer, i);
//     InsertMonteCarloTracerParticleAfter(result, NewNode);
//   }
//   return result;
// }


int Move_MonteCarloTracerParticles_From_CellA_to_CellB(MonteCarloTracerParticle *&headA, MonteCarloTracerParticle *&headB)
{
  // No particles to exchange
  if (headA == NULL)
      return 0;

  // Ensure headA and headB are actually heads 
  bool B_not_head = 0;
  if (headB != NULL)
    if (headB->PrevParticle != NULL)
      B_not_head = 1;
  if (headA->PrevParticle != NULL || B_not_head) {
    printf("%s\n", "headA and headB must MCTracer linked list heads");
    return -1;
  }
  
  // Pointers for tracking the position in the particle list
  MonteCarloTracerParticle *next = NULL;
  MonteCarloTracerParticle *transfer = NULL;
  MonteCarloTracerParticle *current = headA;
  
  while (current != NULL)
  {
    next = current->NextParticle;

    // Exchange particle from cellA to cellB  

    // Reassign this cell's head (headA) before transfering the particle to cellB
    //  if the current particle is at the head of this cell
    if (current->PrevParticle == NULL)
      headA = current->NextParticle;

    transfer = PopMonteCarloTracerParticle(current);
    InsertMonteCarloTracerParticleAfter(headB, transfer);

    current = next;  
  }
  return SUCCESS;
}
