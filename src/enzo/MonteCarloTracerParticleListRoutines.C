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
  if (Node == NULL)
    Node = NewNode;
  else {
    if (NewNode == Node->NextMonteCarloTracerParticle) {
      printf("Node already in list?!\n");
      exit(1);
    }
    NewNode->PrevMonteCarloTracerParticle = Node;
    NewNode->NextMonteCarloTracerParticle = Node->NextMonteCarloTracerParticle;
    if (Node->NextMonteCarloTracerParticle != NULL)
      Node->NextMonteCarloTracerParticle->PrevMonteCarloTracerParticle = NewNode;
    Node->NextMonteCarloTracerParticle = NewNode;
  }
  return;
}

MonteCarloTracerParticle *PopMonteCarloTracerParticle(MonteCarloTracerParticle * &Node)
{
  MonteCarloTracerParticle *result = Node;
  if (Node->PrevMonteCarloTracerParticle != NULL)
    Node->PrevMonteCarloTracerParticle->NextMonteCarloTracerParticle = Node->NextMonteCarloTracerParticle;
  if (Node->NextMonteCarloTracerParticle != NULL)
    Node->NextMonteCarloTracerParticle->PrevMonteCarloTracerParticle = Node->PrevMonteCarloTracerParticle;
  Node = Node->NextMonteCarloTracerParticle;
  result->NextMonteCarloTracerParticle = NULL;
  result->PrevMonteCarloTracerParticle = NULL;
  return result;
}

void DeleteMonteCarloTracerParticle(MonteCarloTracerParticle * &Node)
{
  MonteCarloTracerParticle *Orphan = PopMonteCarloTracerParticle(Node);
  //Node = Node->NextMonteCarloTracerParticle;
  if (Orphan != NULL) delete Orphan;
  return;
}

void DeleteMonteCarloTracerParticleList(MonteCarloTracerParticle * &Node)
{
  MonteCarloTracerParticle *tmp = Node;
  while (tmp)  // delete all linked stars
    DeleteMonteCarloTracerParticle(tmp);
  Node = NULL;
  return;
}

MonteCarloTracerParticle *MonteCarloTracerParticleListToArray(MonteCarloTracerParticle *Node, int n)
{
  int dim, count = 0;
  MonteCarloTracerParticle *result = new MonteCarloTracerParticle[n];
  MonteCarloTracerParticle *tmp = Node;
  while (tmp != NULL) {
    result[count++] = *tmp;
    tmp = tmp->NextMonteCarloTracerParticle;
  }
  return result;
}

/* Since InsertMonteCarloTracerParticleAfter puts the node after the head node.  We insert
   the nodes in a fashion to preserve the order of the array. */

MonteCarloTracerParticle* MonteCarloTracerParticleBufferToList(MonteCarloTracerParticleBuffer buffer)
{
  MonteCarloTracerParticle *result = NULL;
  result = new MonteCarloTracerParticle(buffer);
  return result;
}

MonteCarloTracerParticle* MonteCarloTracerParticleBufferToList(MonteCarloTracerParticleBuffer *buffer, int n) 
{
  int i;
  MonteCarloTracerParticle *result = NULL, *NewNode = NULL;
  if (n > 0) {
    NewNode = new MonteCarloTracerParticle(buffer, 0);
    InsertMonteCarloTracerParticleAfter(result, NewNode);
  }
  for (i = n-1; i > 0; i--) {
    NewNode = new MonteCarloTracerParticle(buffer, i);
    InsertMonteCarloTracerParticleAfter(result, NewNode);
  }
  return result;
}
