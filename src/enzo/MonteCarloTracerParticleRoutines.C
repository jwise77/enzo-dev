/***********************************************************************
/
/  ROUTINES FOR THE MonteCarloTracerParticle PARTICLE CLASS
/
/  written by: Corey Brummel-Smith
/  date:       July, 2021
/
/  PURPOSE: Constructors and destructors and helper functions for 
/           MonteCarloTracerParticles 
/
************************************************************************/
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "phys_constants.h"

// void DeleteMonteCarloTracerParticle(MonteCarloTracerParticle * &Node);
// MonteCarloTracerParticle *PopMonteCarloTracerParticle(MonteCarloTracerParticle * &Node);
// void InsertMonteCarloTracerParticleAfter(MonteCarloTracerParticle * &Node, MonteCarloTracerParticle * &NewNode);

// int GetUnits(float *DensityUnits, float *LengthUnits,
// 	     float *TemperatureUnits, float *TimeUnits,
// 	     float *VelocityUnits, FLOAT Time);

/*******************************

   CONSTRUCTORS AND DESTRUCTOR

 *******************************/

MonteCarloTracerParticle::MonteCarloTracerParticle(void)
{

  int i;

  CurrentGrid = NULL;
  NextParticle = NULL;
  PrevParticle = NULL;
  UniqueID = GroupID = Level = ExchangeCount = 0;
  CreationTime = Mass = 0.0;
  EchangedThisTimestep = false;

  // user defined attributes (e.g. max temperature, group id)
  ParticleAttributes = new float[NumberOfParticleAttributes];
  for (i = 0; i < NumberOfParticleAttributes; i++)
    ParticleAttributes[i] = 0.0;

  for (i = 0; i < MAX_DIMENSION; i++)
    InitialPosition[i] = 0.0;
  
  // history of particle position, including the time when the position was recorded. 
  #ifdef TRACK_MC_HISTORY
  LagrangianHistory = NULL; 
  #endif
}


MonteCarloTracerParticle::MonteCarloTracerParticle(grid *_grid, int _ID, int _groupID, int _level, float _creationTime, FLOAT _pos[MAX_DIMENSION])
{

  int i;

  CurrentGrid = _grid;    
  NextParticle = NULL;
  PrevParticle = NULL;                 
  UniqueID = _ID;
  GroupID  = _groupID;
  Level = _level;
  ExchangeCount = 0;
  CreationTime = _creationTime;
  Mass = 0.0;
  EchangedThisTimestep = false;

  // user defined attributes (e.g. max temperature)
  ParticleAttributes = new float[NumberOfParticleAttributes];
  for (i = 0; i < NumberOfParticleAttributes; i++)
    ParticleAttributes[i] = 0.0;

  for (i = 0; i < MAX_DIMENSION; i++)
    InitialPosition[i] = _pos[i];
  
  // history of particle position, including the time when the position was recorded. 
  #ifdef TRACK_MC_HISTORY
  LagrangianHistory = NULL; /* **** TODO **** */
  #endif
}

// MonteCarloTracerParticle::MonteCarloTracerParticle(MonteCarloTracerParticleBuffer *buffer, int n) 
// {
//   // TODO
// }

// MonteCarloTracerParticle::MonteCarloTracerParticle(MonteCarloTracerParticleBuffer buffer) 
// {
//   // TODO
// }


MonteCarloTracerParticle::~MonteCarloTracerParticle(void)
{
  if (ParticleAttributes != NULL)
      delete [] ParticleAttributes;

  #ifdef TRACK_MC_HISTORY 
  MonteeCarloTracerParticleHistory *temp = &LagrangianHistory;
  MonteeCarloTracerParticleHistory *next = NULL;
  while (temp != NULL) {
    next = temp->NextFrame;
    delete temp;
    temp = next;
  }
  #endif  

  NextParticle = NULL;
  PrevParticle = NULL;
  CurrentGrid  = NULL;  
}

// /***************

//     OPERATORS

//  ***************/

// void MonteCarloTracerParticle::operator=(MonteCarloTracerParticle a)
// {
//   // TODO
//   return;
// }

// /**********************

//    CONVENIENT ROUTINES

//  **********************/

// MonteCarloTracerParticle *MonteCarloTracerParticle::copy(void)
// {
//   int i, dim;
//   MonteCarloTracerParticle *a = new MonteCarloTracerParticle;

//   // TODO: Copy all member data
  
//   return a;
// }

// void MonteCarloTracerParticle::ConvertMassToSolar(void)
// {
//   double dx;
//   float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
//     VelocityUnits, MassConversion;
//   GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
// 	   &TimeUnits, &VelocityUnits, CurrentGrid->Time);
//   dx = LengthUnits * CurrentGrid->CellWidth[0][0];
//   MassConversion = (float) (dx*dx*dx * double(DensityUnits) / SolarMass);
//   this->Mass *= MassConversion;
//   return;
// }

// void MonteCarloTracerParticle::CopyToGrid()
// {

//   // This is cloned from Star::CopyToGrid and is surely wrong!!!
//   // Need to think about this... Is it needed, how do i need to change it?

//   MonteCarloTracerParticle *cParticle;
//   if (CurrentGrid != NULL)   // NULL => On another processor
//     for (cParticle = CurrentGrid->MonteCarloTracerParticles; cParticle; cParticle = cParticle->NextMonteCarloTracerParticle)
//       if (Identifier == cParticle->Identifier) {
// 	     *cParticle = *this;
// 	     break;
//       } // ENDIF match
//   return;
// }

// void MonteCarloTracerParticle::CopyFromParticle(grid *_grid, int _id, int _level)
// {
//   // TODO
//   return;
// }

// void MonteCarloTracerParticle::DeleteCopyInGrid(void)
// {

//   // This is cloned from Star::DeleteCopyInGrid and may be wrong.
//   // This works on a single particle list. Since a grid has lists 
//   // for each cell, this will need to be called for every cell.

//   MonteCarloTracerParticle *cParticle, *MoveParticle;
//   if (CurrentGrid != NULL) {   // NULL => On another processor
//     cParticle = CurrentGrid->MonteCarloTracerParticles;
//     CurrentGrid->MonteCarloTracerParticles = NULL;
//     while (cParticle) {
//       MoveParticle = PopMonteCarloTracerParticle(cParticle);
//       if (Identifier == MoveParticle->Identifier)
// 	delete MoveParticle;
//       else
// 	InsertMonteCarloTracerParticleAfter(CurrentGrid->MonteCarloTracerParticles, MoveParticle);
//     } // ENDWHILE MonteCarloTracerParticles
//   } // ENDIF grid != NULL
//   return;
// }

// void MonteCarloTracerParticle::PrintInfo(void)
// {
//   // TODO
// }


// /**************************************************
//      CONVERSION ROUTINES FROM/TO ARRAY BUFFERS
// **************************************************/

// void MonteCarloTracerParticle::MonteCarloTracerParticleListToBuffer(MonteCarloTracerParticleBuffer *&result, int n)
// {
//   // TODO
//   return;
// }

// void MonteCarloTracerParticle::MonteCarloTracerParticleToBuffer(MonteCarloTracerParticleBuffer *result)
// {
//   // TODO
//   return;
// }

