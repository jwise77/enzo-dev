 /*-*-C++-*-*/
/***********************************************************************
/
/  MOTE CARLO TRACER PARTICLES
/
/  written by: Corey Brummel-Smith
/  date:       May, 2021
/
/  PURPOSE:
/  Tracer particles that live in grid cells and are exchanged from cell
/  to cell based on a probability determined by the cell mass and mass
/  flux. Unlike velocity tracer particles, Monte Carlo tracers do not
/  have lagrangian phase space coordinates. They can be thought of as a
/  property of the grid cell. 
/
/  Monte Carlo Tracers in each cell are stored as a linked list.
/
************************************************************************/

#ifndef MONTE_CARLO_TRACER_PARTICLE_DEFINED__
#define MONTE_CARLO_TRACER_PARTICLE_DEFINED__

#include "typedefs.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

class MonteCarloTracerParticle
{

 private:

  grid      *CurrentGrid;                  
  PINT       UniqueID;
  PINT       GroupID;
  int        Level;
  int        ExchangeCount;          // number of cell exchanges this particle has experienced
  float      CreationTime;
  float      Mass;
  // user defined attributes (e.g. max temperature, max velocity, etc.)
  float     *ParticleAttributes;
  FLOAT      InitialPosition[MAX_DIMENSION];
  FLOAT      Position[MAX_DIMENSION];
  bool       ExchangedThisTimestep;
  int        WillDelete;
  
  // history of particle position, including the time when the position was recorded. (NOT IMPLEMENTED)
  #ifdef TRACK_MC_HISTORY
  struct MonteCarloTracerParticleHistory
  {
      FLOAT Position[3];
      float Timestamp;
      MonteCarloTracerParticleHistory *NextFrame;
  };     
  MonteCarloTracerParticleHistory *LagrangianHistory;
  #endif 
  
  friend class grid;
  friend class ExternalBoundary;
  
 public:

  MonteCarloTracerParticle *NextParticle;
  MonteCarloTracerParticle *PrevParticle; 
  
  // Constructors and destructor
  MonteCarloTracerParticle();
  MonteCarloTracerParticle(grid *_grid, int _ID, int _groupID, int _level, float _creationTime, FLOAT _pos[MAX_DIMENSION]);
  MonteCarloTracerParticle(grid *_grid, int _ID, int _groupID, int _level, float _creationTime, FLOAT _pos[MAX_DIMENSION], 
                           float _mass, int _exchange_count);
  MonteCarloTracerParticle(const MonteCarloTracerParticle* mc); 
  MonteCarloTracerParticle(MonteCarloTracerParticleBuffer buffer);
  ~MonteCarloTracerParticle();

  // Operators
  
  // Routines
  void MonteCarloTracerParticleToBuffer(MonteCarloTracerParticleBuffer *result, FLOAT* pos);
  
  // getters
  
  // setters
  
};

#endif