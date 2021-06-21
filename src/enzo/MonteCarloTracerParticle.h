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
/  Monte Carlo Tracers in each cell are stored as a singly-linked list.
/
************************************************************************/

#include "typedefs.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#

class MonteCarloTracerParticle
{

  private:
    grid      *CurrentGrid;                  
    PINT       Identifier;                    // unique identifier
    int        Level;
    int        CellIndicies[MAX_DIMENSION];   // cell i,j,k where this particle is
    int        ExchangeCount;                 // number of cell exchanges this particle has experienced
    float      CreationTime;
    float      mass;
    // user defined attributes (e.g. max temperature, group id)
    float     *ParticleAttributes[NumberOfMonteCarloTracerAttributes];
    FLOAT      InitialPosition[MAX_DIMENSION];
    bool       EchangedThisTimestep;
    
    // history of particle position, including the time when the position was recorded. 
    // looks like [[x0,y0,z0,t0], [x1,y1,z1,t1], ...]
    #ifdef TRACK_MC_POSITION
    FLOAT     *LagrangianHistory[MAX_DIMENSION + 1]; 
    #endif 

    friend class grid;

  public:

    MonteCarloTracerParticle *NextParticle; 

    // Constructors and destructor
    MonteCarloTracer();
    ~MonteCarloTracer();

    // Operators

    // Routines

    // getters
    // setters
  
};

