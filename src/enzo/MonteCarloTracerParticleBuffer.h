/***********************************************************************
/
/  MONTE CARLO TRACER PARTICLE STRUCTURE FOR COMMUNICATION
/
/  written by: Corey Brummel-Smith
/  date:       November, 2021
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#ifndef __MONTECARLOTRACERPARTICLEBUFFER_H
#define __MONTECARLOTRACERPARTICLEBUFFER_H


struct MonteCarloTracerParticleBuffer {               
  PINT       UniqueID;
  PINT       GroupID;
  int        Level;
  int        ExchangeCount;          // number of cell exchanges this particle has experienced
  float      CreationTime;
  float      Mass;
  // user defined attributes (e.g. max temperature, max velocity, etc.)
  float     *ParticleAttributes;
  FLOAT      InitialPosition[MAX_DIMENSION];
  bool       ExchangedThisTimestep;
  FLOAT      Position[MAX_DIMENSION];
};

#endif
