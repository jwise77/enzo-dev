/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A STAR PARTICLE TEST)
/
/  written by: Greg Bryan
/  date:       June, 2012
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "ActiveParticle.h"
#include "ActiveParticle_SmartStar.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::TestStarParticleInitializeGrid(float TestStarParticleStarMass, 
					 float *Initialdt,
					 FLOAT TestStarParticleStarVelocity[],
           FLOAT TestStarParticleStarPosition[],
           int TestStarParticleUseSmartStar,
           float TestStarParticleSmartStarAge)
{
  /* declarations */

  float CentralMass = 1.0;
  int i, dim;
  float TestInitialdt = *Initialdt;
  
  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  
  /* Get Units. */

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1;
  double MassUnits = 1;
  
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Set Central Mass in simulation units */

  CentralMass = TestStarParticleStarMass*1.99e33* pow(LengthUnits*CellWidth[0][0],-3.0)/DensityUnits;

  printf("Central Mass: %f \n",CentralMass);

  /* Set number of particles for this grid and allocate space. */

  if (TestStarParticleUseSmartStar == 0) {
    NumberOfParticles = 1;
    NumberOfParticleAttributes = 4;
    this->AllocateNewParticles(NumberOfParticles);
    printf("Allocated %d particles\n", NumberOfParticles);

    /* Set particle IDs and types */

    for (i = 0; i < NumberOfParticles; i++) {
      ParticleNumber[i] = i;
      ParticleType[i] = PARTICLE_TYPE_STAR;
    }

    /* Set central particle. */ 
    for (dim = 0; dim < GridRank; dim++) {
      ParticlePosition[dim][0] = TestStarParticleStarPosition[dim]*
        (DomainLeftEdge[dim]+DomainRightEdge[dim]) + 0.5*CellWidth[0][0];
      ParticleVelocity[dim][0] = TestStarParticleStarVelocity[dim]*1e5*TimeUnits/LengthUnits;
    }
    ParticleMass[0] = CentralMass;
    ParticleAttribute[0][0] = Time+1e-7;

    if (STARFEED_METHOD(UNIGRID_STAR)) ParticleAttribute[1][0] = 10.0 * Myr_s/TimeUnits;
    if (STARFEED_METHOD(MOM_STAR))
      if(StarMakerExplosionDelayTime >= 0.0)
        ParticleAttribute[1][0] = 1.0;
      else
        ParticleAttribute[1][0] = 10.0 * Myr_s/TimeUnits;
    
    ParticleAttribute[2][0] = 0.0;  // Metal fraction
    ParticleAttribute[3][0] = 0.0;  // metalfSNIa

  } else {
    /* Smart Star creation */

    
    float AccretionRadius = 3.0;  // in cell widths

    /* Create SmartStar and set all initial properties */
    ActiveParticleType_SmartStar *np = new ActiveParticleType_SmartStar();
    np->level = 0;
    np->GridID = 0;
    np->CurrentGrid = this;
    np->Mass = CentralMass;
    np->oldmass = CentralMass;
    np->BirthTime = Time;
    np->DynamicalTime = 0.0;
    np->type = np->GetEnabledParticleID();
    np->Metallicity = 0.0;
    for (dim = 0; dim < GridRank; dim++) {
      np->pos[dim] = TestStarParticleStarPosition[dim] * (DomainLeftEdge[dim]+DomainRightEdge[dim]) + 
        0.5*CellWidth[0][0];
      np->vel[dim] = TestStarParticleStarVelocity[dim]*1e5*TimeUnits/LengthUnits;
    }
    np->AccretionRadius = AccretionRadius * CellWidth[0][0];
    np->StellarAge = TestStarParticleSmartStarAge * Myr_s / TimeUnits;
    np->NotEjectedMass = 0.0;
    for (i = 0; i < 2; i++) {
      np->AccretionRate[i] = 0.0;
      np->AccretionRateTime[i] = Time + i*1e-6;
    }
    np->TimeIndex = 1;
    np->RadiationLifetime = 20 * Myr_s / TimeUnits;  // 20 Myr

    // Add SmartStar to grid AP list
    this->AddActiveParticle(np);
    delete np;  // AddActiveParticle copies it to the grid, so we delete this copy

  }

  return SUCCESS;
}

