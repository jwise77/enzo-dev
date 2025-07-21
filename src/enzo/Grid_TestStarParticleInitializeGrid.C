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

int FindField(int field, int farray[], int numfields);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::TestStarParticleInitializeGrid(float TestStarParticleStarMass, 
					 float *Initialdt,
					 FLOAT TestStarParticleStarVelocity[],
           FLOAT TestStarParticleStarPosition[],
           int TestStarParticleIsothermalSphere,
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

  /* Set isothermal sphere density profile.  Default diameter of the box size. */

  if (TestStarParticleIsothermalSphere) {
    const float SphereRadius = 0.5;
    const float min_r2 = 0.25*CellWidth[0][0]*CellWidth[0][0];
    const float max_r2 = 0.25*(DomainRightEdge[0] - DomainLeftEdge[0])*(DomainRightEdge[0] - DomainLeftEdge[0]);
    int i, j, k, index, DensNum, field;
    index = GRIDINDEX(0,0,0);
    DensNum = FindField(Density, this->FieldType, this->NumberOfBaryonFields);
    const float BackgroundDensity = BaryonField[DensNum][index];
    float dx, dy, dz, dx2, dy2, dz2, dr2, dr2_0, factor;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      dz = TestStarParticleStarPosition[2] - (CellLeftEdge[2][k] + 0.5*CellWidth[2][k]);
      dz2 = dz*dz;
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        dy = TestStarParticleStarPosition[1] - (CellLeftEdge[1][j] + 0.5*CellWidth[1][j]);
        dy2 = dy*dy;
        dr2_0 = dy2 + dz2;
        index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
          dx = TestStarParticleStarPosition[0] - (CellLeftEdge[0][i] + 0.5*CellWidth[0][i]);
          dx2 = dx*dx;
          dr2 = max(dr2_0 + dx2, min_r2);
          if (dr2 < max_r2) {
            factor = max_r2/dr2;
            for (field = 0; field < NumberOfBaryonFields; field++) {
              if (FieldTypeIsDensity(field)) {
                BaryonField[field][index] *= factor;
              }
            }
          }
        } // ENDFOR i
      } // ENDFOR j
    } // ENDFOR k
  } // ENDIF isothermal sphere

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
    np->InfluenceRadius = AccretionRadius * CellWidth[0][0];
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

