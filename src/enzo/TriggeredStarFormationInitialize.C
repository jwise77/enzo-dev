/***********************************************************************
/
/  INITIALIZE A TRIGGERED STAR FORMATION SIMULATION
/  
/  based on PhotonTestInitialize(1) and TestStarParticleInitialize(2)
/  (1) written by: Tom Abel
/  date:       Oct 2003
/  (2) written by: Greg Bryan
/  date:       June, 2012

/  modified:   Corey Brummel-Smith
/  date:       October, 2019
/
/  PURPOSE:
/    Initialize a Pop III star particle at the end of its life, near 
/    a sphereical cloud, embedded in a uniform medium. The supernova 
/    from the Pop III star may potentially crush the cloud and induce 
/    gravitational collapse.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
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

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);

static float TSF_InitialFractionHII   = 1.2e-5;
static float TSF_InitialFractionHeII  = 1.0e-14;
static float TSF_InitialFractionHeIII = 1.0e-17;
static float TSF_InitialFractionHM    = 2.0e-9;
static float TSF_InitialFractionH2I   = 2.0e-20;
static float TSF_InitialFractionH2II  = 3.0e-14;

int TriggeredStarFormationInitialize(FILE *fptr, FILE *Outfptr,
       HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  const char *DensName      = "Density";
  const char *TEName        = "TotalEnergy";
  const char *GEName        = "GasEnergy";
  const char *Vel1Name      = "x-velocity";
  const char *Vel2Name      = "y-velocity";
  const char *Vel3Name      = "z-velocity";
  const char *ColourName    = "colour";
  const char *ElectronName  = "Electron_Density";
  const char *HIName        = "HI_Density";
  const char *HIIName       = "HII_Density";
  const char *HeIName       = "HeI_Density";
  const char *HeIIName      = "HeII_Density";
  const char *HeIIIName     = "HeIII_Density";
  const char *HMName        = "HM_Density";
  const char *H2IName       = "H2I_Density";
  const char *H2IIName      = "H2II_Density";
  const char *DIName        = "DI_Density";
  const char *DIIName       = "DII_Density";
  const char *HDIName       = "HDI_Density";
  const char *kphHIName     = "HI_kph";
  const char *gammaName     = "PhotoGamma";
  const char *kphHeIName    = "HeI_kph";   
  const char *kphHeIIName   = "HeII_kph";
  const char *kphHMName     = "HM_kphs";   
  const char *kdissH2IIName = "H2II_kdiss";
  const char *kdissH2IName  = "H2I_kdiss"; 
  const char *RadAccel1Name = "x-RadPressure";
  const char *RadAccel2Name = "y-RadPressure";
  const char *RadAccel3Name = "z-RadPressure";


  /* declarations */

  char  line[MAX_LINE_LENGTH];
  char *dummy = new char[MAX_LINE_LENGTH];
  int   dim, ret, i, source;
  dummy[0] = 0;

  char *TSF_DensityFilename       = NULL,
       *TSF_HIIFractionFilename   = NULL,
       *TSF_HeIIFractionFilename  = NULL,
       *TSF_HeIIIFractionFilename = NULL,
       *TSF_TemperatureFilename   = NULL;
  int   TSF_UseColour       = TRUE; 
 

  /* uniform background params */

  float TSF_UniformDensity;
  float TSF_UniformTemperature;
  float TSF_UniformVelocity[MAX_DIMENSION];
  float TSF_UniformBField[MAX_DIMENSION];

  /* star params */

  float TSF_StarVelocity[MAX_DIMENSION];
  FLOAT TSF_StarPosition[MAX_DIMENSION];
  float TSF_StarMass;
  int   TSF_TimeToExplosion;

  /* spherical cloud params */

  int   TSF_SphereType,
        TSF_SphereConstantPressure,
        TSF_SphereSmoothSurface,
        TSF_SphereNumShells;
  float TSF_SphereDensity,
        TSF_SphereTemperature,
        TSF_SphereVelocity[MAX_DIMENSION],
        TSF_FracKeplerianRot,
        TSF_SphereTurbulence,
        TSF_SphereCutOff,
        TSF_SphereAng1,
        TSF_SphereAng2,
        TSF_SphereSmoothRadius,
        TSF_SphereRadius,
        TSF_SphereCoreRadius,
        TSF_SphereHIIFraction,
        TSF_SphereHeIIFraction,
        TSF_SphereHeIIIFraction,
        TSF_SphereH2IFraction;
  FLOAT TSF_SpherePosition[MAX_DIMENSION];

  /* set default parameters */

  rewind(fptr);

  if (debug) {
    fprintf(stderr, "TSF_Initialize: Set up test problem.\n");
  }

  // Set default values

  TSF_UniformDensity = 1.0;
  TSF_UniformTemperature = 1000;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    TSF_UniformVelocity[dim] = 0;
    TSF_UniformBField[dim] = 0;  
    TSF_StarVelocity[dim] = 0;     
  }
  TSF_SphereRadius     = 0.5;
  TSF_SphereCoreRadius = 0.1;
  TSF_SphereDensity    = 1.0;
  TSF_SphereTemperature = 1.0;
  TSF_FracKeplerianRot = 0.0;
  TSF_SphereTurbulence = 0.0;
  TSF_SphereCutOff = 6.5;
  TSF_SphereAng1 = 0;
  TSF_SphereAng2 = 0;
  TSF_SphereNumShells = 1;
  TSF_SphereSmoothRadius = 1.2;
  TSF_SphereHIIFraction = TSF_InitialFractionHII;
  TSF_SphereHeIIFraction = TSF_InitialFractionHeII;
  TSF_SphereHeIIIFraction = TSF_InitialFractionHeIII;
  TSF_SphereH2IFraction = TSF_InitialFractionH2I;

  TSF_SpherePosition[0] = 0.5;
  TSF_SpherePosition[1] = 0.5;
  TSF_SpherePosition[2] = 0.5;

  TSF_SphereType       = 0;
  TSF_SphereConstantPressure = FALSE;
  TSF_SphereSmoothSurface = FALSE;
  
  TSF_StarPosition[0] = 0.5;
  TSF_StarPosition[1] = 0.5;
  TSF_StarPosition[2] = 0.0;
  TSF_StarMass        = 100.0; // Msun
  TSF_TimeToExplosion = 0.0;   // kyr

  if (!StarParticleFeedback)
    StarParticleFeedback = pow(2, POP3_STAR);


  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "TSF_UniformDensity     = %"FSYM, &TSF_UniformDensity);
    ret += sscanf(line, "TSF_UniformTemperature = %"FSYM,  &TSF_UniformTemperature);    
    ret += sscanf(line, "TSF_UniformVelocity    = %"FSYM" %"FSYM" %"FSYM, 
          &TSF_UniformVelocity[0],
          &TSF_UniformVelocity[1],
          &TSF_UniformVelocity[2]);
    ret += sscanf(line,"TSF_StarVelocity = %"PSYM" %"PSYM" %"PSYM, 
          &TSF_StarVelocity[0],
          &TSF_StarVelocity[1],
          &TSF_StarVelocity[2]);
    ret += sscanf(line,"TSF_StarPosition = %"PSYM" %"PSYM" %"PSYM, 
          &TSF_StarPosition[0],
          &TSF_StarPosition[1],
          &TSF_StarPosition[2]);
    ret += sscanf(line, "TSF_StarMass           = %"FSYM, &TSF_StarMass);
    ret += sscanf(line, "TSF_TimeToExplosion    = %"ISYM,  &TSF_UseColour);    
    ret += sscanf(line, "TSF_UseColour          = %"ISYM,  &TSF_UseColour);

    /* read cloud parameters */

    if (sscanf(line, "TSF_DensityFilename = %s", dummy) == 1) {
      ret++;
      TSF_DensityFilename = new char[MAX_LINE_LENGTH];
      strcpy(TSF_DensityFilename, dummy);
    }
    if (sscanf(line, "TSF_HIIFractionFilename = %s", dummy) == 1) {
      ret++;
      TSF_HIIFractionFilename = new char[MAX_LINE_LENGTH];
      strcpy(TSF_HIIFractionFilename, dummy);
    }
    if (sscanf(line, "TSF_HeIIFractionFilename = %s", dummy) == 1) {
      ret++;
      TSF_HeIIFractionFilename = new char[MAX_LINE_LENGTH];
      strcpy(TSF_HeIIFractionFilename, dummy);
    }
    if (sscanf(line, "TSF_HeIIIFractionFilename = %s", dummy) == 1) {
      ret++;
      TSF_HeIIIFractionFilename = new char[MAX_LINE_LENGTH];
      strcpy(TSF_HeIIIFractionFilename, dummy);
    }
    if (sscanf(line, "TSF_TemperatureFilename = %s", dummy) == 1) {
      ret++;
      TSF_TemperatureFilename = new char[MAX_LINE_LENGTH];
      strcpy(TSF_TemperatureFilename, dummy);
    }
    ret += sscanf(line, "TSF_SphereType             = %"ISYM, &TSF_SphereType);
    ret += sscanf(line, "TSF_SphereConstantPressure = %"ISYM, &TSF_SphereConstantPressure);
    ret += sscanf(line, "TSF_SphereSmoothSurface    = %"ISYM, &TSF_SphereSmoothSurface);
    ret += sscanf(line, "TSF_SphereSmoothRadius     = %"FSYM, &TSF_SphereSmoothRadius);
    ret += sscanf(line, "TSF_SphereRadius           = %"FSYM, &TSF_SphereRadius);
    ret += sscanf(line, "TSF_SphereCoreRadius       = %"FSYM, &TSF_SphereCoreRadius);
    ret += sscanf(line, "TSF_SphereDensity          = %"FSYM, &TSF_SphereDensity);
    ret += sscanf(line, "TSF_SphereTemperature      = %"FSYM, &TSF_SphereTemperature);
    ret += sscanf(line, "TSF_FracKeplerianRot       = %"FSYM, &TSF_FracKeplerianRot);
    ret += sscanf(line, "TSF_SphereTurbulence       = %"FSYM, &TSF_SphereTurbulence);
    ret += sscanf(line, "TSF_SphereCutOff           = %"FSYM, &TSF_SphereCutOff);
    ret += sscanf(line, "TSF_SphereAng1             = %"FSYM, &TSF_SphereAng1);
    ret += sscanf(line, "TSF_SphereAng2             = %"FSYM, &TSF_SphereAng2);
    ret += sscanf(line, "TSF_SphereNumShells        = %"ISYM, &TSF_SphereNumShells);
    ret += sscanf(line, "TSF_SphereHIIFraction      = %"FSYM, &TSF_SphereHIIFraction);
    ret += sscanf(line, "TSF_SphereHeIIFraction     = %"FSYM, &TSF_SphereHeIIFraction);
    ret += sscanf(line, "TSF_SphereHeIIIFraction    = %"FSYM, &TSF_SphereHeIIIFraction);
    ret += sscanf(line, "TSF_SphereH2IFraction      = %"FSYM, &TSF_SphereH2IFraction);
    ret += sscanf(line, "TSF_InitialFractionHII     = %"FSYM, &TSF_InitialFractionHII);
    ret += sscanf(line, "TSF_InitialFractionHeII    = %"FSYM, &TSF_InitialFractionHeII);
    ret += sscanf(line, "TSF_InitialFractionHeIII   = %"FSYM, &TSF_InitialFractionHeIII);
    ret += sscanf(line, "TSF_InitialFractionHM      = %"FSYM, &TSF_InitialFractionHM);
    ret += sscanf(line, "TSF_InitialFractionH2I     = %"FSYM, &TSF_InitialFractionH2I);
    ret += sscanf(line, "TSF_InitialFractionH2II    = %"FSYM, &TSF_InitialFractionH2II);
    ret += sscanf(line, "TSF_SpherePosition         = %"PSYM" %"PSYM" %"PSYM, 
      &TSF_SpherePosition[0],
      &TSF_SpherePosition[1],
      &TSF_SpherePosition[2]);
    ret += sscanf(line, "TSF_SphereVelocity         = %"FSYM" %"FSYM" %"FSYM, 
      &TSF_SphereVelocity[0],
      &TSF_SphereVelocity[1],
      &TSF_SphereVelocity[2]);    

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "StarParticle") 
    && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up uniform grid, star particle, and spherical cloud */

  if (TopGrid.GridData->TriggeredStarFormationInitializeGrid(
      TSF_UniformDensity, TSF_UniformTemperature,
      TSF_UniformVelocity, TSF_UniformBField,
      TSF_SphereRadius, TSF_SphereCoreRadius, TSF_SphereDensity,
      TSF_SphereTemperature, TSF_SpherePosition, TSF_SphereVelocity,
      TSF_FracKeplerianRot, TSF_SphereTurbulence, TSF_SphereCutOff,
      TSF_SphereAng1, TSF_SphereAng2, TSF_SphereNumShells,
      TSF_SphereType, TSF_SphereConstantPressure, TSF_SphereSmoothSurface,
      TSF_SphereSmoothRadius, TSF_SphereHIIFraction, TSF_SphereHeIIFraction,
      TSF_SphereHeIIIFraction, TSF_SphereH2IFraction, TSF_UseColour,
      TSF_InitialFractionHII, TSF_InitialFractionHeII, TSF_InitialFractionHeIII,
      TSF_InitialFractionHM, TSF_InitialFractionH2I, TSF_InitialFractionH2II,
      TSF_DensityFilename, TSF_HIIFractionFilename, TSF_HeIIFractionFilename,
      TSF_HeIIIFractionFilename, TSF_TemperatureFilename, 
      TSF_StarMass, TSF_StarPosition, TSF_StarVelocity, TSF_TimeToExplosion))
        ENZO_FAIL("Error in TriggeredStarFormationInitializeGrid.\n");

  
 /* set up field names and units */

  int count = 0;
  DataLabel[count++] = (char*) DensName;
  DataLabel[count++] = (char*) TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = (char*) GEName;
  DataLabel[count++] = (char*) Vel1Name;
  DataLabel[count++] = (char*) Vel2Name;
  DataLabel[count++] = (char*) Vel3Name;
  if (MultiSpecies) {
    DataLabel[count++] = (char*) ElectronName;
    DataLabel[count++] = (char*) HIName;
    DataLabel[count++] = (char*) HIIName;
    DataLabel[count++] = (char*) HeIName;
    DataLabel[count++] = (char*) HeIIName;
    DataLabel[count++] = (char*) HeIIIName;
    if (MultiSpecies > 1) {
      DataLabel[count++] = (char*) HMName;
      DataLabel[count++] = (char*) H2IName;
      DataLabel[count++] = (char*) H2IIName;
    }
    if (MultiSpecies > 2) {
      DataLabel[count++] = (char*) DIName;
      DataLabel[count++] = (char*) DIIName;
      DataLabel[count++] = (char*) HDIName;
    }
  }  // if Multispecies
  if (TSF_UseColour)
    DataLabel[count++] = (char*) ColourName;
  
  if (RadiativeTransfer)
    if (MultiSpecies) {
      DataLabel[count++]  = (char*) kphHIName;
      DataLabel[count++]  = (char*) gammaName;
      DataLabel[count++]  = (char*) kphHeIName;
      DataLabel[count++]  = (char*) kphHeIIName;
      if (MultiSpecies > 1) {
  DataLabel[count++]= (char*) kdissH2IName;
  DataLabel[count++]= (char*) kdissH2IIName;
  DataLabel[count++]= (char*) kphHMName;
      }
    } // if RadiativeTransfer

  if (RadiationPressure) {
    DataLabel[count++]  = (char*) RadAccel1Name;
    DataLabel[count++]  = (char*) RadAccel2Name;
    DataLabel[count++]  = (char*) RadAccel3Name;
  }

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "TSF_UseColour          = %"ISYM"\n", TSF_UseColour);
    fprintf(Outfptr, "TSF_SphereType             = %"ISYM"\n",             TSF_SphereType);
    fprintf(Outfptr, "TSF_SphereConstantPressure = %"ISYM"\n", TSF_SphereConstantPressure);
    fprintf(Outfptr, "TSF_SphereSmoothSurface    = %"ISYM"\n",    TSF_SphereSmoothSurface);
    fprintf(Outfptr, "TSF_SphereSmoothRadius     = %"GOUTSYM"\n",  TSF_SphereSmoothRadius);
    fprintf(Outfptr, "TSF_SphereRadius           = %"GOUTSYM"\n",        TSF_SphereRadius);
    fprintf(Outfptr, "TSF_SphereCoreRadius       = %"GOUTSYM"\n",    TSF_SphereCoreRadius);
    fprintf(Outfptr, "TSF_SphereDensity          = %"FSYM"\n",          TSF_SphereDensity);
    fprintf(Outfptr, "TSF_SphereTemperature      = %"FSYM"\n",      TSF_SphereTemperature);
    fprintf(Outfptr, "TSF_SpherePosition         = ");
    WriteListOfFloats(Outfptr, MetaData.TopGridRank, TSF_SpherePosition);
    fprintf(Outfptr, "TSF_SphereVelocity         = ");
    WriteListOfFloats(Outfptr, MetaData.TopGridRank, TSF_SphereVelocity);
    fprintf(Outfptr, "TSF_SphereHIIFraction      = %"GOUTSYM"\n",   TSF_SphereHIIFraction);
    fprintf(Outfptr, "TSF_SphereHeIIFraction     = %"GOUTSYM"\n",  TSF_SphereHeIIFraction);
    fprintf(Outfptr, "TSF_SphereHeIIIFraction    = %"GOUTSYM"\n", TSF_SphereHeIIIFraction);
    fprintf(Outfptr, "TSF_SphereH2IFraction      = %"GOUTSYM"\n",   TSF_SphereH2IFraction);
    fprintf(Outfptr, "TSF_FracKeplerianRot       = %"GOUTSYM"\n",    TSF_FracKeplerianRot);
    fprintf(Outfptr, "TSF_SphereTurbulence       = %"GOUTSYM"\n",    TSF_SphereTurbulence);
    fprintf(Outfptr, "TSF_SphereCutOff           = %"GOUTSYM"\n",        TSF_SphereCutOff);
    fprintf(Outfptr, "TSF_SphereAng1             = %"GOUTSYM"\n",          TSF_SphereAng1);
    fprintf(Outfptr, "TSF_SphereAng2             = %"GOUTSYM"\n",          TSF_SphereAng2);
    fprintf(Outfptr, "TSF_SphereNumShells        = %"ISYM"\n\n",      TSF_SphereNumShells);
    fprintf(Outfptr, "TSF_UniformDensity         = %"FSYM"\n",         TSF_UniformDensity);
    fprintf(Outfptr, "TSF_UniformTemperature     = %"FSYM"\n",     TSF_UniformTemperature);
    fprintf(Outfptr, "TSF_UniformVelocity        = %"FSYM" %"FSYM" %"FSYM"\n");
    WriteListOfFloats(Outfptr, MetaData.TopGridRank, TSF_UniformVelocity); 
    fprintf(Outfptr, "TSF_StarPosition           = %"FSYM" %"FSYM" %"FSYM"\n");
    WriteListOfFloats(Outfptr, MetaData.TopGridRank, TSF_StarPosition);     
    fprintf(Outfptr, "TSF_StarVelocity           = %"FSYM" %"FSYM" %"FSYM"\n");
    WriteListOfFloats(Outfptr, MetaData.TopGridRank, TSF_StarVelocity);     
    fprintf(Outfptr, "TSF_StarMass               = %"FSYM"\n",               TSF_StarMass);
    fprintf(Outfptr, "TSF_TimeToExplosion        = %"FSYM"\n",        TSF_TimeToExplosion);

  }

  delete[] TSF_DensityFilename;
  delete[] TSF_HIIFractionFilename;
  delete[] TSF_HeIIFractionFilename;
  delete[] TSF_HeIIIFractionFilename;
  delete[] TSF_TemperatureFilename;
  delete[] dummy;  

  return SUCCESS;
}

