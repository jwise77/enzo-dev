/***********************************************************************
/
/  INITIALIZE A TRIGGERED STAR FORMATION SIMULATION
/
/  written by: Corey Brummel-Smith
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


static float TSFInitialFractionHII   = 1.2e-5;
static float TSFInitialFractionHeII  = 1.0e-14;
static float TSFInitialFractionHeIII = 1.0e-17;
static float TSFInitialFractionHM    = 2.0e-9;
static float TSFInitialFractionH2I   = 2.0e-20;
static float TSFInitialFractionH2II  = 3.0e-14;

int TriggeredStarFormationInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
       HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  const char *DensName = "Density";
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  const char *Vel3Name = "z-velocity";
  const char *ColourName = "colour";
  const char *ElectronName = "Electron_Density";
  const char *HIName    = "HI_Density";
  const char *HIIName   = "HII_Density";
  const char *HeIName   = "HeI_Density";
  const char *HeIIName  = "HeII_Density";
  const char *HeIIIName = "HeIII_Density";
  const char *HMName    = "HM_Density";
  const char *H2IName   = "H2I_Density";
  const char *H2IIName  = "H2II_Density";
  const char *DIName    = "DI_Density";
  const char *DIIName   = "DII_Density";
  const char *HDIName   = "HDI_Density";
  const char *kphHIName    = "HI_kph";
  const char *gammaName  = "PhotoGamma";
  const char *kphHeIName   = "HeI_kph";   
  const char *kphHeIIName  = "HeII_kph";
  const char *kphHMName   = "HM_kphs";   
  const char *kdissH2IIName  = "H2II_kdiss";
  const char *kdissH2IName = "H2I_kdiss"; 
  const char *RadAccel1Name = "x-RadPressure";
  const char *RadAccel2Name = "y-RadPressure";
  const char *RadAccel3Name = "z-RadPressure";


  /* declarations */

  char  line[MAX_LINE_LENGTH];
  char *dummy = new char[MAX_LINE_LENGTH];
  int   dim, ret, level, i, source;
  dummy[0] = 0;

  char *TSFDensityFilename = NULL,
    *TSFHIIFractionFilename = NULL,
    *TSFHeIIFractionFilename = NULL,
    *TSFHeIIIFractionFilename = NULL,
    *TSFTemperatureFilename = NULL;
  int TSFUseParticles    = FALSE; /*** LOOK AT WHAT THIS DOES ***/
  int TSFUseColour       = FALSE; 
  int TestProblemUseMetallicityField;
  float TestProblemInitialMetallicityFraction;
 

  /* uniform background params */

  float TSFUniformDensity;
  float TSFUniformEnergy;
  float TSFUniformVelocity[MAX_DIMENSION];
  float TSFUniformBField[3];


  /* star params */

  float TSFStarParticleStarVelocity[3];
  FLOAT TSFStarPosition[3];
  float TSFStarMass;

  /* spherical cloud params */

  float TSFInitialTemperature = 1000;
  int   TSFSphereType,
        TSFSphereConstantPressure,
        TSFSphereSmoothSurface,
        TSFSphereNumShells;
  float TSFSphereDensity,
        TSFSphereTemperature,
        TSFSphereVelocity[MAX_DIMENSION],
        TSFFracKeplerianRot,
        TSFSphereTurbulence,
        TSFSphereCutOff,
        TSFSphereAng1,
        TSFSphereAng2,
        TSFSphereSmoothRadius,
        TSFSphereRadius,
        TSFSphereCoreRadius,
        TSFSphereHIIFraction,
        TSFSphereHeIIFraction,
        TSFSphereHeIIIFraction,
        TSFSphereH2IFraction;
  FLOAT TSFSpherePosition[3];

  float TSFOmegaBaryonNow=0.05; // *** NOT NEEDED *** 

  /* set default parameters */

  rewind(fptr);

  if (debug) {
    fprintf(stderr, "TSFInitialize: Set up test problem.\n");
  }

  // Set default values

  TSFUniformDensity     = 1.0;
  TSFUniformEnergy      = 1.0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    TSFUniformVelocity[dim] = 0;
    TSFUniformBField[dim] = 0;  
    TSFStarParticleStarVelocity[dim] = 0;     
  }
  TSFSphereRadius     = 0.5;
  TSFSphereCoreRadius = 0.1;
  TSFSphereDensity    = 1.0;
  TSFSphereTemperature = 1.0;
  TSFFracKeplerianRot = 0.0;
  TSFSphereTurbulence = 0.0;
  TSFSphereCutOff = 6.5;
  TSFSphereAng1 = 0;
  TSFSphereAng2 = 0;
  TSFSphereNumShells = 1;
  TSFSphereSmoothRadius = 1.2;
  TSFSphereHIIFraction = TSFInitialFractionHII;
  TSFSphereHeIIFraction = TSFInitialFractionHeII;
  TSFSphereHeIIIFraction = TSFInitialFractionHeIII;
  TSFSphereH2IFraction = TSFInitialFractionH2I;

  TSFSpherePosition[0] = 0.5;
  TSFSpherePosition[1] = 0.5;
  TSFSpherePosition[2] = 0.5;

  TSFSphereType       = 0;
  TSFSphereConstantPressure = FALSE;
  TSFSphereSmoothSurface = FALSE;
  
  TSFStarPosition[0] = 0.5;
  TSFStarPosition[1] = 0.5;
  TSFStarPosition[2] = 0.0;
  TSFStarMass    = 100.0;

  TestProblemUseMetallicityField = 1;
  TestProblemInitialMetallicityFraction = 1e-11; // Zero metallicity
  

  TestProblemData.MultiSpecies =    MultiSpecies;
  TestProblemData.UseMetallicityField = TestProblemUseMetallicityField;
  TestProblemData.MetallicityField_Fraction = TestProblemInitialMetallicityFraction;

  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "TSFUniformDensity  = %"FSYM, &TSFUniformDensity);
    ret += sscanf(line, "TSFUniformEnergy   = %"FSYM, &TSFUniformEnergy);
    ret += sscanf(line, "TSFUniformVelocity = %"FSYM" %"FSYM" %"FSYM, 
          &TSFUniformVelocity[0],
          &TSFUniformVelocity[1],
          &TSFUniformVelocity[2]);
    ret += sscanf(line,"TSFStarVelocity = %"PSYM" %"PSYM" %"PSYM, 
          &TSFStarVelocity[0],
          &TSFStarVelocity[1],
          &TSFStarVelocity[2]);
    ret += sscanf(line,"TSFStarPosition = %"PSYM" %"PSYM" %"PSYM, 
          &TSFStarPosition[0],
          &TSFStarPosition[1],
          &TSFStarPosition[2]);
    ret += sscanf(line, "TSFStarMass = %"FSYM, &TSFStarMass);

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction); 
    ret += sscanf(line, "TestProblemInitialHIFraction      = %"FSYM, &TestProblemData.HI_Fraction);
    ret += sscanf(line, "TestProblemInitialHIIFraction     = %"FSYM, &TestProblemData.HII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIFraction     = %"FSYM, &TestProblemData.HeI_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIFraction    = %"FSYM, &TestProblemData.HeII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIIIFraction  = %"FSYM, &TestProblemData.HeIII_Fraction);
    ret += sscanf(line, "TestProblemHydrogenFractionByMass = %"FSYM, &TestProblemData.HydrogenFractionByMass);


    /* read cloud parameters */

    ret += sscanf(line, "TSFUseParticles = %"ISYM,  &TSFUseParticles);
    ret += sscanf(line, "TSFUseColour    = %"ISYM,  &TSFUseColour);
    ret += sscanf(line, "TSFInitialTemperature = %"FSYM,  &TSFInitialTemperature);
    if (sscanf(line, "TSFDensityFilename = %s", dummy) == 1) {
      ret++;
      TSFDensityFilename = new char[MAX_LINE_LENGTH];
      strcpy(TSFDensityFilename, dummy);
    }
    if (sscanf(line, "TSFHIIFractionFilename = %s", dummy) == 1) {
      ret++;
      TSFHIIFractionFilename = new char[MAX_LINE_LENGTH];
      strcpy(TSFHIIFractionFilename, dummy);
    }
    if (sscanf(line, "TSFHeIIFractionFilename = %s", dummy) == 1) {
      ret++;
      TSFHeIIFractionFilename = new char[MAX_LINE_LENGTH];
      strcpy(TSFHeIIFractionFilename, dummy);
    }
    if (sscanf(line, "TSFHeIIIFractionFilename = %s", dummy) == 1) {
      ret++;
      TSFHeIIIFractionFilename = new char[MAX_LINE_LENGTH];
      strcpy(TSFHeIIIFractionFilename, dummy);
    }
    if (sscanf(line, "TSFTemperatureFilename = %s", dummy) == 1) {
      ret++;
      TSFTemperatureFilename = new char[MAX_LINE_LENGTH];
      strcpy(TSFTemperatureFilename, dummy);
    }
    ret += sscanf(line, "TSFSphereType             = %"ISYM, &TSFSphereType);
    ret += sscanf(line, "TSFSphereConstantPressure = %"ISYM, &TSFSphereConstantPressure);
    ret += sscanf(line, "TSFSphereSmoothSurface    = %"ISYM, &TSFSphereSmoothSurface);
    ret += sscanf(line, "TSFSphereSmoothRadius     = %"FSYM, &TSFSphereSmoothRadius);
    ret += sscanf(line, "TSFSphereRadius           = %"FSYM, &TSFSphereRadius);
    ret += sscanf(line, "TSFSphereCoreRadius       = %"FSYM, &TSFSphereCoreRadius);
    ret += sscanf(line, "TSFSphereDensity          = %"FSYM, &TSFSphereDensity);
    ret += sscanf(line, "TSFSphereTemperature      = %"FSYM, &TSFSphereTemperature);
    ret += sscanf(line, "TSFFracKeplerianRot       = %"FSYM, &TSFFracKeplerianRot);
    ret += sscanf(line, "TSFSphereTurbulence       = %"FSYM, &TSFSphereTurbulence);
    ret += sscanf(line, "TSFSphereCutOff           = %"FSYM, &TSFSphereCutOff);
    ret += sscanf(line, "TSFSphereAng1             = %"FSYM, &TSFSphereAng1);
    ret += sscanf(line, "TSFSphereAng2             = %"FSYM, &TSFSphereAng2);
    ret += sscanf(line, "TSFSphereNumShells        = %"ISYM, &TSFSphereNumShells);
    ret += sscanf(line, "TSFSphereHIIFraction      = %"FSYM, &TSFSphereHIIFraction);
    ret += sscanf(line, "TSFSphereHeIIFraction     = %"FSYM, &TSFSphereHeIIFraction);
    ret += sscanf(line, "TSFSphereHeIIIFraction    = %"FSYM, &TSFSphereHeIIIFraction);
    ret += sscanf(line, "TSFSphereH2IFraction      = %"FSYM, &TSFSphereH2IFraction);
    ret += sscanf(line, "PhotonTimeStep            = %"FSYM, &dtPhoton);
    ret += sscanf(line, "TSFOmegaBaryonNow         = %"FSYM, &TSFOmegaBaryonNow);
    ret += sscanf(line, "TSFInitialFractionHII     = %"FSYM, &TSFInitialFractionHII);
    ret += sscanf(line, "TSFInitialFractionHeII    = %"FSYM, &TSFInitialFractionHeII);
    ret += sscanf(line, "TSFInitialFractionHeIII   = %"FSYM, &TSFInitialFractionHeIII);
    ret += sscanf(line, "TSFInitialFractionHM      = %"FSYM, &TSFInitialFractionHM);
    ret += sscanf(line, "TSFInitialFractionH2I     = %"FSYM, &TSFInitialFractionH2I);
    ret += sscanf(line, "TSFInitialFractionH2II    = %"FSYM, &TSFInitialFractionH2II);
    ret += sscanf(line, "TSFSpherePosition         = %"PSYM" %"PSYM" %"PSYM, 
      &TSFSpherePosition[0],
      &TSFSpherePosition[1],
      &TSFSpherePosition[2]);
    ret += sscanf(line, "TSFSphereVelocity         = %"FSYM" %"FSYM" %"FSYM, 
      &TSFSphereVelocity[0],
      &TSFSphereVelocity[1],
      &TSFSphereVelocity[2]);    

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "StarParticle") 
    && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up uniform grid */

  if (TopGrid.GridData->InitializeUniformGrid(TSFUniformDensity, 
                          TSFUniformEnergy,
                          TSFUniformEnergy,
                          TSFUniformVelocity,
                          TSFUniformBField) == FAIL)
    ENZO_FAIL("Error in InitializeUniformGrid.");

  /* set up star and spherical cloud */

  if (TopGrid.GridData->TSFInitializeGrid(
      TSFSphereRadius, TSFSphereCoreRadius, TSFSphereDensity,
      TSFSphereTemperature, TSFSpherePosition, TSFSphereVelocity,
      TSFFracKeplerianRot, TSFSphereTurbulence, TSFSphereCutOff,
      TSFSphereAng1, TSFSphereAng2, TSFSphereNumShells,
      TSFSphereType, TSFSphereConstantPressure, TSFSphereSmoothSurface,
      TSFSphereSmoothRadius, TSFSphereHIIFraction, TSFSphereHeIIFraction,
      TSFSphereHeIIIFraction, TSFSphereH2IFraction, TSFUseParticles,
      TSFUseColour, TSFInitialTemperature, 0,
      TSFInitialFractionHII, TSFInitialFractionHeII, TSFInitialFractionHeIII,
      TSFInitialFractionHM, TSFInitialFractionH2I, TSFInitialFractionH2II,
      TSFDensityFilename, TSFHIIFractionFilename, TSFHeIIFractionFilename,
      TSFHeIIIFractionFilename, TSFTemperatureFilename, 
      TSFStarMass, TSFStarVelocity, TSFStarPosition,
      Initialdt)
        ENZO_FAIL("Error in TriggeredStarFormationInitializeGrid.\n"));

  
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
  if (TSFUseColour)
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
    fprintf(Outfptr, "PhotonTestUseParticles       = %"ISYM"\n",
      PhotonTestUseParticles);
    fprintf(Outfptr, "PhotonTestUseColour          = %"ISYM"\n",
      PhotonTestUseColour);
    fprintf(Outfptr, "PhotonTestInitialTemperature = %"FSYM"\n",
      PhotonTestInitialTemperature);

      fprintf(Outfptr, "PhotonTestSphereType             = %"ISYM"\n",             PhotonTestSphereType);
      fprintf(Outfptr, "PhotonTestSphereConstantPressure = %"ISYM"\n", PhotonTestSphereConstantPressure);
      fprintf(Outfptr, "PhotonTestSphereSmoothSurface    = %"ISYM"\n",    PhotonTestSphereSmoothSurface);
      fprintf(Outfptr, "PhotonTestSphereSmoothRadius     = %"GOUTSYM"\n",  PhotonTestSphereSmoothRadius);
      fprintf(Outfptr, "PhotonTestSphereRadius           = %"GOUTSYM"\n",        PhotonTestSphereRadius);
      fprintf(Outfptr, "PhotonTestSphereCoreRadius       = %"GOUTSYM"\n",    PhotonTestSphereCoreRadius);
      fprintf(Outfptr, "PhotonTestSphereDensity          = %"FSYM"\n",          PhotonTestSphereDensity);
      fprintf(Outfptr, "PhotonTestSphereTemperature      = %"FSYM"\n",      PhotonTestSphereTemperature);
      fprintf(Outfptr, "PhotonTestSpherePosition         = ");
      WriteListOfFloats(Outfptr, MetaData.TopGridRank, PhotonTestSpherePosition);
      fprintf(Outfptr, "PhotonTestSphereVelocity         = ");
      WriteListOfFloats(Outfptr, MetaData.TopGridRank, PhotonTestSphereVelocity);
      fprintf(Outfptr, "PhotonTestSphereHIIFraction      = %"GOUTSYM"\n",   PhotonTestSphereHIIFraction);
      fprintf(Outfptr, "PhotonTestSphereHeIIFraction     = %"GOUTSYM"\n",  PhotonTestSphereHeIIFraction);
      fprintf(Outfptr, "PhotonTestSphereHeIIIFraction    = %"GOUTSYM"\n", PhotonTestSphereHeIIIFraction);
      fprintf(Outfptr, "PhotonTestSphereH2IFraction      = %"GOUTSYM"\n",   PhotonTestSphereH2IFraction);
      fprintf(Outfptr, "PhotonTestFracKeplerianRot       = %"GOUTSYM"\n",    PhotonTestFracKeplerianRot);
      fprintf(Outfptr, "PhotonTestSphereTurbulence       = %"GOUTSYM"\n",    PhotonTestSphereTurbulence);
      fprintf(Outfptr, "PhotonTestSphereCutOff           = %"GOUTSYM"\n",        PhotonTestSphereCutOff);
      fprintf(Outfptr, "PhotonTestSphereAng1             = %"GOUTSYM"\n",          PhotonTestSphereAng1);
      fprintf(Outfptr, "PhotonTestSphereAng2             = %"GOUTSYM"\n",          PhotonTestSphereAng2);
      fprintf(Outfptr, "PhotonTestSphereNumShells        = %"ISYM"\n\n",      PhotonTestSphereNumShells);

    fprintf(Outfptr, "TSFUniformDensity = %"FSYM"\n",
      TSFUniformDensity);
    fprintf(Outfptr, "TSFUniformEnergy = %"FSYM"\n",
      TSFUniformEnergy);
    fprintf(Outfptr, "MetallicityField_Fraction = %"FSYM"\n",
            TestProblemData.MetallicityField_Fraction);

    fprintf(Outfptr, "PhotonTestUniformVelocity          = %"FSYM" %"FSYM" %"FSYM"\n");
    WriteListOfFloats(Outfptr, MetaData.TopGridRank, PhotonTestUniformVelocity);
    
  }

  return SUCCESS;
}

