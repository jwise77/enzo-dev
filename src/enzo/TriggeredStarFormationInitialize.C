/***********************************************************************
/
/  INITIALIZE A TRIGGERED STAR FORMATION SIMULATION
/
/  written by: Greg Bryan
/  date:       June, 2012
/  modified1:
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

int TriggeredStarFormationInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                   TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *MetalName = "Metal_Density";
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";


  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret;

  /* Error check. */


  /* set default parameters */

  float TSFUniformDensity     = 1.0;
  float TSFUniformEnergy      = 1.0;
  float TSFUniformVelocity[3] = {0.0, 0.0, 0TSFStarParticleStarVelocity[3] = {0.0, 0.0, 0.0};
  FLOAT TSFStarPosition[3] = {0.5, 0.5, 0.1};
  TSFUniformBField[3]   = {0.0, 0.0, 0.0};
  float TSFStarMass    = 100.0;
  int TestProblemUseMetallicityField = 1;
  float TestProblemInitialMetallicityFraction = 1e-11; // Zeri metallicity
  




  TestProblemData.MultiSpecies =    MultiSpecies;
  TestProblemData.UseMetallicityField = TestProblemUseMetallicityField;
  TestProblemData.MetallicityField_Fraction = TestProblemInitialMetallicityFraction;

  /* read input from file */


  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "TSFUniformDensity = %"FSYM,
          &TSFUniformDensity);
    ret += sscanf(line, "TSFUniformEnergy = %"FSYM,
          &TSFUniformEnergy);
    ret += sscanf(line, "TSFStarMass = %"FSYM,
          &TSFStarMass);
    ret += sscanf(line,"TSFStarVelocity = %"PSYM" %"PSYM" %"PSYM, 
          &TSFStarVelocity[0],
          &TSFStarVelocity[1],
          &TSFStarVelocity[2]);
    ret += sscanf(line,"TSFStarPosition = %"PSYM" %"PSYM" %"PSYM, 
          &TSFStarPosition[0],
          &TSFStarPosition[1],
          &TSFStarPosition[2]);
    

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction); 

    ret += sscanf(line, "TestProblemInitialHIFraction  = %"FSYM, &TestProblemData.HI_Fraction);
    ret += sscanf(line, "TestProblemInitialHIIFraction  = %"FSYM, &TestProblemData.HII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIFraction  = %"FSYM, &TestProblemData.HeI_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIFraction  = %"FSYM, &TestProblemData.HeII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIIIFraction  = %"FSYM, &TestProblemData.HeIII_Fraction);
    ret += sscanf(line, "TestProblemHydrogenFractionByMass = %"FSYM, &TestProblemData.HydrogenFractionByMass);


    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "StarParticle") 
    && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  /* set up uniform grid as of before explosion */

  
  
  if (TopGrid.GridData->InitializeUniformGrid(TSFUniformDensity, 
                          TSFUniformEnergy,
                          TSFUniformEnergy,
                          TSFUniformVelocity,
                          TSFUniformBField) == FAIL)
    ENZO_FAIL("Error in InitializeUniformGrid.");
 
  if (TopGrid.GridData->
      TriggeredStarFormationInitializeGrid(TSFStarMass,
                     Initialdt, 
                     TSFStarVelocity,
                     TSFStarPosition) == FAIL)
    ENZO_FAIL("Error in TriggeredStarFormationInitializeGrid.\n");

  /* set up field names and units */
  
  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[count++] = GEName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  if (TestProblemData.MultiSpecies){
    DataLabel[count++] = ElectronName;
    DataLabel[count++] = HIName;
    DataLabel[count++] = HIIName;
    DataLabel[count++] = HeIName;
    DataLabel[count++] = HeIIName;
    DataLabel[count++] = HeIIIName;
  }
  if (TestProblemData.UseMetallicityField)
    DataLabel[count++] = MetalName;


  int j;
  for(j=0; j < count; j++)
    DataUnits[j] = NULL;

   
  /* Write parameters to parameter output file */
  
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "TSFUniformDensity = %"FSYM"\n",
        TSFUniformDensity);
    fprintf(Outfptr, "TSFUniformEnergy = %"FSYM"\n",
        TSFUniformEnergy);
    fprintf(Outfptr, "MetallicityField_Fraction = %"FSYM"\n",
            TestProblemData.MetallicityField_Fraction);
  }

  fprintf(stderr, "TSFUniformDensity = %"FSYM"\n",
      TSFUniformDensity);
  fprintf(stderr, "TSFUniformEnergy = %"FSYM"\n",
      TSFUniformEnergy);
  fprintf(stderr, "MetallicityField_Fraction = %"FSYM"\n",
      TestProblemData.MetallicityField_Fraction);



  return SUCCESS;

}

