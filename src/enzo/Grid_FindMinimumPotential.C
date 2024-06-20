/***********************************************************************
/
/  RETURN MINIMUM OF GRAVITATIONAL POTENTIAL IN A CONTROL VOLUME
/
/  written by: John Regan
/  date:       January 2017
/  modified1:
/
/  PURPOSE: Used to test for minimum of the potential
/
************************************************************************/
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "phys_constants.h"

float grid::FindMinimumPotential(FLOAT *cellpos, FLOAT radius, float *PotentialField)
{

  if (PotentialField == NULL)
	ENZO_FAIL("PotentialField not defined.");

  int i = 0, j = 0, k = 0;
  float dx, dy, dz;
  float dx2, dy2, dz2;
  FLOAT r2 = radius*radius;
  FLOAT pos[3];
  float GravitationalMinimum = 1e20;
  int index;
 
  /* Need to find gravitational minimum of this grid in advance */
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	pos[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	dz = cellpos[2] - pos[2];
	dz2 = dz*dz;
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
		pos[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
		dy = cellpos[1] - pos[1];
		dy2 = dy*dy;
		index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
      	for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
			pos[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
			dx = cellpos[0] - pos[0];
			dx2 = dx*dx;
			if (PotentialField[index] != PotentialField[index])
	  			PotentialField[index] = 0.0;
			/* Use only values within the control volume */
			if ((dx2 + dy2 + dz2 <= r2) && (PotentialField[index] < GravitationalMinimum))
	    		GravitationalMinimum = PotentialField[index];
    	}
    }
  }

  return GravitationalMinimum;
}

/*
 * Calculate the Potential Field in this grid by direct summation
 * This is an N^2 problem - be careful how frequently you call it. 
 * Phi = -G m_j / (x - x_j)
 */
void grid::CalculatePotentialField(float *PotentialField, int DensNum, float DensityUnits,
				   float TimeUnits, float LengthUnits)
{
  int i = 0, j = 0, k = 0, ii = 0, jj = 0, kk = 0;

  FLOAT dx, dy, dz, dx2, dy2, dz2;
  FLOAT sep = 0;
  FLOAT dV = pow(CellWidth[0][0], 3);
  FLOAT pos[3], pos_j[3];
  //float MassUnits = DensityUnits*LengthUnits*LengthUnits*LengthUnits;
  //float G = 4*M_PI*GravConst*DensityUnits*TimeUnits*TimeUnits;
  float G = GravitationalConstant;
  float mass_j = 0;
  int index = 0, lindex = 0;
  float *density = BaryonField[DensNum];
  
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
	pos[2] = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	  pos[1] = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  index = GRIDINDEX_NOGHOST(GridStartIndex[0], j, k);
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
		pos[0] = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	

	for (kk = GridStartIndex[2]; kk <= GridEndIndex[2]; kk++) {
      pos_j[2] = CellLeftEdge[2][kk] + 0.5*CellWidth[2][kk];
	  dz = pos[2] - pos_j[2];
	  dz2 = dz*dz;
	  for (jj = GridStartIndex[1]; jj <= GridEndIndex[1]; jj++) {
        pos_j[1] = CellLeftEdge[1][jj] + 0.5*CellWidth[1][jj];
		dy = pos[1] - pos_j[1];
		dy2 = dy*dy;
        lindex = GRIDINDEX_NOGHOST(GridStartIndex[0], jj, kk);
	    for (ii = GridStartIndex[0]; ii <= GridEndIndex[0]; ii++, lindex++) {
	      mass_j = density[lindex]*dV;
	      pos_j[0] = CellLeftEdge[0][ii] + 0.5*CellWidth[0][ii];
		  dx = pos[0] - pos_j[0];
		  dx2 = dx*dx;
	      sep = sqrt(dx2 + dy2 + dz2);
	      if(index == lindex) {
			sep = CellWidth[0][ii]/2.0;
	      }
	      PotentialField[index] += -G * mass_j/sep;
	    }
	  }
	}
    }
    }
  }
  return;
}
