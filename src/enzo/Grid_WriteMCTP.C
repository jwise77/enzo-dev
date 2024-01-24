#include <hdf5.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "h5utilities.h"
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CommunicationUtilities.h"

void my_exit(int status);

void grid::WriteMCTP(char* prename)
{
  int i, j, k, dim, field, size, active_size, ActiveDim[MAX_DIMENSION];

  float *temp;
  hid_t       group_id, dset_id;
  hid_t       file_id;
  hsize_t     OutDims[MAX_DIMENSION];
  hsize_t     FullOutDims[MAX_DIMENSION];
  hsize_t     TempIntArray[1];
  herr_t      h5_status;
  herr_t      h5_error = -1;

  char pid[MAX_TASK_TAG_SIZE];
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
 
  char gpid[MAX_TASK_TAG_SIZE];
  sprintf(gpid, "%"TASK_TAG_FORMAT""ISYM, ProcessorNumber);

  char id[MAX_GROUP_TAG_SIZE];
  sprintf(id, "%"GROUP_TAG_FORMAT""ISYM, this->ID);  
 
  char *base_name = new char[MAX_LINE_LENGTH];
  strcpy(base_name, prename);
  strcat(base_name, "_MCTP_data");

  char *groupfilename = new char[MAX_LINE_LENGTH];
  strcpy(groupfilename, base_name);
  strcat(groupfilename, ".pid");
  strcat(groupfilename, pid);
  strcat(groupfilename, ".gpid");
  strcat(groupfilename, gpid);

 
  char *name = new char[MAX_LINE_LENGTH];
  strcpy(name, "/Grid");
  strcat(name, id);

  file_id = H5Fcreate(groupfilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  group_id = H5Gcreate(file_id, name, 0);


   for (dim = 0; dim < GridRank; dim++) {
     OutDims[GridRank-dim-1] = ActiveDim[dim];
     FullOutDims[GridRank-dim-1] = GridDimension[dim];
   }

   size = 1;
 
   for (dim = 0; dim < GridRank; dim++) size *= GridDimension[dim];
 
   /* create temporary buffer */
 
   temp = new float[size];


   int NumberOfMCTracers = this->CountMonteCarloTracerParticles();
   printf("\nproc%0d: WriteMCTP: grid %d, ProcessorNumber %d, NMCTP: %d\n", MyProcessorNumber, this->ID, ProcessorNumber, NumberOfMCTracers);
   fflush(stdout);

   if (NumberOfMCTracers > 0) {
   
   int i, j, k, index, dim, n_per_cell, size = 1, n = 0;
   FLOAT pos[GridRank];
   MonteCarloTracerParticle *mc;

   for (dim = 0; dim < GridRank; dim++)
     size *= GridDimension[dim];   

   temp = new float[size];
   float NumberOfMCTracersPerCell[size]; /* float because write_dataset needs HDF5_REAL for 3d 
                                            data when ACTIVE_ONLY == 1 */ 

   TempIntArray[0] = NumberOfMCTracers;

   /* Compute particle position (all particles in cell center) */
   for (k = 0; k < GridDimension[2]; k++) {      
     pos[2] = CellLeftEdge[2][k] + 0.5 * CellWidth[2][0];
   for (j = 0; j < GridDimension[1]; j++) {
     pos[1] = CellLeftEdge[1][j] + 0.5 * CellWidth[1][0];
   for (i = 0; i < GridDimension[0]; i++) {
     pos[0] = CellLeftEdge[0][i] + 0.5 * CellWidth[0][0];
     index = i + GridDimension[0]*(j + GridDimension[1]*k);
     mc = MonteCarloTracerParticles[index];
     n_per_cell = 0;

     /* Store particle properties in single arrays */
     while(mc != NULL){

       mc = mc->NextParticle;
       n++;
       n_per_cell++;
     }
     NumberOfMCTracersPerCell[index] = n_per_cell;
   } // end i
   } // end j
   } // end k

   this->write_dataset(GridRank, FullOutDims, "NumberOfMonteCarloTracersPerCell",
       group_id, HDF5_REAL, (VOIDP) NumberOfMCTracersPerCell,
       FALSE, temp);           
}
  h5_status = H5Gclose(group_id); 
  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);} 
  h5_status = H5Fclose(file_id);
  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}
}
