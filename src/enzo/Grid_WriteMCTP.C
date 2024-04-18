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

// void grid::WriteMCTP(char* prename)
// {
//   int i, j, k, dim, field, size, active_size, ActiveDim[MAX_DIMENSION];

//   float *temp;
//   hid_t       group_id, dset_id;
//   hid_t       file_id;
//   hsize_t     OutDims[MAX_DIMENSION];
//   hsize_t     FullOutDims[MAX_DIMENSION];
//   hsize_t     TempIntArray[1];
//   herr_t      h5_status;
//   herr_t      h5_error = -1;

//   char pid[MAX_TASK_TAG_SIZE];
//   sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
 
//   char gpid[MAX_TASK_TAG_SIZE];
//   sprintf(gpid, "%"TASK_TAG_FORMAT""ISYM, ProcessorNumber);

//   char id[MAX_GROUP_TAG_SIZE];
//   sprintf(id, "%"GROUP_TAG_FORMAT""ISYM, this->ID);  
 
//   char *base_name = new char[MAX_LINE_LENGTH];
//   strcpy(base_name, prename);
//   strcat(base_name, "_MCTP_data");

//   char *groupfilename = new char[MAX_LINE_LENGTH];
//   strcpy(groupfilename, base_name);
//   strcat(groupfilename, ".pid");
//   strcat(groupfilename, pid);
//   strcat(groupfilename, ".gpid");
//   strcat(groupfilename, gpid);

 
//   char *name = new char[MAX_LINE_LENGTH];
//   strcpy(name, "/Grid");
//   strcat(name, id);

//   file_id = H5Fcreate(groupfilename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
//   group_id = H5Gcreate(file_id, name, 0);


//    for (dim = 0; dim < GridRank; dim++) {
//      OutDims[GridRank-dim-1] = ActiveDim[dim];
//      FullOutDims[GridRank-dim-1] = GridDimension[dim];
//    }

//    size = 1;
 
//    for (dim = 0; dim < GridRank; dim++) size *= GridDimension[dim];
 
//    /* create temporary buffer */
 
//    temp = new float[size];


//    int NumberOfMCTracers = this->CountMonteCarloTracerParticles();
//    printf("\nproc%0d: WriteMCTP: grid %d, ProcessorNumber %d, NMCTP: %d\n", MyProcessorNumber, this->ID, ProcessorNumber, NumberOfMCTracers);
//    fflush(stdout);

//    if (NumberOfMCTracers > 0) {
   
//    int i, j, k, index, dim, n_per_cell, size = 1, n = 0;
//    FLOAT index_position[GridRank];
//    MonteCarloTracerParticle *mc;

//    for (dim = 0; dim < GridRank; dim++)
//      size *= GridDimension[dim];   

//    temp = new float[size];
//    float NumberOfMCTracersPerCell[size]; /* float because write_dataset needs HDF5_REAL for 3d 
//                                             data when ACTIVE_ONLY == 1 */ 
//    float position_x[size], position_y[size], position_z[size];
//    float part_pos_x[size], part_pos_y[size], part_pos_z[size];
//    float init_position_x[size], init_position_y[size], init_position_z[size];

//    TempIntArray[0] = NumberOfMCTracers;

//    /* Compute particle position (all particles in cell center) */
//    for (k = 0; k < GridDimension[2]; k++) {      
//      index_position[2] = CellLeftEdge[2][k] + 0.5 * CellWidth[2][0];
//    for (j = 0; j < GridDimension[1]; j++) {
//      index_position[1] = CellLeftEdge[1][j] + 0.5 * CellWidth[1][0];
//    for (i = 0; i < GridDimension[0]; i++) {
//      index_position[0] = CellLeftEdge[0][i] + 0.5 * CellWidth[0][0];
//      index = i + GridDimension[0]*(j + GridDimension[1]*k);
//      mc = MonteCarloTracerParticles[index];
//      n_per_cell = 0;

//      if (mc != NULL) {
//        position_x[index] = index_position[0];
//        position_y[index] = index_position[1];
//        position_z[index] = index_position[2];
//        part_pos_x[index] = mc->Position[0];
//        part_pos_y[index] = mc->Position[1];
//        part_pos_z[index] = mc->Position[2];         
       
//     }
//     else {
//        position_x[index] = NULL;
//        position_y[index] = NULL;
//        position_z[index] = NULL; 
//        part_pos_x[index] = NULL;
//        part_pos_y[index] = NULL;
//        part_pos_z[index] = NULL;          
//     }     

//      /* Store particle properties in single arrays */
//      while(mc != NULL){

//        mc = mc->NextParticle;
//        n++;
//        n_per_cell++;
//      }
//      NumberOfMCTracersPerCell[index] = n_per_cell;


//    } // end i
//    } // end j
//    } // end k






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
  strcat(groupfilename, ".");
  strcat(groupfilename, id);

 
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

   int NumberOfMCTracers = this->CountMonteCarloTracerParticles();
   //int NumberOfMCTracers = GetNumberOfMonteCarloTracerParticles();
   printf("\nWMC pid%s gpid%s NMC %d\n", pid, gpid, NumberOfMCTracers);

    if (NumberOfMCTracers > 0) {

    // allocate space all particles on this grid
    int   MonteCarloTracerExchangeCount[NumberOfMCTracers];
    PINT  MonteCarloTracerGroupID[NumberOfMCTracers];
    PINT  MonteCarloTracerUniqueID[NumberOfMCTracers];
    float MonteCarloTracerMass[NumberOfMCTracers];
    float MonteCarloTracerCreationTime[NumberOfMCTracers];
    FLOAT MonteCarloTracerInitialPosition_x[NumberOfMCTracers];
    FLOAT MonteCarloTracerInitialPosition_y[NumberOfMCTracers];
    FLOAT MonteCarloTracerInitialPosition_z[NumberOfMCTracers];
    FLOAT IndexPosition_x[NumberOfMCTracers];
    FLOAT IndexPosition_y[NumberOfMCTracers];
    FLOAT IndexPosition_z[NumberOfMCTracers];   
    FLOAT MCPosition_x[NumberOfMCTracers];
    FLOAT MCPosition_y[NumberOfMCTracers];
    FLOAT MCPosition_z[NumberOfMCTracers];       
    
    int i, j, k, index, dim, n_per_cell, size = 1, n = 0;
    FLOAT pos[GridRank];
    MonteCarloTracerParticle *mc;

    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];   

   /* create temporary buffer */

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

        // if (i < GridStartIndex[0] || i > GridEndIndex[0]
        //     j < GridStartIndex[1] || j > GridEndIndex[1]
        //     k < GridStartIndex[2] || k > GridEndIndex[2]){
        //   // DEBUG
        //   printf("\nMCinGZ (i,j,k)=(%d,%d,%d), pos=(%f,%f,%f)", i,j,k, pos[0], pos[1], pos[2]);
        // }

        MonteCarloTracerExchangeCount[n]     =  mc->ExchangeCount;
        MonteCarloTracerGroupID[n]           =  mc->GroupID;
        MonteCarloTracerUniqueID[n]          =  mc->UniqueID;
        MonteCarloTracerMass[n]              =  mc->Mass;
        MonteCarloTracerCreationTime[n]      =  mc->CreationTime;
        MonteCarloTracerInitialPosition_x[n] =  mc->InitialPosition[0];
        MonteCarloTracerInitialPosition_y[n] =  mc->InitialPosition[1];
        MonteCarloTracerInitialPosition_z[n] =  mc->InitialPosition[2];
        IndexPosition_x[n]        =  pos[0];
        IndexPosition_y[n]        =  pos[1];
        IndexPosition_z[n]        =  pos[2];
        MCPosition_x[n]        =  mc->Position[0];
        MCPosition_y[n]        =  mc->Position[1];
        MCPosition_z[n]        =  mc->Position[2];                   
        mc = mc->NextParticle;
        n++;
        n_per_cell++;
      }
      NumberOfMCTracersPerCell[index] = n_per_cell;
    } // end i
    } // end j
    } // end k

    /* Write MC particle data to hdf5 */

    // this->write_dataset(1, TempIntArray, "MonteCarloTracerExchangeCount",
    //     group_id, HDF5_PINT, (VOIDP) MonteCarloTracerExchangeCount, FALSE);
    // this->write_dataset(1, TempIntArray, "MonteCarloTracerGroupID",
    //     group_id, HDF5_PINT, (VOIDP) MonteCarloTracerGroupID, FALSE);
    // this->write_dataset(1, TempIntArray, "MonteCarloTracerUniqueID",
    //     group_id, HDF5_PINT, (VOIDP) MonteCarloTracerUniqueID, FALSE);
    // this->write_dataset(1, TempIntArray, "MonteCarloTracerMass",
    //     group_id, HDF5_REAL, (VOIDP) MonteCarloTracerMass, FALSE);
    // this->write_dataset(1, TempIntArray, "MonteCarloTracerCreationTime",
    //     group_id, HDF5_REAL, (VOIDP) MonteCarloTracerCreationTime, FALSE);
    this->write_dataset(1, TempIntArray, "MonteCarloTracerInitialPosition_x",
        group_id, HDF5_REAL, (VOIDP) MonteCarloTracerInitialPosition_x, FALSE);
    this->write_dataset(1, TempIntArray, "MonteCarloTracerInitialPosition_y",
        group_id, HDF5_REAL, (VOIDP) MonteCarloTracerInitialPosition_y, FALSE);
    this->write_dataset(1, TempIntArray, "MonteCarloTracerInitialPosition_z",
        group_id, HDF5_REAL, (VOIDP) MonteCarloTracerInitialPosition_z, FALSE);
    this->write_dataset(1, TempIntArray, "IndexPosition_x",
        group_id, HDF5_REAL, (VOIDP) IndexPosition_x, FALSE);
    this->write_dataset(1, TempIntArray, "IndexPosition_y",
        group_id, HDF5_REAL, (VOIDP) IndexPosition_y, FALSE);
    this->write_dataset(1, TempIntArray, "IndexPosition_z",
        group_id, HDF5_REAL, (VOIDP) IndexPosition_z, FALSE);  
    this->write_dataset(1, TempIntArray, "MCPosition_x",
        group_id, HDF5_REAL, (VOIDP) MCPosition_x, FALSE);
    this->write_dataset(1, TempIntArray, "MCPosition_y",
        group_id, HDF5_REAL, (VOIDP) MCPosition_y, FALSE);
    this->write_dataset(1, TempIntArray, "MCPosition_z",
        group_id, HDF5_REAL, (VOIDP) MCPosition_z, FALSE);      
    // this->write_dataset(GridRank, OutDims, "NumberOfMonteCarloTracersPerCell",
    //     group_id, HDF5_REAL, (VOIDP) NumberOfMCTracersPerCell,
    //     CopyOnlyActive, temp); 
    this->write_dataset(GridRank, FullOutDims, "NumberOfMonteCarloTracersPerCell",
        group_id, HDF5_REAL, (VOIDP) NumberOfMCTracersPerCell,
        FALSE, temp);                                      
}
  h5_status = H5Gclose(group_id); 
  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);} 
  h5_status = H5Fclose(file_id);
  if( h5_status == h5_error ){my_exit(EXIT_FAILURE);}

}




