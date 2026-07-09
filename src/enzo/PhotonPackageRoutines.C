/***********************************************************************
/
/  PHOTON PACKAGE ROUTINES
/
/  written by: John Wise
/  date:       February, 2010
/  modified1:  
/
/  PURPOSE: Constructs a Linked List of Photon Packages including data
/
************************************************************************/
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "PhotonPackage.h"

PhotonPackageEntry::PhotonPackageEntry(void)
{
  NextPackage = NULL;
  PreviousPackage = NULL;
  CurrentSource = NULL;
  Photons = 0.0;
  Type = 0;
  Energy = 0.0;
  CrossSection = 0.0;
  EmissionTimeInterval = 0.0;
  EmissionTime = 0.0;
  CurrentTime = 0.0;
  Radius = 0.0;
  ColumnDensity = 0.0;
  ipix = 0;
  level = 0;
  SourcePosition[0] = 0.0;
  SourcePosition[1] = 0.0;
  SourcePosition[2] = 0.0;
  SourcePositionDiff = 0.0;
}

/**********************************************************************/

#ifdef MEMORY_POOL
void* PhotonPackageEntry::operator new(size_t object_size)
{
  return PhotonMemoryPool->GetMemory(object_size);
}

void PhotonPackageEntry::operator delete(void* object)
{
  PhotonMemoryPool->FreeMemory(object);
}
#endif

/**********************************************************************/
/*                     PhotonPackageSoA IMPLEMENTATION                */
/**********************************************************************/

PhotonPackageSoA::PhotonPackageSoA(void)
{
  numPackages = 0;
  capacity = 0;
  Flux = NULL;
  Type = NULL;
  Energy = NULL;
  CrossSection = NULL;
  TimeInterval = NULL;
  EmissionTime = NULL;
  CurrentTime = NULL;
  Radius = NULL;
  ColumnDensity = NULL;
  PixelNum = NULL;
  Level = NULL;
  SourceX = NULL;
  SourceY = NULL;
  SourceZ = NULL;
  SourcePositionDiff = NULL;
  CurrentSource = NULL;
}

PhotonPackageSoA::~PhotonPackageSoA(void)
{
  this->free_arrays();
}

void PhotonPackageSoA::free_arrays(void)
{
  if (capacity > 0) {
    delete [] Flux;
    delete [] Type;
    delete [] Energy;
    delete [] CrossSection;
    delete [] TimeInterval;
    delete [] EmissionTime;
    delete [] CurrentTime;
    delete [] Radius;
    delete [] ColumnDensity;
    delete [] PixelNum;
    delete [] Level;
    delete [] SourceX;
    delete [] SourceY;
    delete [] SourceZ;
    delete [] SourcePositionDiff;
    delete [] CurrentSource;
  }
  Flux = NULL;
  Type = NULL;
  Energy = NULL;
  CrossSection = NULL;
  TimeInterval = NULL;
  EmissionTime = NULL;
  CurrentTime = NULL;
  Radius = NULL;
  ColumnDensity = NULL;
  PixelNum = NULL;
  Level = NULL;
  SourceX = NULL;
  SourceY = NULL;
  SourceZ = NULL;
  SourcePositionDiff = NULL;
  CurrentSource = NULL;
  numPackages = 0;
  capacity = 0;
}

void PhotonPackageSoA::initialize(int initial_capacity)
{
  this->free_arrays();
  if (initial_capacity > 0) {
    this->resize(initial_capacity);
  }
}

void PhotonPackageSoA::resize(int new_capacity)
{
  if (new_capacity <= capacity) return;

  float            *new_Flux = new float[new_capacity];
  int              *new_Type = new int[new_capacity];
  float            *new_Energy = new float[new_capacity];
  double           *new_CrossSection = new double[new_capacity];
  FLOAT            *new_TimeInterval = new FLOAT[new_capacity];
  FLOAT            *new_EmissionTime = new FLOAT[new_capacity];
  FLOAT            *new_CurrentTime = new FLOAT[new_capacity];
  FLOAT            *new_Radius = new FLOAT[new_capacity];
  float            *new_ColumnDensity = new float[new_capacity];
  int64_t          *new_PixelNum = new int64_t[new_capacity];
  int              *new_Level = new int[new_capacity];
  FLOAT            *new_SourceX = new FLOAT[new_capacity];
  FLOAT            *new_SourceY = new FLOAT[new_capacity];
  FLOAT            *new_SourceZ = new FLOAT[new_capacity];
  float            *new_SourcePositionDiff = new float[new_capacity];
  SuperSourceEntry **new_CurrentSource = new SuperSourceEntry*[new_capacity];

  if (numPackages > 0) {
    for (int i = 0; i < numPackages; i++) {
      new_Flux[i] = Flux[i];
      new_Type[i] = Type[i];
      new_Energy[i] = Energy[i];
      new_CrossSection[i] = CrossSection[i];
      new_TimeInterval[i] = TimeInterval[i];
      new_EmissionTime[i] = EmissionTime[i];
      new_CurrentTime[i] = CurrentTime[i];
      new_Radius[i] = Radius[i];
      new_ColumnDensity[i] = ColumnDensity[i];
      new_PixelNum[i] = PixelNum[i];
      new_Level[i] = Level[i];
      new_SourceX[i] = SourceX[i];
      new_SourceY[i] = SourceY[i];
      new_SourceZ[i] = SourceZ[i];
      new_SourcePositionDiff[i] = SourcePositionDiff[i];
      new_CurrentSource[i] = CurrentSource[i];
    }
  }

  if (capacity > 0) {
    delete [] Flux;
    delete [] Type;
    delete [] Energy;
    delete [] CrossSection;
    delete [] TimeInterval;
    delete [] EmissionTime;
    delete [] CurrentTime;
    delete [] Radius;
    delete [] ColumnDensity;
    delete [] PixelNum;
    delete [] Level;
    delete [] SourceX;
    delete [] SourceY;
    delete [] SourceZ;
    delete [] SourcePositionDiff;
    delete [] CurrentSource;
  }

  Flux = new_Flux;
  Type = new_Type;
  Energy = new_Energy;
  CrossSection = new_CrossSection;
  TimeInterval = new_TimeInterval;
  EmissionTime = new_EmissionTime;
  CurrentTime = new_CurrentTime;
  Radius = new_Radius;
  ColumnDensity = new_ColumnDensity;
  PixelNum = new_PixelNum;
  Level = new_Level;
  SourceX = new_SourceX;
  SourceY = new_SourceY;
  SourceZ = new_SourceZ;
  SourcePositionDiff = new_SourcePositionDiff;
  CurrentSource = new_CurrentSource;

  capacity = new_capacity;
}

void PhotonPackageSoA::append(const PhotonPackageEntry &PP)
{
  this->append(PP.Photons, PP.Type, PP.Energy, PP.CrossSection,
               PP.EmissionTimeInterval, PP.EmissionTime, PP.CurrentTime,
               PP.Radius, PP.ColumnDensity, PP.ipix, PP.level,
               PP.SourcePosition[0], PP.SourcePosition[1], PP.SourcePosition[2],
               PP.SourcePositionDiff, PP.CurrentSource);
}

void PhotonPackageSoA::append(float flux, int type, float energy, double cross_section,
                             FLOAT time_interval, FLOAT emission_time, FLOAT current_time,
                             FLOAT radius, float column_density, int64_t pixel_num, int level,
                             FLOAT source_x, FLOAT source_y, FLOAT source_z, float source_pos_diff,
                             SuperSourceEntry *source)
{
  if (numPackages >= capacity) {
    int new_capacity = (capacity == 0) ? 64 : capacity * 2;
    this->resize(new_capacity);
  }

  int idx = numPackages;
  Flux[idx] = flux;
  Type[idx] = type;
  Energy[idx] = energy;
  CrossSection[idx] = cross_section;
  TimeInterval[idx] = time_interval;
  EmissionTime[idx] = emission_time;
  CurrentTime[idx] = current_time;
  Radius[idx] = radius;
  ColumnDensity[idx] = column_density;
  PixelNum[idx] = pixel_num;
  Level[idx] = level;
  SourceX[idx] = source_x;
  SourceY[idx] = source_y;
  SourceZ[idx] = source_z;
  SourcePositionDiff[idx] = source_pos_diff;
  CurrentSource[idx] = source;

  numPackages++;
}

void PhotonPackageSoA::DeletePackage(int idx)
{
  if (idx < 0 || idx >= numPackages) return;
  if (idx < numPackages - 1) {
    int last = numPackages - 1;
    Flux[idx] = Flux[last];
    Type[idx] = Type[last];
    Energy[idx] = Energy[last];
    CrossSection[idx] = CrossSection[last];
    TimeInterval[idx] = TimeInterval[last];
    EmissionTime[idx] = EmissionTime[last];
    CurrentTime[idx] = CurrentTime[last];
    Radius[idx] = Radius[last];
    ColumnDensity[idx] = ColumnDensity[last];
    PixelNum[idx] = PixelNum[last];
    Level[idx] = Level[last];
    SourceX[idx] = SourceX[last];
    SourceY[idx] = SourceY[last];
    SourceZ[idx] = SourceZ[last];
    SourcePositionDiff[idx] = SourcePositionDiff[last];
    CurrentSource[idx] = CurrentSource[last];
  }
  numPackages--;
}

void PhotonPackageSoA::print_info(int idx)
{
  if (idx < 0 || idx >= numPackages) return;
  FLOAT r[3];
  double u[3];
  pix2vec_nest64((int64_t)(1 << Level[idx]), PixelNum[idx], u);
  for (int dim = 0; dim < 3; dim++) {
    FLOAT source_pos_dim = (dim == 0) ? SourceX[idx] : ((dim == 1) ? SourceY[idx] : SourceZ[idx]);
    r[dim] = source_pos_dim + u[dim] * Radius[idx];
  }
  printf("Photons = %g, Type = %" ISYM ", Radius = %" PSYM "\n", Flux[idx], Type[idx], Radius[idx]);
  printf("ipix = %lld, level = %d\n", (long long)PixelNum[idx], Level[idx]);
  printf("normal = %" PSYM " %" PSYM " %" PSYM "\n", u[0], u[1], u[2]);
  printf("SourcePosition = %" PSYM " %" PSYM " %" PSYM "\n",
	   SourceX[idx], SourceY[idx], SourceZ[idx]);
  printf("RayPosition = %" PSYM " %" PSYM " %" PSYM "\n", r[0], r[1], r[2]);
}

