 /*-*-C++-*-*/
/***********************************************************************
/
/  MOTE CARLO TRACER PARTICLES
/
/  written by: Corey Brummel-Smith
/  date:       May, 2021
/
/  PURPOSE:
/  Tracer particles that live in grid cells and are exchanged from cell
/  to cell based on a probability determined by the cell mass and mass
/  flux. Unlike velocity tracer particles, Monte Carlo tracers do not
/  have lagrangian phase space coordinates. They can be thought of as a
/  property of the grid cell. 
/
/  Monte Carlo Tracers in each cell are stored as a linked list.
/
************************************************************************/

#ifndef MONTE_CARLO_TRACER_PARTICLE_DEFINED__
#define MONTE_CARLO_TRACER_PARTICLE_DEFINED__

#include "typedefs.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

class MonteCarloTracerParticle
{

 private:

  grid      *CurrentGrid;                  
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
  int        WillDelete;
  
  // history of particle position, including the time when the position was recorded. (NOT IMPLEMENTED)
  #ifdef TRACK_MC_HISTORY
  struct MonteCarloTracerParticleHistory
  {
      FLOAT Position[3];
      float Timestamp;
      MonteCarloTracerParticleHistory *NextFrame;
  };     
  MonteCarloTracerParticleHistory *LagrangianHistory;
  #endif 
  
  friend class grid;
  friend class ExternalBoundary;
  friend class MonteCarloTracerParticleList;
  
 public:

  MonteCarloTracerParticle *NextParticle;
  MonteCarloTracerParticle *PrevParticle; 
  
  // Constructors and destructor
  MonteCarloTracerParticle();
  MonteCarloTracerParticle(grid *_grid, int _ID, int _groupID, int _level, float _creationTime, FLOAT _pos[MAX_DIMENSION]);
  MonteCarloTracerParticle(grid *_grid, int _ID, int _groupID, int _level, float _creationTime, FLOAT _pos[MAX_DIMENSION], 
                            float _mass, int _exchange_count);
  ~MonteCarloTracerParticle();

  // Operators
  
  // Routines
  
  // getters
  int ReturnGridID(void){ return CurrentGrid->ID };
  int ReturnID(void){ return CurrentGrid->UniqueID };
  // setters
  
};

/* Comparer functions for sorting particle buffers with std::sort */

struct cmp_mc_grid {
  bool operator()(MonteCarloTracerParticle* const& a, MonteCarloTracerParticle* const& b) const {
    if (a->ReturnGridID() < b->ReturnGridID()) return true;
    else return false;
  }
};

struct cmp_mc_number {
  bool operator()(MonteCarloTracerParticle* const& a, MonteCarloTracerParticle* const& b) const {
    if (a->ReturnID() < b->ReturnID()) return true;
    else return false;
  }
};

class MonteCarloTracerParticleList
{
private:
  std::vector<MonteCarloTracerParticle*> internalBuffer;

public:
  MonteCarloTracerParticleList(void) {};
  MonteCarloTracerParticleList(const int inputParticleCount);
  MonteCarloTracerParticleList(const MonteCarloTracerParticleList &OtherList);
  ~MonteCarloTracerParticleList(void);
  MonteCarloTracerParticle*& operator[] (const int nIndex);
  MonteCarloTracerParticleList& operator=(const MonteCarloTracerParticleList& OtherList);
  void copy_and_insert(MonteCarloTracerParticle& input_particle);
  void insert(MonteCarloTracerParticle& input_particle);
  void clear(void);
  void erase(int index);
  void reserve(int size);
  void move_to_end(int index);
  void mark_for_deletion(int index);
  void delete_marked_particles(void);
  int size(void);
  void sort_grid(const int first, const int last);
  void sort_number(const int first, const int last);
};

// Constructors 
MonteCarloTracerParticleList::MonteCarloTracerParticleList(const int inputParticleCount)
{
  this->internalBuffer.reserve(inputParticleCount);
}

MonteCarloTracerParticleList::MonteCarloTracerParticleList(
    const MonteCarloTracerParticleList &OtherList)
{
  this->internalBuffer = OtherList.internalBuffer;
}

// Destructor
MonteCarloTracerParticleList::~MonteCarloTracerParticleList(void)
{
  for (typename std::vector<MonteCarloTracerParticle*>::iterator 
         it=this->internalBuffer.begin(); 
       it != this->internalBuffer.end(); ++it)
  {
    delete *it;
  }
  
}

MonteCarloTracerParticle*& MonteCarloTracerParticleList::operator[](const int nIndex)
{
  return this->internalBuffer.at(nIndex);
}

void MonteCarloTracerParticleList::copy_and_insert(MonteCarloTracerParticle& input_particle)
{
  this->internalBuffer.push_back(static_cast<MonteCarloTracerParticle*>(input_particle.clone()));
}

void MonteCarloTracerParticleList::insert(MonteCarloTracerParticle& input_particle)
{
  this->internalBuffer.push_back(static_cast<MonteCarloTracerParticle*>(&input_particle));
}

void MonteCarloTracerParticleList::clear(void)
{

  if (this->size() > 0) {
    for (typename std::vector<MonteCarloTracerParticle*>::iterator it=
           this->internalBuffer.begin(); it != this->internalBuffer.end(); ++it)
      {
        delete *it;
      }

    this->internalBuffer.clear();
  }
}

void MonteCarloTracerParticleList::erase(int index)
{
  delete this->internalBuffer[index];
  this->internalBuffer.erase(this->internalBuffer.begin() + index);
}

void MonteCarloTracerParticleList::reserve(int size)
{
  this->internalBuffer.reserve(size);
}

void MonteCarloTracerParticleList::move_to_end(int index)
{
  typename std::vector<MonteCarloTracerParticle*>::iterator it = 
    this->internalBuffer.begin() + index;
  std::rotate(it, it+1, this->internalBuffer.end());
}

void MonteCarloTracerParticleList::mark_for_deletion(int index)
{
  (*this)[index]->WillDelete = 1;
}

bool should_delete(MonteCarloTracerParticle* item) // deletes the MC particle but doesn't affect the list
{
  bool will_delete = item->ShouldDelete();
  if (will_delete) {
    delete item;
  }
  return will_delete;
} 

void MonteCarloTracerParticleList::delete_marked_particles(void)
{
  this->internalBuffer.erase(
      std::remove_if(
          this->internalBuffer.begin(),
          this->internalBuffer.end(),
          should_delete), 
      this->internalBuffer.end());
}

int MonteCarloTracerParticleList::size(void)
{
  return this->internalBuffer.size();
}

MonteCarloTracerParticleList& MonteCarloTracerParticleList::operator=(
    const MonteCarloTracerParticleList& OtherList)
{
  for (typename std::vector<MonteCarloTracerParticle*>::const_iterator it = 
         OtherList.internalBuffer.begin(); 
       it != OtherList.internalBuffer.end(); ++it)
  {
    this->internalBuffer.push_back((*it)->clone());
  }

  return *this;

}

void MonteCarloTracerParticleList::sort_grid(
    int first, int last)
{
  struct cmp_mc_grid comparator = cmp_mc_grid();
  if (this->size() > 0) {
    std::sort(
        this->internalBuffer.begin() + first, 
        this->internalBuffer.begin() + last,
        comparator);
  }
}

void MonteCarloTracerParticleList::sort_number(
    int first, int last)
{
  struct cmp_mc_number comparator = cmp_mc_number();
  if (this->size() > 0) {
    std::sort(
        this->internalBuffer.begin() + first, 
        this->internalBuffer.begin() + last,
        comparator);
  }
}

#endif