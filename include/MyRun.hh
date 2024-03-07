#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "globals.hh"

class MyRun : public G4Run
{
  public:
    MyRun();
    virtual ~MyRun();
    virtual void RecordEvent(const G4Event*);
    virtual void Merge(const G4Run*);
  private:
    G4int hitsCollID1;
};
