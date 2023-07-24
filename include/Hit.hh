#ifndef Hit_h
#define Hit_h 1

#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4THitsCollection.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4VHit.hh"

class G4AttDef;
class G4AttValue;

// create a hit class that inherits from the G4VHit base class
class MyHit : public G4VHit {
public:
  MyHit() {}
  virtual ~MyHit() {}

  // G4int operator==(const MyHit &right) const;
  //  for memory allocation
  inline void *operator new(size_t);
  inline void operator delete(void *aHit);

  void Print();

private:
  // attributes?
  G4ThreeVector position;
  G4ThreeVector momentum;
  G4double energy;
  G4int particleID;
  G4String particlename;

public:
  // all kinds of fun stuff here
  inline void SetPosition(G4ThreeVector pos) { position = pos; }

  inline G4ThreeVector GetPosition() { return position; }

  inline void SetMomentum(G4ThreeVector mom) { momentum = mom; }

  inline G4ThreeVector GetMomentum() { return momentum; }

  inline void SetEnergy(G4double enep) { energy = enep; }

  inline G4double GetEnergy() { return energy; }

  inline void SetID(G4int hid) { particleID = hid; }

  inline G4int GetID() { return particleID; }

  inline void SetName(G4String pname) { particlename = pname; }

  inline G4String GetName() { return particlename; }
};

// defining what the hit collection is - not a class
typedef G4THitsCollection<MyHit> MyHitCollection;

// memory things
extern G4ThreadLocal G4Allocator<MyHit> *MyHitAllocator;
inline void *MyHit::operator new(size_t) {
  if (!MyHitAllocator)
    MyHitAllocator = new G4Allocator<MyHit>;
  return (void *)MyHitAllocator->MallocSingle();
}
inline void MyHit::operator delete(void *aHit) {
  MyHitAllocator->FreeSingle((MyHit *)aHit);
}

#endif
