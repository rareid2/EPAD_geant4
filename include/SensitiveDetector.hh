#ifndef SensitiveDetector_h
#define SensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "Hit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

// create a SD class that inherits from the SD base class in geant
class MySensitiveDetector : public G4VSensitiveDetector {
public:
  MySensitiveDetector(G4String SDname);
  virtual ~MySensitiveDetector();

  virtual void Initialize(G4HCofThisEvent *HCE);
  virtual G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROhist);
  // virtual void EndOfEvent(G4HCofThisEvent *HCE);
private:
  // assign hits collection to it
  MyHitCollection *hitsCollection;
  G4int collectionID;
};

#endif