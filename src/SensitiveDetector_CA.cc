#include "SensitiveDetector.hh"
#include "Hit.hh"

#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4ios.hh"

// constructor
MySensitiveDetector::MySensitiveDetector(G4String SDname)
    : G4VSensitiveDetector(SDname), hitsCollection(0), collectionID(-1) {
  G4String HCname;
  collectionName.insert(HCname = "hitsCollection");
}

MySensitiveDetector::~MySensitiveDetector() {}

// initialize - just get the collection ID once
void MySensitiveDetector::Initialize(G4HCofThisEvent *HCE) {
  // create a hit collection!
  hitsCollection =
      new MyHitCollection(SensitiveDetectorName, collectionName[0]);

  if (collectionID < 0) {
    collectionID =
        G4SDManager::GetSDMpointer()->GetCollectionID(hitsCollection);
  }
  HCE->AddHitsCollection(collectionID, hitsCollection);
}

G4bool MySensitiveDetector::ProcessHits(G4Step *aStep, G4TouchableHistory *) {
  G4StepPoint *preStep = aStep->GetPreStepPoint();
  G4TouchableHistory *touchable =
      (G4TouchableHistory *)(preStep->GetTouchable());

  // ensure only incoming tracks
  if (preStep->GetStepStatus() == fGeomBoundary) {

    MyHit *newHit = new MyHit();
    // gets the particles info from pre point
    newHit->SetPosition(aStep->GetPreStepPoint()->GetPosition());
    newHit->SetMomentum(aStep->GetPreStepPoint()->GetMomentum());
    // double check this is right
    newHit->SetEnergy(aStep->GetPreStepPoint()->GetKineticEnergy());

    G4Track *aTrack = aStep->GetTrack();
    newHit->SetID(aTrack->GetParentID());
    newHit->SetName(aTrack->GetParticleDefinition()->GetParticleName());

    hitsCollection->insert(newHit);

    aTrack->SetTrackStatus(fStopAndKill);
  }

  return true;
}