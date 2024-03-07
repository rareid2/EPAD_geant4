#include "G4SDManager.hh"
#include "MyRun.hh"
#include "G4HCofThisEvent.hh"
#include "Hit.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

MyRun::MyRun()
:hitsCollID1(0)
{ 
  hitsCollID1 = -1;
  // gran the hits collection names for each detector
  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  if (hitsCollID1 < 0) {

    G4String colNam;

    hitsCollID1 = SDman->GetCollectionID(colNam = "SD1/hitsCollection");
  }
}

MyRun::~MyRun() {}


void MyRun::RecordEvent(const G4Event* evt)
{
  if (hitsCollID1 < 0)
    return;
  G4HCofThisEvent *HCE = evt->GetHCofThisEvent();

  // otherwise proceed with getting the hits and saving to the hits file
  MyHitCollection *HC1 = 0;

  if (HCE)

  {

    HC1 = (MyHitCollection *)(HCE->GetHC(hitsCollID1));
  }

  if (HC1) {

    int n_hit = HC1->entries();

    for (int i = 0; i < n_hit; i++) {

      G4ThreeVector position = (*HC1)[i]->GetPosition();

      G4ThreeVector momentum = (*HC1)[i]->GetMomentum();

      G4double energy = (*HC1)[i]->GetEnergy();

      G4int particleID = (*HC1)[i]->GetID();

      G4String particlename = (*HC1)[i]->GetName();

      G4ThreeVector vertex = (*HC1)[i]->GetVertex();

      int thr=G4Threading::G4GetThreadId();
      std::stringstream ss;
      ss << thr;
      std::string numberAsString = ss.str();
      std::string myhitsfilename = "../simulation-data/hits_" + numberAsString + ".csv";
      
      std::ofstream hitFile;
      hitFile.open(myhitsfilename, std::ios_base::app);
      hitFile << "\n"
              << 1 << "," << position.x() / cm << "," << position.y() / cm
              << "," << position.z() / cm << "," << energy / keV << "," << particleID << "," << particlename
              << "," << vertex.x() / cm << "," << vertex.y() / cm << "," << vertex.z() / cm;
      hitFile.close();
    }
  }
  G4Run::RecordEvent(evt);
}

void MyRun::Merge(const G4Run* aRun)
{
  const MyRun* localRun = static_cast<const MyRun*>(aRun);
  G4Run::Merge(aRun);
} 