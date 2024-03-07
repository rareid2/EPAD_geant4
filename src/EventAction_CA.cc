//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: EventAction.cc 93886 2015-11-03 08:28:26Z gcosmo $
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "Hit.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction *runAction)
    : G4UserEventAction(), fRunAction(runAction) {

  // initialize ids to -1 to reconginze the start of the run
  hitsCollID1 = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event *event) {

  // gran the hits collection names for each detector
  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  if (hitsCollID1 < 0) {

    G4String colNam;

    hitsCollID1 = SDman->GetCollectionID(colNam = "SD1/hitsCollection");
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event *event) {
  // if -1 need to iniiailze
  if (hitsCollID1 < 0)
    return;
  G4HCofThisEvent *HCE = event->GetHCofThisEvent();

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

      // G4cout << "---- Hit # " << i << G4endl;

      // G4cout << " Position " << position / cm << " [cm] " << G4endl;

      // G4cout << " Momentum " << momentum / keV << " [keV] " << G4endl;

      // G4cout << " Energy   " << energy / keV << " [keV] " << G4endl;

      std::ofstream hitFile;
      hitFile.open("../simulation-data/hits.csv", std::ios_base::app);
      hitFile << "\n"
              << 1 << "," << position.x() / cm << "," << position.y() / cm
              << "," << position.z() / cm << "," << energy / keV << "," << particleID << "," << particlename
              << "," << vertex.x() / cm << "," << vertex.y() / cm << "," << vertex.z() / cm;
      hitFile.close();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
