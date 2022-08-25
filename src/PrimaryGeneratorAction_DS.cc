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
// $Id: PrimaryGeneratorAction.cc 94307 2015-11-11 13:42:46Z gcosmo $
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction_DS.hh"

#include "PrimaryGeneratorMessenger.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include <fstream>
#include <stdexcept>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction(), fParticleGun(0), fE_folding(100.),
      fPI(3.14159265358979323846), fSphereR(15. * cm), fPitchAngleDist(0),
      fElectronParticle(0), fPrimaryGeneratorMessenger(0) {

  fParticleGun = new G4ParticleGun();

  fPrimaryGeneratorMessenger = new PrimaryGeneratorMessenger(this);

  fElectronParticle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fParticleGun;
  delete fPrimaryGeneratorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateElectrons(ParticleSample *r) {

  G4double pitchAngle, gyroPhase, u1, u2;

  // Switch-case for the background spatial distribution type
  switch (fPitchAngleDist) {

  case (0): // 90 deg (all perpendicular) distribution
  {

    // generate two random uniformly dist vars -1 to 1
    u1 = G4UniformRand() * 2. - 1.;
    u2 = G4UniformRand() * 2. - 1.;

    // sample spatially using rejection method
    // enter while loop if not less than 1
    while ((u1 * u1) + (u2 * u2) >= 1) {
      u1 = G4UniformRand() * 2. - 1.;
      u2 = G4UniformRand() * 2. - 1.;
    }

    // exit loop, condition met
    r->z = fSphereR * 2 * u1 * std::sqrt(1 - (u1 * u1) - (u2 * u2));
    r->x = fSphereR * 2 * u2 * std::sqrt(1 - (u1 * u1) - (u2 * u2));
    r->y = fSphereR * (1 - 2 * ((u1 * u1) + (u2 * u2)));

    // PITCH ANGLE DIST

    // pitchAngle ~ all perpendicular
    pitchAngle = fPI / 2;

    // set unifrom gyrophase and direction
    gyroPhase = G4UniformRand() * 2. * fPI;

    r->zDir = -std::sin(pitchAngle) * std::cos(gyroPhase);
    r->xDir = -std::sin(pitchAngle) * std::sin(gyroPhase);
    r->yDir = -std::cos(pitchAngle);

    break;
  }

  case (1): // sin(alpha) distribution
  {

    // generate two random uniformly dist vars -1 to 1
    u1 = G4UniformRand() * 2. - 1.;
    u2 = G4UniformRand() * 2. - 1.;

    // sample spatially using rejection method
    // enter while loop if not less than 1
    while ((u1 * u1) + (u2 * u2) >= 1) {
      u1 = G4UniformRand() * 2. - 1.;
      u2 = G4UniformRand() * 2. - 1.;
    }

    // exit loop, condition met
    r->z = fSphereR * 2 * u1 * std::sqrt(1 - (u1 * u1) - (u2 * u2));
    r->x = fSphereR * 2 * u2 * std::sqrt(1 - (u1 * u1) - (u2 * u2));
    r->y = fSphereR * (1 - 2 * ((u1 * u1) + (u2 * u2)));

    // PITCH ANGLE DIST
    // pitchAngle ~ sine -- inversion method
    pitchAngle = std::acos(G4UniformRand() * 2. - 1.);

    // set unifrom gyrophase and direction
    gyroPhase = G4UniformRand() * 2. * fPI;

    r->zDir = -std::sin(pitchAngle) * std::cos(gyroPhase);
    r->xDir = -std::sin(pitchAngle) * std::sin(gyroPhase);
    r->yDir = -std::cos(pitchAngle);

    break;
  }

  case (2): // sin^2(alpha) distribution
  {
    G4double u1_sinsq, u2_sinsq;

    // generate two random uniformly dist vars -1 to 1
    u1 = G4UniformRand() * 2. - 1.;
    u2 = G4UniformRand() * 2. - 1.;

    // sample spatially using rejection method
    // enter while loop if not less than 1
    while ((u1 * u1) + (u2 * u2) >= 1) {
      u1 = G4UniformRand() * 2. - 1.;
      u2 = G4UniformRand() * 2. - 1.;
    }

    // exit loop, condition met
    r->z = fSphereR * 2 * u1 * std::sqrt(1 - (u1 * u1) - (u2 * u2));
    r->x = fSphereR * 2 * u2 * std::sqrt(1 - (u1 * u1) - (u2 * u2));
    r->y = fSphereR * (1 - 2 * ((u1 * u1) + (u2 * u2)));

    // PITCH ANGLE DIST
    //  pitchAngle ~ sine^2 -- rejection sampling
    u1_sinsq = G4UniformRand();
    u2_sinsq = G4UniformRand() * fPI;

    while (u1_sinsq < std::sin(u2_sinsq) * std::sin(u2_sinsq)) {
      u1_sinsq = G4UniformRand();
      u2_sinsq = G4UniformRand() * fPI;
    }

    pitchAngle = u2_sinsq;

    // set unifrom gyrophase and direction
    gyroPhase = G4UniformRand() * 2. * fPI;

    r->zDir = -std::sin(pitchAngle) * std::cos(gyroPhase);
    r->xDir = -std::sin(pitchAngle) * std::sin(gyroPhase);
    r->yDir = -std::cos(pitchAngle);

    break;
  }

  case (3): // triangle distribution
  {
    G4double fc, bound1, bound2, tri_center, u_tri, tri_s;

    // generate two random uniformly dist vars -1 to 1
    u1 = G4UniformRand() * 2. - 1.;
    u2 = G4UniformRand() * 2. - 1.;

    // sample spatially using rejection method
    // enter while loop if not less than 1
    while ((u1 * u1) + (u2 * u2) >= 1) {
      u1 = G4UniformRand() * 2. - 1.;
      u2 = G4UniformRand() * 2. - 1.;
    }

    // exit loop, condition met
    r->z = fSphereR * 2 * u1 * std::sqrt(1 - (u1 * u1) - (u2 * u2));
    r->x = fSphereR * 2 * u2 * std::sqrt(1 - (u1 * u1) - (u2 * u2));
    r->y = fSphereR * (1 - 2 * ((u1 * u1) + (u2 * u2)));

    // pitch angle ~ triangle -- inversion method
    // in degrees
    fc = 0.5;
    bound1 = 75;
    bound2 = 105;
    tri_center = 90;
    u_tri = G4UniformRand();

    if (u_tri < fc) {
      tri_s =
          bound1 + std::sqrt(u_tri * (bound2 - bound1) * (tri_center - bound1));
    } else if (u_tri >= fc) {
      tri_s = bound2 - std::sqrt((1 - u_tri) * (bound2 - bound1) *
                                 (tri_center - bound1));
    } else {
      tri_s = 0.;
    }
    // convret to radians
    pitchAngle = tri_s * fPI / 180;

    // set unifrom gyrophase and direction
    gyroPhase = G4UniformRand() * 2. * fPI;

    r->zDir = -std::sin(pitchAngle) * std::cos(gyroPhase);
    r->xDir = -std::sin(pitchAngle) * std::sin(gyroPhase);
    r->yDir = -std::cos(pitchAngle);

    break;
  }

  default:
    throw std::invalid_argument("Enter a pitch angle distribution!");
  }

  // exponential with folding energy fE_folding
  r->energy = -fE_folding * std::log(G4UniformRand()) * keV;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent) {
  // Method called to populate member variables with number of
  // particles to generate

  // Selects electron for particle type
  fParticleGun->SetParticleDefinition(fElectronParticle);

  // Constant sphere offsets
  G4double xShift = 0.;
  G4double yShift = 0. * cm;
  G4double zShift = 500. * cm;

  // Struct that holds position, momentum direction, and energy
  ParticleSample *r = new ParticleSample();
  // r->x ~ x position
  // r->y ~ y position
  // r->z ~ z position
  // r->xDir ~ x momentum direction
  // r->yDir ~ y momentum direction
  // r->zDir ~ z momentum direction
  // r->energy ~ particle energy

  GenerateElectrons(r);

  fParticleGun->SetParticlePosition(
      G4ThreeVector(r->x + xShift, r->y + yShift, r->z + zShift));
  fParticleGun->SetParticleMomentumDirection(
      G4ThreeVector(r->xDir, r->yDir, r->zDir));
  fParticleGun->SetParticleEnergy(r->energy);
  fParticleGun->GeneratePrimaryVertex(anEvent);

  // Free particle utility struct
  delete r;
}
