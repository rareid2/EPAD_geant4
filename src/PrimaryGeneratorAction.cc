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

#include "PrimaryGeneratorAction.hh"

#include "PrimaryGeneratorMessenger.hh"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

#include <fstream>
#include <stdexcept>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

unsigned int PrimaryGeneratorAction::fPhotonFileLineCounter = 1;

PrimaryGeneratorAction::PrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fParticleGun(0),
  fE_folding(100.),
  fPI(3.14159265358979323846),
  fDeg2Rad(3.14159265358979323846 / 180.),
  // TODO: CHECK THE RADIUS
  fSphereR(15.*cm),
  fLossConeAngleDeg(64.),
  fPhotonPhiLimitDeg(40.), 
  fDirectionTheta(0.), 
  fThetaSigma(0.),
  fDirectionPhi(0.), 
  fPhiSigma(0.),
  fSpatialSignalDist(3),
  fBackgroundEnergyDist(0),
  fBackgroundSpatialDist(2),
  fWhichParticle(0),
  fPhotonFilename(),
  //fPhotonFileLineCounter(0),
  fElectronParticle(0),
  fPhotonParticle(0),
  fPrimaryGeneratorMessenger(0)
{

  fParticleGun  = new G4ParticleGun();
 
  fPrimaryGeneratorMessenger = new PrimaryGeneratorMessenger(this);

  fElectronParticle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
  
  fPhotonParticle = G4ParticleTable::GetParticleTable()->FindParticle("gamma");

  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fPrimaryGeneratorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GenerateLossConeElectrons(ParticleSample* r)
{
  // not implemented
  
  // This method generates electrons that were in the loss cone, but have
  // backscattered off of the atmosphere and are directed anti-Earthward
  // and towards the satellite.


  // Fits:
  // f = a*exp(b*x)
  // In order to sample from this distribution, we take
  //
  // => X := b*ln(b/a*U) 
  //
  // we want X to be bounded on [0, phi_lossCone], so we'll choose a U such that
  // 
  //  U ~ U( U_min, Umax ) = U( a/b , a/b*exp(phi_lossCone/b) )
  //
  // we need to scale phi_lossCone due to numerical issues, so 64 deg -> 0.064 [unitless]
  //
  // therefore U = Rand(0,1) * (U_max - U_min) + U_min
  //
  //

}

void PrimaryGeneratorAction::GenerateSignalSource(ParticleSample* r)
{
  // NOT IMPLEMENTED
  // This method generates a photon signal from a energetic particle
  // precipitation event that occurs at approximately 50 - 100 km 
  // altitude. Variables are photon energy and the 
  // zenith angle distribution of the photons.

}

void PrimaryGeneratorAction::GenerateOtherDistributions(ParticleSample* r)
{
  // NOT IMPLEMENTED
  
}



void PrimaryGeneratorAction::GenerateTrappedElectrons(ParticleSample* r)
{
  
    G4double pitchAngle, gyroPhase, u1, u2;

    // Switch-case for the background spatial distribution type
    switch(fBackgroundSpatialDist)
    {
    
      case(0):{ // sin(alpha) distribution

        // generate two random uniformly dist vars -1 to 1
        u1 = G4UniformRand()*2.-1.;
        u2 = G4UniformRand()*2.-1.;
  
        // sample spatially using rejection method
        // enter while loop if not less than 1
        while ((u1*u1) + (u2*u2) >=1) {
          u1 = G4UniformRand()*2.-1.;
          u2 = G4UniformRand()*2.-1.;
          }
        
        // exit loop, condition met
        r->z = fSphereR * 2 * u1 * std::sqrt(1- (u1*u1) - (u2*u2)); 
        r->x = fSphereR * 2 * u2 * std::sqrt(1- (u1*u1) - (u2*u2));
        r->y = fSphereR * (1 - 2 * ((u1*u1) + (u2*u2))) ;

        // pitchAngle ~ sine
        // cos^-1(U)
        // gyroPhase ~ U[0 , 2 pi]
        pitchAngle = std::acos(G4UniformRand()*2.-1.);
        gyroPhase  = G4UniformRand() * 2. *fPI;

        r->zDir = -std::sin(pitchAngle) * std::cos(gyroPhase);	 
        r->xDir = -std::sin(pitchAngle) * std::sin(gyroPhase);  
        r->zDir = -std::cos(pitchAngle);
        

        break;
      }
    
      case(1): // sin(2 alpha) distribution
      // NOT IMPLEMENTED

	      
      case(2): 
      // NOT IMPLEMENTED

      default:
        throw std::invalid_argument("Enter a background spatial distribution!");
    }

    // Switch-case for the energy distribution type
    switch(fBackgroundEnergyDist)
    {
      case(0): // exponential with folding energy fE_folding
        r->energy = -fE_folding * std::log( G4UniformRand() )*keV;
        break;

      case(1): // monoenergetic with energy fE_folding
	r->energy = fE_folding * keV;
	break;
      
      default:
	throw std::invalid_argument("Enter a background energy distribution type!");
    
    }
}


void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  // Method called to populate member variables with number of 
  // particles to generate
  //CalculateParticlesToGenerate();

  // Selects electron for particle type
  fParticleGun->SetParticleDefinition(fElectronParticle);

  // Constant sphere offsets
  G4double xShift = 0.;
  G4double yShift = 0.*cm;
  G4double zShift = 500.*cm;

  // Struct that holds position, momentum direction, and energy
  ParticleSample* r = new ParticleSample();
  // r->x ~ x position
  // r->y ~ y position
  // r->z ~ z position
  // r->xDir ~ x momentum direction
  // r->yDir ~ y momentum direction
  // r->zDir ~ z momentum direction
  // r->energy ~ particle energy 

  switch(fWhichParticle){
    case(0): // Background electron, outside loss cone

    GenerateTrappedElectrons(r);
    
    fParticleGun->SetParticlePosition(
		    G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(
		    G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    break;

    case(1): // Background electron, inside loss cone

    GenerateLossConeElectrons(r);
    
    fParticleGun->SetParticlePosition(
		    G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(
		    G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    break;

    case(2): // Signal photon


    case(3): // Signal photon, other distribution types

    
    case(4): // Signal photon, other distribution types

    default: 
       throw std::invalid_argument("Need to chose particle type!");
  }

  // Free particle utility struct
  delete r;

}

