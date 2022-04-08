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
// $Id: PrimaryGeneratorAction.hh 90623 2015-06-05 09:24:30Z gcosmo $
//
/// \file PrimaryGeneratorAction.hh
/// \brief Definition of the PrimaryGeneratorAction class

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
//#include "G4GeneralParticleSource.hh"
#include "globals.hh"

class G4ParticleGun;
class PrimaryGeneratorMessenger;
class G4Event;
class G4Box;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued 
/// in front of the phantom across 80% of the (X,Y) phantom size.

struct ParticleSample{
G4double x, y, z;
G4double xDir, yDir, zDir;
G4double energy;
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event* anEvent);         

    // Statistical selections of part. energy, position, and direction
    void GenerateLossConeElectrons(ParticleSample*);
    void GenerateTrappedElectrons(ParticleSample*);
    void GenerateSignalSource(ParticleSample*);
    void GenerateOtherDistributions(ParticleSample*);
    void GenerateSignalFromFile(ParticleSample*);
    
    // Methods for messenger class
    void SetWhichParticle(G4int partSelection) {fWhichParticle = partSelection;};

    // Sets background and signal folding energies
    void SetFoldingEnergy(G4double E0) { fE_folding = E0; };
    
    // Under GenerateSignalPhotons, selects limiting zenith angle
    void SetEventAngle(G4double ang) { fPhotonPhiLimitDeg = ang; };
    
    // Sets spatial distribution type under GenerateOtherDistributions
    void SetSpatialDistribution(G4int type) {fSpatialSignalDist = type; };    
    
    // Sets background spatial distribution
    void SetBackgroundSpatialDistribution(G4int type){fBackgroundSpatialDist=type;};


    void SetPhotonFileName(G4String fname){fPhotonFilename = fname;};

    // Under SpatialDistributions case 3, sets Gaussian parameters
    void SetThetaDirection(G4double ang) { fDirectionTheta = ang; }; 
    void SetThetaSigma(G4double stddev) { fThetaSigma = stddev; }; 
    void SetPhiDirection(G4double ang) { fDirectionPhi = ang; }; 
    void SetPhiSigma(G4double stddev) { fPhiSigma = stddev; }; 
    
    // Method to access particle gun
    const G4ParticleGun* GetParticleGun() const { return fParticleGun; }

  private:
    G4ParticleGun*  fParticleGun; // pointer a to G4 gun class
    
    G4double fE_folding;
    G4double fPI;
    G4double fDeg2Rad;
    G4double fSphereR;
    G4double fLossConeAngleDeg;
    G4double fPhotonPhiLimitDeg;
    
    G4double fDirectionTheta;
    G4double fThetaSigma;
    G4double fDirectionPhi;
    G4double fPhiSigma;

    G4int    fSpatialSignalDist;
    G4int    fBackgroundEnergyDist;    
    G4int    fBackgroundSpatialDist;
    G4int    fWhichParticle;

    G4String fPhotonFilename;
    unsigned static int fPhotonFileLineCounter;
    G4ParticleDefinition* fElectronParticle; 
    G4ParticleDefinition* fPhotonParticle; 
    PrimaryGeneratorMessenger* fPrimaryGeneratorMessenger;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
