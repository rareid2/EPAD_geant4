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
  fSphereR(25.*cm),
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

  G4double phi_lossCone = 64./1000.; // scaled degrees, upper limit on returning electrons

  // Selects parameter set to be used per distribution
  G4double E0;
  G4double a100, b100;
  G4double a200, b200;
  G4double a300, b300;
  switch(fBackgroundSpatialDist)
  {
    case(0): // sine distribution
      a100 = 0.00705;
      b100 = 0.02356;
      a200 = 0.00759;
      b200 = 0.02296;
      a300 = 0.00978;
      b300 = 0.01496;
      break;
    case(1): // sine 2 distribution
      a100 = 0.00862;
      b100 = 0.02078;
      a200 = 0.00908;
      b200 = 0.02140;
      a300 = 0.01069;
      b300 = 0.01742;
      break;
    case(2): // isotropic distribution
      a100 = 0.00990;
      b100 = 0.01732;
      a200 = 0.01031;
      b200 = 0.01657;
      a300 = 0.01131;
      b300 = 0.01417;
      break;
    default:
	throw std::invalid_argument("Loss cone electron PAD parameter distribution selection is incorrect!");
  }

  // Selects parameters to be used with each energy range
  G4double a, b;

  // NB: E0's are for the sin distribution ONLY!!!!!

  if     (fE_folding == 100.){a = a100; b = b100; E0 = 85.;}
  else if(fE_folding == 200.){a = a200; b = b200; E0 = 147.;} 
  else if(fE_folding == 300.){a = a300; b = b300; E0 = 196.;} 
  else {throw std::invalid_argument("Loss cone electron PAD parameter selection energy range is incorrect!");}


  // Draw sample from exponential fit to pitch angle distribution
  G4double U_min = a / b;			       // U lower bound
  G4double U_max = a / b * std::exp(phi_lossCone / b); // U upper bound
  
  G4double U = G4UniformRand() * (U_max - U_min) + U_min; // shifts uniform random
  
  G4double lossConePhi = 1000*b*std::log(b/a * U); // 1000 due to scaling factor used above
 

  // Selects exponential folding energy E0 based on backscattered 
  // pitch angle range
  // TODO: write in energy-pitch angle coupling

  /*
  if(     lossConePhi < 20) {E0 = 159.;}
  else if(lossConePhi < 30) {E0 = 141.;}
  else if(lossConePhi < 40) {E0 = 177.;}
  else if(lossConePhi < 50) {E0 = 194.;}
  else if(lossConePhi < 64) {E0 = 230.;}
  */




  // NB: Mathematics spherical coordinates definition used below
  
  // Uniform sampling of angle about the field line theta on [0 , 2pi)
  G4double theta = G4UniformRand() * 2. * fPI; 
  
  // Zenith angle (pitch angle), converted to radians
  G4double phi = lossConePhi * fDeg2Rad;
 
  G4double phiDirection = phi + G4UniformRand() * 10. * fDeg2Rad;

  // We want our Y direction to be "up," otherwise standard
  // spherical to cartesian coordinate transform
  r->x = fSphereR * std::cos(theta) * std::sin(phi);
  r->y = fSphereR * std::cos(phi);
  r->z = fSphereR * std::sin(theta) * std::sin(phi);
  

  // Uniform random numbers on [0, 1) for particle direction
  G4double s1 = (G4UniformRand()*2.-1.)*fPI/8.;
  G4double s2 = (G4UniformRand()*2.-1.)*fPI/8.;
  G4double s3 = (G4UniformRand()*2.-1.)*fPI/8.;
  
  r->xDir = -std::sin(phiDirection+s1) * std::cos(theta);	   
  r->yDir = -std::cos(phiDirection+s2);
  r->zDir = -std::sin(phiDirection+s3) * std::sin(theta);


  // Enforces inward directionality to particles
  //if(r->x > 0) {r->xDir = -r->xDir;}
  //if(r->y > 0) {r->yDir = -r->yDir;}
  //if(r->z > 0) {r->zDir = -r->zDir;}


  // Continous inverse CDF sampling for exponential energy distribution
  // (only samples electrons with energy 50 keV or above)
  r->energy = -E0*std::log( G4UniformRand() )*keV;
 
  if(r->energy < 0.) {throw std::invalid_argument("E < 0 :(((((((");}
}

void PrimaryGeneratorAction::GenerateSignalFromFile(ParticleSample* r)
{
  std::ifstream photonFile;
  photonFile.open(fPhotonFilename, std::ios_base::in);

  G4String line, word;

  for(unsigned int i=0; i<fPhotonFileLineCounter; i++)
  {
    getline(photonFile, line);
  }

  std::stringstream s_ptr(line);
 
  getline(s_ptr, word, ',');
  r->x = std::stod(word) * cm;

  getline(s_ptr, word, ',');
  r->y = std::stod(word) * cm;
  
  getline(s_ptr, word, ',');
  r->z = std::stod(word) * cm;
  
  getline(s_ptr, word, ',');
  r->xDir = std::stod(word);
  
  getline(s_ptr, word, ',');
  r->yDir = std::stod(word);
  
  getline(s_ptr, word, '\n');
  r->zDir = std::stod(word);

  photonFile.close();

  fPhotonFileLineCounter++;
  
  // Power law sampling for photon energy
  // f(E) = C E^-k , where k > 0 and E > E_min; 
  //
  // Inverse CDF sampling:
  // E = E_min * (1 - U[0,1])^(-1/(k-1))
  r->energy = 50.*keV*std::pow((1 - G4UniformRand()),-1/(fE_folding-1));
}

void PrimaryGeneratorAction::GenerateSignalSource(ParticleSample* r)
{
  // This method generates a photon signal from a energetic particle
  // precipitation event that occurs at approximately 50 - 100 km 
  // altitude. Variables are photon energy and the 
  // zenith angle distribution of the photons.

  G4double theta;
  
  //Sample from exponential, energy dist. fitted to Wei's results
  // (valid for E0,source = 100 keV) 
  G4double shiftThreshold = 0.;    // keV
  G4double meanEnergy     = 241.4;  // keV

    
  r->energy = -(meanEnergy - shiftThreshold) *
	    (std::log( G4UniformRand() ) + shiftThreshold) * keV;
  
  

  // Uniformly distributed around azimuthal direction ab. fieldline 
  // theta ~ U[0, 2 pi]
  theta = G4UniformRand() * 2. * fPI;  
  
  // Phi (half angle) takes values in a cone determined by 
  // spacecraft altitude, precipitation event altitude and size
  G4double photonPhiLimitRad = fPhotonPhiLimitDeg * fDeg2Rad;       
  
  // u ~ U[0, 1-1/2*cos(phi_limit)] s.t. phi is distributed uniformly
  G4double u = G4UniformRand()*(1.-std::cos(photonPhiLimitRad))/2.;
  
  // phi ~ U[0, phi_limit] = cos^-1(1 - 2 u)
  G4double phi = std::acos(1 - 2 * u);


  // We want our Y direction to be "up"
  r->x = (fSphereR + 15.*cm) * std::sin(phi) * std::sin(theta);
  r->y = (fSphereR + 15.*cm) * std::cos(phi);
  r->z = (fSphereR + 15.*cm) * std::sin(phi) * std::cos(theta);


  // Inward unit normal
  r->xDir = -std::sin(theta) * std::sin(phi);
  r->zDir = -std::cos(theta) * std::sin(phi); 
  r->yDir = -std::cos(phi);


  // Uniform nudges s.t. each particle isn't aimed at the origin
  G4double nudgeFactor = 0.3;

  r->xDir += (G4UniformRand()*2. - 1.) * nudgeFactor; 
  r->yDir += (G4UniformRand()*2. - 1.) * nudgeFactor; 
  r->zDir += (G4UniformRand()*2. - 1.) * nudgeFactor; 
}

void PrimaryGeneratorAction::GenerateOtherDistributions(ParticleSample* r)
{
  
  G4double theta, phi, R;
  G4double narrowingOffset;
  G4double detectorSize;
  G4double N1, N2, U1, U2;
  
  G4double photonPhiLimitRad = fPhotonPhiLimitDeg * fDeg2Rad;

  detectorSize  = 250*2;  // mm , units assigned in switch-case

  G4double meanEnergyShift = 0.;
  r->energy = -( ( fE_folding-meanEnergyShift ) * 
		  std::log( G4UniformRand() ) + meanEnergyShift) * keV;
	
 
  G4double z_centers[4] = {-6.*cm, -2.*cm, 2.*cm, 6.*cm};
  G4double x_centers[4] = {-6.*cm, -2.*cm, 2.*cm, 6.*cm};
  G4int randIndex1, randIndex2;
  switch(fSpatialSignalDist)
  {
        case 0: // point source, near

		// Randomly selects integer from {0,1,2,3}
		// but not the sets {(2,2),(3,3),(2,3),(3,2)}
		// (that's where detectors aren't)
		do{
                  randIndex1 = std::floor(G4UniformRand()*4.);
                  randIndex2 = std::floor(G4UniformRand()*4.);
		} while((randIndex1 == 2 && randIndex2 == 2) ||
			(randIndex1 == 2 && randIndex2 == 3) ||
			(randIndex1 == 3 && randIndex2 == 2) ||
		        (randIndex1 == 3 && randIndex2 == 3));

		r->x = x_centers[randIndex1];
		r->z = z_centers[randIndex2];
                r->y = 20.*cm;

                narrowingOffset = 0.002;
                r->xDir = G4UniformRand()*narrowingOffset-narrowingOffset/2.;
                r->zDir = G4UniformRand()*narrowingOffset-narrowingOffset/2.;
                r->yDir = -1;
                break;

        case 1: // point source, infinitely far
		r->x = G4UniformRand()*detectorSize - detectorSize/2.; 
                r->x *= mm;
                r->z = G4UniformRand()*detectorSize - detectorSize/2.;
                r->z *= mm;
                r->y = 20.*cm;

                r->xDir = r->zDir = 0;
                r->yDir = -1;
                
		break;
        
	case 2: // structured circle
		
		do{
                  randIndex1 = std::floor(G4UniformRand()*4.);
                  randIndex2 = std::floor(G4UniformRand()*4.);
		} while((randIndex1 == 2 && randIndex2 == 2) ||
			(randIndex1 == 2 && randIndex2 == 3) ||
			(randIndex1 == 3 && randIndex2 == 2) ||
		        (randIndex1 == 3 && randIndex2 == 3));

		//Filled circle
                R = std::sqrt(G4UniformRand() * 10.) * cm;
                
		// Empty circle
		//R = 2. * cm;
		theta = G4UniformRand() * 2. * fPI;
                r->x = R * std::cos(theta) + x_centers[randIndex1] + 1.5*cm;
                r->z = R * std::sin(theta) + z_centers[randIndex2] + 1.5*cm;
                r->y = 20.*cm;
                
                narrowingOffset = 0.2;
		r->xDir = (G4UniformRand()*2.-1.)*narrowingOffset;
                r->zDir = (G4UniformRand()*2.-1.)*narrowingOffset;
 		r->yDir = -1;

                break;

        case 3: // Gaussian spot (default spatial distribution)
	
		// REMEMBER this is two rotations, zenith angle away from
		// the Y-axis, THEN about Y with theta

		// Box-Mueller transform for standard normals, N1 & N2
		U1 = G4UniformRand();
		U2 = G4UniformRand();
		N1 = std::sqrt(-2*std::log(U1)) * std::cos(2*fPI*U2);
		N2 = std::sqrt(-2*std::log(U1)) * std::sin(2*fPI*U2);

		// X ~ N(x_mean, sigma^2)
		theta = (fDirectionTheta + N1 * fThetaSigma) * fDeg2Rad;
		phi   = (fDirectionPhi   + N2 * fPhiSigma)   * fDeg2Rad;

                // Define the random position on a plane above AXIS
		r->x = G4UniformRand()*detectorSize - detectorSize/2.;
                r->x *= mm;
                r->x -= 5*cm;

                r->z = G4UniformRand()*detectorSize - detectorSize/2.;
                r->z *= mm;
		
		// Height of plane
		r->y = 20.*cm;

                // Directionality of particle
		r->xDir = std::sin(phi) * std::cos(theta);
                r->zDir = std::sin(phi) * std::sin(theta);
                r->yDir = -std::cos(phi);
                break;

        case 4: // Uniformly distributed phi angle

                R = 20.*cm;
 		
		// Angle of particle about field line
                theta = G4UniformRand()*2.*fPI;

                // Zenith angle of particles wrt to detector
                phi = std::acos(1 - 2 * G4UniformRand()*
                  (1.-std::cos(photonPhiLimitRad))/2.);

                // Position on surface of sphere
                r->x = R * std::sin(phi) * std::cos(theta);
                r->y = R * std::sin(phi) * std::sin(theta);
                r->z = -R * std::cos(phi);

                // Particle direction 
                r->xDir = G4UniformRand()*2. - 1.; // |x| ~ U[-1, 1]
                r->yDir = G4UniformRand()*2. - 1.; // |y| ~ U[-1, 1]

                // |z| ~ U[0, phi_limit]
                r->zDir = std::sqrt(r->xDir * r->xDir + r->yDir * r->yDir)/
                std::tan(G4UniformRand()*photonPhiLimitRad);
                break;

	case 5: // distributed source (for real this time)

		phi   = (G4UniformRand() * 40.)* fDeg2Rad; // ~U[0,10] deg
		theta = G4UniformRand() * 2. * fPI;       // ~U[0,2pi]rad

		r->x = (G4UniformRand()*2.-1.)*10.*cm; // ~U[-5, 5]
		r->y = 20.*cm;
		r->z = (G4UniformRand()*2.-1.)*10.*cm; // ~U[-5, 5]

		r->xDir = std::sin(phi) * std::cos(theta);
                r->zDir = std::sin(phi) * std::sin(theta);
                r->yDir = -std::cos(phi);

		break;
        default:
                throw std::invalid_argument("Choose distribution type!");
  }

  }


void PrimaryGeneratorAction::GenerateTrappedElectrons(ParticleSample* r)
{
  
    // Loss cone angle (same as polar angle, phi) at 500 km, in radians
    G4double theta_exclusion = 64.*fPI/180.;
    G4double maxPitchAngle   = (90.-64.)*fPI/180.; // 26 deg 
    G4double pitchAngle, gyroPhase, sign;
    G4double heightAdj; 

    // Switch-case for the background spatial distribution type
    switch(fBackgroundSpatialDist)
    {
    
      case(0): // sin(alpha) distribution
   
        // pitchAngle ~ sine[0, maxPitchAngle]

	// sign is +/- 1 with 50/50 probability
	sign = G4UniformRand();
	if(sign > 0.5) sign = 1;
	else sign = -1;

	// U ~ Unif[-1, 1]
	// 90 deg - cos^-1(U) /(pi / angular distance of distribution) 
        pitchAngle = 90.*fDeg2Rad - sign * std::acos(G4UniformRand()*2.-1.)/(fPI/maxPitchAngle);
      
        // gyroPhase ~ U[0 , 2 pi]
        gyroPhase  = G4UniformRand() * 2. *fPI;
	
	heightAdj = (G4UniformRand() * 2. - 1.) * 10.;
	//heightAdj = 0.;

        r->x = fSphereR * std::cos(gyroPhase); 
        r->y = heightAdj*cm + fSphereR * std::cos(pitchAngle);
        r->z = fSphereR * std::sin(gyroPhase);

        r->xDir = -std::sin(pitchAngle) * std::cos(gyroPhase);	   
        r->yDir = -std::cos(pitchAngle);
        r->zDir = -std::sin(pitchAngle) * std::sin(gyroPhase);

        break;
    
      case(1): // sin(2 alpha) distribution
      
        // pitchAngle ~ sine[0, maxPitchAngle / 2]

	// sign is +/- 1 with 50/50 probability
	sign = G4UniformRand();
	if(sign > 0.5) sign = 1;
	else sign = -1;

	// U ~ Unif[-1, 1]
	// 90 deg - cos^-1(U) /(pi / angular distance of distribution) 
        pitchAngle = 90.*fDeg2Rad - sign * std::acos(G4UniformRand()*2.-1.)/(fPI/(2.*maxPitchAngle));
      
        // gyroPhase ~ U[0 , 2 pi]
        gyroPhase  = G4UniformRand() * 2. *fPI;
	   
        r->x = fSphereR * std::cos(gyroPhase); 
        r->y = 2.5*cm + fSphereR * std::cos(pitchAngle);
        r->z = fSphereR * std::sin(gyroPhase);

        r->xDir = -std::sin(pitchAngle) * std::cos(gyroPhase);	   
        r->yDir = -std::cos(pitchAngle);
        r->zDir = -std::sin(pitchAngle) * std::sin(gyroPhase);
        break;
	      
      case(2): // isotropically distributed electrons outside of the 
	       // upward-going (anti-Earthward) loss cone

        // Calculate random particle position on sphere via rejection 
        // sampling, excluding spherical cap that makes up the loss cone
        do{
          
	  // Rand on [-1, 1)
          G4double u = G4UniformRand()*2.-1.;

          // Rand on [0, 2*pi)
          G4double theta = G4UniformRand()*2.*fPI;
          r->x = fSphereR * std::sqrt(1 - u * u) * std::cos(theta); 
          r->y = fSphereR * u; 
          r->z = fSphereR * std::sqrt(1 - u * u) * std::sin(theta); 
          }
        // exits when y position falls below spherical cap
        while(r->y > fSphereR * std::cos(theta_exclusion));

        // Uniform random numbers on [0, 1)
        r->xDir = G4UniformRand();
        r->yDir = G4UniformRand();
        r->zDir = G4UniformRand();

        // Enforces inward directionality to particles
        if(r->x > 0) {r->xDir = -r->xDir;}
        if(r->y > 0) {r->yDir = -r->yDir;}
        if(r->z > 0) {r->zDir = -r->zDir;}

        break;

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
  G4double yShift = 2.*cm;
  G4double zShift = 0.;

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

    // Selects photon for particle type
    fParticleGun->SetParticleDefinition(fPhotonParticle);
  
    GenerateSignalSource(r);
    
    fParticleGun->SetParticlePosition(
		    G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(
		    G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    break;

    case(3): // Signal photon, other distribution types

    // Selects photon for particle type
    fParticleGun->SetParticleDefinition(fPhotonParticle);
  
    GenerateOtherDistributions(r);
    
    fParticleGun->SetParticlePosition(
		    G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(
		    G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    break;
    
    case(4): // Signal photon, other distribution types

    // Selects photon for particle type
    fParticleGun->SetParticleDefinition(fPhotonParticle);
  
    GenerateSignalFromFile(r);
    
    fParticleGun->SetParticlePosition(
		    G4ThreeVector(r->x+xShift, r->y+yShift, r->z+zShift));
    fParticleGun->SetParticleMomentumDirection(
		    G4ThreeVector(r->xDir, r->yDir, r->zDir));
    fParticleGun->SetParticleEnergy(r->energy);
    fParticleGun->GeneratePrimaryVertex(anEvent);
    break;
    
    default: 
       throw std::invalid_argument("Need to chose particle type!");
  }

  // Free particle utility struct
  delete r;

}

