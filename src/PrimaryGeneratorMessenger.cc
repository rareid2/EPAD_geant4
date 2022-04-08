

#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIdirectory.hh"


PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* prim):G4UImessenger(),fPrimaryGenerator(prim) 
{
  fPrimDir = new G4UIdirectory("/particleSource/");
  fPrimDir->SetGuidance("Select trapped background, loss cone background, or photon signal source.");
  
  // Sets fWhichParticle 
  fcmd = new G4UIcmdWithAnInteger("/particleSource/setBackgroundType",this);
  fcmd->SetParameterName("Background Type {0,1,2,3} (Trapped, Backscattered Loss Cone, Signal Photons, Signal photons-other distributions)",true);
  fcmd->SetDefaultValue(0);
  fcmd->AvailableForStates(G4State_PreInit, G4State_Idle);
 
  // Sets fDistType
  fcmd2 = new G4UIcmdWithAnInteger("/particleSource/setSpatialDistribution",this);
  fcmd2->SetParameterName("Special Distribution type.",true);
  fcmd2->SetDefaultValue(3);
  fcmd2->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  // Sets fBackgroundDistribution
  fcmd3 = new G4UIcmdWithAnInteger("/particleSource/setBackgroundSpatialDistribution",this);
  fcmd3->SetParameterName("Trapped electron spatial distribution {0,1,2} (sine, sine 2 theta, isotropic).",true);
  fcmd3->SetDefaultValue(2);
  fcmd3->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Sets E_folding
  fDcmd = new G4UIcmdWithADouble("/particleSource/setFoldingEnergy",this);
  fDcmd->SetParameterName("Folding Energy {100, 200, 300} keV",true);
  fDcmd->SetDefaultValue(100.);
  fDcmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Sets photonPhiLimitDeg
  fD2cmd= new G4UIcmdWithADouble("/particleSource/setEventAngularSize",this);
  fD2cmd->SetParameterName("Enter an angle within (0, 45] degrees.",true);
  fD2cmd->SetDefaultValue(100.);
  fD2cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
 
  // Sets fDirectionTheta
  fD3cmd= new G4UIcmdWithADouble("/particleSource/setThetaDirection",this);
  fD3cmd->SetParameterName("Enter a direction angle for theta [deg].",true);
  fD3cmd->SetDefaultValue(0.);
  fD3cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
 
  // Sets fThetaSigma
  fD4cmd= new G4UIcmdWithADouble("/particleSource/setThetaSigma",this);
  fD4cmd->SetParameterName("Enter an angular std dev for theta [deg].",true);
  fD4cmd->SetDefaultValue(0.);
  fD4cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Sets fDirectionPhi
  fD5cmd= new G4UIcmdWithADouble("/particleSource/setPhiDirection",this);
  fD5cmd->SetParameterName("Enter a direction angle for phi [deg].",true);
  fD5cmd->SetDefaultValue(0.);
  fD5cmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  // Sets fPhiSigma
  fD6cmd= new G4UIcmdWithADouble("/particleSource/setPhiSigma",this);
  fD6cmd->SetParameterName("Enter an angular std dev for phi [deg].",true);
  fD6cmd->SetDefaultValue(0.);
  fD6cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Sets fPhotonFilename
  fScmd1=new G4UIcmdWithAString("/dataCollection/setPhotonFileName",this);
  fScmd1->SetParameterName("Enter file name to draw photons from.",true);
  fScmd1->SetDefaultValue("test1_photons.csv");
  fScmd1->AvailableForStates(G4State_PreInit, G4State_Idle);
  
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fPrimDir;
  delete fcmd;
  delete fcmd2;
  delete fcmd3;
  delete fDcmd;
  delete fD2cmd;
  delete fD3cmd;
  delete fD4cmd;
  delete fD5cmd;
  delete fD6cmd;
  delete fScmd1; 
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, 
					    G4String newValue)
{

  if(command == fcmd){
    fPrimaryGenerator->SetWhichParticle(std::stoi(newValue));
  }    	  
  
  if(command == fcmd2){
    fPrimaryGenerator->SetSpatialDistribution(std::stoi(newValue));
  }    	  
  
  if(command == fcmd3){
    fPrimaryGenerator->SetBackgroundSpatialDistribution(std::stoi(newValue));
  }    	  

  if(command == fDcmd){
    fPrimaryGenerator->SetFoldingEnergy(std::stod(newValue));
  }
  
  if(command == fD2cmd){
    fPrimaryGenerator->SetEventAngle(std::stod(newValue));
  }
  
  if(command == fD3cmd){
    fPrimaryGenerator->SetThetaDirection(std::stod(newValue));
  }
  
  if(command == fD4cmd){
    fPrimaryGenerator->SetThetaSigma(std::stod(newValue));
  }
  
  if(command == fD5cmd){
    fPrimaryGenerator->SetPhiDirection(std::stod(newValue));
  }
  
  if(command == fD6cmd){
    fPrimaryGenerator->SetPhiSigma(std::stod(newValue));
  }

  if(command == fScmd1){
    G4String pathAndName = "../analysis/photonFiles/";
    pathAndName += newValue;
    fPrimaryGenerator->SetPhotonFileName(pathAndName);
  }
}
