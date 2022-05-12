

#include "PrimaryGeneratorMessenger.hh"

#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIdirectory.hh"
#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
    PrimaryGeneratorAction *prim)
    : G4UImessenger(), fPrimaryGenerator(prim) {
  fPrimDir = new G4UIdirectory("/particleSource/");

  // Sets fBackgroundDistribution
  fcmd = new G4UIcmdWithAnInteger("/particleSource/setPitchAngleDistribution",
                                  this);
  fcmd->SetParameterName("pitch angle distribution {0,1,2,3} "
                         "(90,sin,sin^2,tri).",
                         true);
  fcmd->SetDefaultValue(0);
  fcmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  // Sets E_folding
  fDcmd = new G4UIcmdWithADouble("/particleSource/setFoldingEnergy", this);
  fDcmd->SetParameterName("Folding Energy keV", true);
  fDcmd->SetDefaultValue(100.);
  fDcmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger() {
  delete fPrimDir;
  delete fcmd;
  delete fDcmd;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand *command,
                                            G4String newValue) {

  if (command == fcmd) {
    fPrimaryGenerator->SetPitchAngleDistribution(std::stoi(newValue));
  }

  if (command == fDcmd) {
    fPrimaryGenerator->SetFoldingEnergy(std::stod(newValue));
  }
}
