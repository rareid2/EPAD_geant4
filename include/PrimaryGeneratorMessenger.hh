

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;

class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

class PrimaryGeneratorMessenger : public G4UImessenger {
public:
  PrimaryGeneratorMessenger(PrimaryGeneratorAction *);
  virtual ~PrimaryGeneratorMessenger();

  virtual void SetNewValue(G4UIcommand *, G4String);

private:
  PrimaryGeneratorAction *fPrimaryGenerator;
  G4UIdirectory *fPrimDir;
  G4UIcmdWithAnInteger *fcmd;
  G4UIcmdWithAnInteger *fcmd2;
  G4UIcmdWithAnInteger *fcmd3;
  G4UIcmdWithADouble *fDcmd;
  G4UIcmdWithADouble *fD2cmd;
  G4UIcmdWithADouble *fD3cmd;
  G4UIcmdWithADouble *fD4cmd;
  G4UIcmdWithADouble *fD5cmd;
  G4UIcmdWithADouble *fD6cmd;
  G4UIcmdWithAString *fScmd1;
};

#endif
