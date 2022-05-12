

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
  G4UIcmdWithADouble *fDcmd;
};

#endif
